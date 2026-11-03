// gpu_probe.cu — cross-platform compute-capability probe (incl. 5090D AI-TOPS lock detection)
// Measures: FP32 FMA peak, tensor-core dense peaks (FP16/BF16/INT8/FP8), b1 BMMA rate, HBM bandwidth.
// Usage: gpu_probe.exe [out.tsv]   (default: probe.tsv in cwd)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cuda_runtime.h>
#include <mma.h>

#define CK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return 1; } } while (0)

// ---------------- FP32 FMA peak ----------------
__global__ void fp32_kernel(float* out, int iters) {
    float a = 1.000001f, b = 0.999999f, c = 0.000001f;
    float d = 1.000002f, e = 0.999998f, f = 0.000002f;
    float g = 1.000003f, h = 0.999997f, i2 = 0.000003f;
    float j = 1.000004f, k = 0.999996f, l = 0.000004f;
    for (int t = 0; t < iters; t++) {
        a = fmaf(a, b, c); d = fmaf(d, e, f);
        g = fmaf(g, h, i2); j = fmaf(j, k, l);
        a = fmaf(a, b, c); d = fmaf(d, e, f);
        g = fmaf(g, h, i2); j = fmaf(j, k, l);
    }
    if (out) out[threadIdx.x + blockIdx.x * blockDim.x] = a + d + g + j;
}

// ---------------- FP16 tensor mma peak (m16n8k16) ----------------
__global__ void fp16mma_kernel(float* out, int iters) {
    unsigned a0 = threadIdx.x, a1 = threadIdx.x + 1, a2 = threadIdx.x + 2, a3 = threadIdx.x + 3;
    unsigned b0 = blockIdx.x, b1 = blockIdx.x + 1;
    float c0 = 0.f, c1 = 0.f, c2 = 0.f, c3 = 0.f;
    float d0 = 0.f, d1 = 0.f, d2 = 0.f, d3 = 0.f;
    float e0 = 0.f, e1 = 0.f, e2 = 0.f, e3 = 0.f;
    float f0 = 0.f, f1 = 0.f, f2 = 0.f, f3 = 0.f;
    for (int t = 0; t < iters; t++) {
#define MMA16(c0,c1,c2,c3) asm volatile( \
        "mma.sync.aligned.m16n8k16.row.col.f32.f16.f16.f32 " \
        "{%0,%1,%2,%3}, {%4,%5,%6,%7}, {%8,%9}, {%0,%1,%2,%3};" \
        : "+f"(c0), "+f"(c1), "+f"(c2), "+f"(c3) \
        : "r"(a0), "r"(a1), "r"(a2), "r"(a3), "r"(b0), "r"(b1))
        MMA16(c0,c1,c2,c3); MMA16(d0,d1,d2,d3);
        MMA16(e0,e1,e2,e3); MMA16(f0,f1,f2,f3);
        MMA16(c0,c1,c2,c3); MMA16(d0,d1,d2,d3);
        MMA16(e0,e1,e2,e3); MMA16(f0,f1,f2,f3);
    }
    if (out) out[threadIdx.x + blockIdx.x * blockDim.x] = c0+c1+c2+c3+d0+d1+d2+d3+e0+e1+e2+e3+f0+f1+f2+f3;
}

// ---------------- INT8 tensor mma peak (m16n8k32) ----------------
__global__ void int8mma_kernel(int* out, int iters) {
    unsigned a0 = threadIdx.x, a1 = threadIdx.x + 1, a2 = threadIdx.x + 2, a3 = threadIdx.x + 3;
    unsigned b0 = blockIdx.x, b1 = blockIdx.x + 1;
    int c0 = 0, c1 = 0, c2 = 0, c3 = 0;
    int d0 = 0, d1 = 0, d2 = 0, d3 = 0;
    int e0 = 0, e1 = 0, e2 = 0, e3 = 0;
    int f0 = 0, f1 = 0, f2 = 0, f3 = 0;
    for (int t = 0; t < iters; t++) {
#define MMA8(c0,c1,c2,c3) asm volatile( \
        "mma.sync.aligned.m16n8k32.row.col.s32.s8.s8.s32 " \
        "{%0,%1,%2,%3}, {%4,%5,%6,%7}, {%8,%9}, {%0,%1,%2,%3};" \
        : "+r"(c0), "+r"(c1), "+r"(c2), "+r"(c3) \
        : "r"(a0), "r"(a1), "r"(a2), "r"(a3), "r"(b0), "r"(b1))
        MMA8(c0,c1,c2,c3); MMA8(d0,d1,d2,d3);
        MMA8(e0,e1,e2,e3); MMA8(f0,f1,f2,f3);
        MMA8(c0,c1,c2,c3); MMA8(d0,d1,d2,d3);
        MMA8(e0,e1,e2,e3); MMA8(f0,f1,f2,f3);
    }
    if (out) out[threadIdx.x + blockIdx.x * blockDim.x] = c0+c1+c2+c3+d0+d1+d2+d3+e0+e1+e2+e3+f0+f1+f2+f3;
}

// ---------------- FP8 tensor mma peak (e4m3, m16n8k32; sm_89+) ----------------
__global__ void fp8mma_kernel(float* out, int iters) {
#if __CUDA_ARCH__ >= 890
    unsigned a0 = threadIdx.x, a1 = threadIdx.x + 1, a2 = threadIdx.x + 2, a3 = threadIdx.x + 3;
    unsigned b0 = blockIdx.x, b1 = blockIdx.x + 1;
    float c0 = 0.f, c1 = 0.f, c2 = 0.f, c3 = 0.f;
    float d0 = 0.f, d1 = 0.f, d2 = 0.f, d3 = 0.f;
    float e0 = 0.f, e1 = 0.f, e2 = 0.f, e3 = 0.f;
    float f0 = 0.f, f1 = 0.f, f2 = 0.f, f3 = 0.f;
    for (int t = 0; t < iters; t++) {
#define MMAF8(c0,c1,c2,c3) asm volatile( \
        "mma.sync.aligned.m16n8k32.row.col.f32.e4m3.e4m3.f32 " \
        "{%0,%1,%2,%3}, {%4,%5,%6,%7}, {%8,%9}, {%0,%1,%2,%3};" \
        : "+f"(c0), "+f"(c1), "+f"(c2), "+f"(c3) \
        : "r"(a0), "r"(a1), "r"(a2), "r"(a3), "r"(b0), "r"(b1))
        MMAF8(c0,c1,c2,c3); MMAF8(d0,d1,d2,d3);
        MMAF8(e0,e1,e2,e3); MMAF8(f0,f1,f2,f3);
        MMAF8(c0,c1,c2,c3); MMAF8(d0,d1,d2,d3);
        MMAF8(e0,e1,e2,e3); MMAF8(f0,f1,f2,f3);
    }
    if (out) out[threadIdx.x + blockIdx.x * blockDim.x] = c0+c1+c2+c3+d0+d1+d2+d3+e0+e1+e2+e3+f0+f1+f2+f3;
#else
    if (out) out[threadIdx.x + blockIdx.x * blockDim.x] = 0.f;
#endif
}

// ---------------- b1 BMMA rate (same wmma path as the paper's router) ----------------
__global__ void b1bmma_kernel(int* out, int iters) {
    namespace wx = nvcuda::wmma::experimental;
    nvcuda::wmma::fragment<nvcuda::wmma::accumulator, 8, 8, 128, int> acc0, acc1, acc2, acc3;
    nvcuda::wmma::fragment<nvcuda::wmma::matrix_a, 8, 8, 128, wx::precision::b1, nvcuda::wmma::row_major> af;
    nvcuda::wmma::fragment<nvcuda::wmma::matrix_b, 8, 8, 128, wx::precision::b1, nvcuda::wmma::col_major> bf;
    __shared__ __align__(32) int smem[128];
    for (int i = threadIdx.x; i < 128; i += blockDim.x) smem[i] = (int)(i * 2654435761u);
    __syncthreads();
    nvcuda::wmma::load_matrix_sync(af, reinterpret_cast<unsigned*>(&smem[0]), 4);
    nvcuda::wmma::load_matrix_sync(bf, reinterpret_cast<unsigned*>(&smem[32]), 4);
    nvcuda::wmma::fill_fragment(acc0, 0);
    nvcuda::wmma::fill_fragment(acc1, 0);
    nvcuda::wmma::fill_fragment(acc2, 0);
    nvcuda::wmma::fill_fragment(acc3, 0);
    for (int t = 0; t < iters; t++) {
#define BMMA(A) nvcuda::wmma::bmma_sync(A, af, bf, A, wx::bmmaBitOpAND, wx::bmmaAccumulateOpPOPC)
        BMMA(acc0); BMMA(acc1); BMMA(acc2); BMMA(acc3);
        BMMA(acc0); BMMA(acc1); BMMA(acc2); BMMA(acc3);
    }
    int r = 0;
    int* p;
    nvcuda::wmma::store_matrix_sync(&smem[64], acc0, 8, nvcuda::wmma::mem_row_major);
    nvcuda::wmma::store_matrix_sync(&smem[96], acc1, 8, nvcuda::wmma::mem_row_major);
    p = &smem[64]; r += p[threadIdx.x & 31] + p[(threadIdx.x + 32) & 63];
    p = &smem[96]; r += p[threadIdx.x & 31];
    r += acc2.x[0] + acc3.x[0];
    if (out) out[threadIdx.x + blockIdx.x * blockDim.x] = r;
}

// ---------------- HBM bandwidth (float4 saxpy) ----------------
__global__ void bw_kernel(float4* dst, const float4* src, size_t n, float a) {
    size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    size_t stride = (size_t)gridDim.x * blockDim.x;
    for (; i < n; i += stride) {
        float4 v = src[i];
        v.x = a * v.x + v.x; v.y = a * v.y + v.y; v.z = a * v.z + v.z; v.w = a * v.w + v.w;
        dst[i] = v;
    }
}

static float time_kernel(void (*launch)(int), int iters) { (void)launch; return 0.f; } // unused placeholder

struct Timer {
    cudaEvent_t a, b;
    Timer() { cudaEventCreate(&a); cudaEventCreate(&b); }
    ~Timer() { cudaEventDestroy(a); cudaEventDestroy(b); }
};

template <typename F>
static float bench(F&& f, int reps) {
    f(); // warmup
    CK(cudaDeviceSynchronize());
    float best = 1e30f;
    Timer t;
    for (int r = 0; r < reps; r++) {
        cudaEventRecord(t.a);
        f();
        cudaEventRecord(t.b);
        CK(cudaEventSynchronize(t.b));
        float ms; cudaEventElapsedTime(&ms, t.a, t.b);
        best = std::min(best, ms);
    }
    return best;
}

int main(int argc, char** argv) {
    const char* outp = (argc >= 2) ? argv[1] : "probe.tsv";
    FILE* f = fopen(outp, "w");
    if (!f) { fprintf(stderr, "cannot open %s\n", outp); return 1; }
    fprintf(f, "probe\tvalue\tunit\n");
    const char* GPUNAME = nullptr;

    cudaDeviceProp prop;
    int dev = 0;
    CK(cudaGetDevice(&dev));
    CK(cudaGetDeviceProperties(&prop, dev));
    int driver = 0; cudaDriverGetVersion(&driver);
    fprintf(f, "gpu_name\t%s\t\n", prop.name);
    GPUNAME = prop.name;
    fprintf(f, "compute_capability\t%d.%d\t\n", prop.major, prop.minor);
    fprintf(f, "sm_count\t%d\t\n", prop.multiProcessorCount);
    fprintf(f, "clock_rate_mhz\t%d\t\n", prop.clockRate / 1000);
    fprintf(f, "mem_clock_rate_mhz\t%d\t\n", prop.memoryClockRate / 1000);
    fprintf(f, "mem_bus_width_bit\t%d\t\n", prop.memoryBusWidth);
    fprintf(f, "mem_total_gb\t%.1f\t\n", prop.totalGlobalMem / 1073741824.0);
    fprintf(f, "driver_version\t%d\t\n", driver);
    fprintf(f, "l2_cache_mb\t%.1f\t\n", prop.l2CacheSize / 1048576.0);

    double v_fp32 = 0, v_fp16 = 0, v_int8 = 0, v_fp8 = -1, v_b1 = 0, v_bw = 0;
    const int blocks = prop.multiProcessorCount * 8;
    const int threads = 256;
    printf("gpu=%s cc=%d.%d sm=%d blocks=%d\n", prop.name, prop.major, prop.minor, prop.multiProcessorCount, blocks);

    { // FP32 peak
        int iters = 400000;
        float* d = nullptr;
        float ms = bench([&]{ fp32_kernel<<<blocks, threads>>>(d, iters); }, 3);
        double flops = 2.0 * 8.0 * iters * (double)blocks * threads;
        double tf = flops / (ms / 1000.0) / 1e12;
        v_fp32 = tf; fprintf(f, "fp32_peak_tflops\t%.2f\tTFLOPS\n", tf);
        printf("fp32 = %.2f TFLOPS\n", tf);
    }
    { // FP16 mma peak
        int iters = 300000;
        float* d = nullptr;
        float ms = bench([&]{ fp16mma_kernel<<<blocks, threads>>>(d, iters); }, 3);
        double mma_per_warp = (double)iters * 8.0;               // 8 mma per loop iter
        double warps = (double)blocks * threads / 32.0;
        double flops = mma_per_warp * warps * (16.0 * 8 * 16 * 2);
        double tf = flops / (ms / 1000.0) / 1e12;
        v_fp16 = tf; fprintf(f, "fp16_mma_tflops\t%.2f\tTFLOPS\n", tf);
        printf("fp16 mma = %.2f TFLOPS\n", tf);
    }
    { // INT8 mma peak
        int iters = 300000;
        int* d = nullptr;
        float ms = bench([&]{ int8mma_kernel<<<blocks, threads>>>(d, iters); }, 3);
        double mma_per_warp = (double)iters * 8.0;
        double warps = (double)blocks * threads / 32.0;
        double tops = mma_per_warp * warps * (16.0 * 8 * 32 * 2) / (ms / 1000.0) / 1e12;
        v_int8 = tops; fprintf(f, "int8_mma_tops\t%.2f\tTOPS\n", tops);
        printf("int8 mma = %.2f TOPS\n", tops);
    }
    if (prop.major >= 8 && prop.minor >= 9 || prop.major >= 9) { // FP8 mma (sm_89+)
        int iters = 300000;
        float* d = nullptr;
        float ms = bench([&]{ fp8mma_kernel<<<blocks, threads>>>(d, iters); }, 3);
        double mma_per_warp = (double)iters * 8.0;
        double warps = (double)blocks * threads / 32.0;
        double flops = mma_per_warp * warps * (16.0 * 8 * 32 * 2);
        double tf = flops / (ms / 1000.0) / 1e12;
        v_fp8 = tf; fprintf(f, "fp8_mma_tflops\t%.2f\tTFLOPS\n", tf);
        printf("fp8 mma = %.2f TFLOPS\n", tf);
    } else {
        fprintf(f, "fp8_mma_tflops\tskip_cc<89\t\n");
    }
    { // b1 BMMA rate (paper's own engine path)
        int iters = 300000;
        int* d = nullptr;
        float ms = bench([&]{ b1bmma_kernel<<<blocks, threads>>>(d, iters); }, 3);
        double calls_per_warp = (double)iters * 8.0;
        double warps = (double)blocks * threads / 32.0;
        double bitops = calls_per_warp * warps * (8.0 * 8 * 128);
        double tops = bitops / (ms / 1000.0) / 1e12;
        v_b1 = tops; fprintf(f, "b1_bmma_tops\t%.2f\tTOPS(AND+POPC bit-ops)\n", tops);
        printf("b1 bmma = %.2f TOPS\n", tops);
    }
    { // HBM bandwidth
        size_t n = (size_t)128 * 1024 * 1024 / sizeof(float4);   // 128 MB per array
        float4 *src, *dst;
        CK(cudaMalloc(&src, n * sizeof(float4)));
        CK(cudaMalloc(&dst, n * sizeof(float4)));
        CK(cudaMemset(src, 1, n * sizeof(float4)));
        CK(cudaMemset(dst, 0, n * sizeof(float4)));
        int grid = prop.multiProcessorCount * 16;
        float ms = bench([&]{ bw_kernel<<<grid, 256>>>(dst, src, n, 1.5f); }, 5);
        double bytes = 2.0 * (double)n * sizeof(float4);
        double gbs = bytes / (ms / 1000.0) / 1e9;
        v_bw = gbs; fprintf(f, "hbm_bandwidth_gbs\t%.1f\tGB/s\n", gbs);
        printf("hbm = %.1f GB/s\n", gbs);
        cudaFree(src); cudaFree(dst);
    }

    // ---------------- reference table + lock verdict ----------------
    // Public spec sheet values (vendor marketing pages; FP32 & bandwidth are NOT
    // subject to AI-TOPS export limits, so they double as harness sanity checks).
    struct Ref { const char* key; double fp32, bw, note; };
    const Ref refs[] = {
        {"RTX 5090",   104.8, 1792.0},   // 5090D reports the same FP32/BW class silicon
        {"RTX 5080",    56.3,  960.0},
        {"RTX 4090",    82.6, 1008.0},
        {"RTX 4080",    48.7,  716.8},
        {"RTX 4070",    29.1,  504.2},
        {"RTX 3060 Ti", 16.2,  448.0},
        {"RTX 3060",    12.7,  360.0},
    };
    const Ref* R = nullptr;
    for (const Ref& r : refs) if (strstr(GPUNAME, r.key)) { R = &r; break; }
    if (!R) {
        fprintf(f, "spec_match\tnone\t(no reference entry for this GPU; raw numbers still valid)\n");
        fprintf(f, "lock_verdict\tINCONCLUSIVE-UNKNOWN-GPU\t(indicative only)\n");
    } else {
        double r_fp32 = v_fp32 / R->fp32, r_bw = v_bw / R->bw;
        double tensor_health = v_fp16 / (v_fp32 * 2.5);   // unlocked consumer ratio ~2.4-2.6 (f32-acc)
        double int8_health = v_int8 / (v_fp32 * 9.5);     // unlocked ratio ~9.6-10 (measured 3060 Ti)
        fprintf(f, "spec_ref\t%s fp32=%.1fTF bw=%.0fGB/s\t(public spec)\n", R->key, R->fp32, R->bw);
        fprintf(f, "ratio_fp32_of_spec\t%.3f\t\n", r_fp32);
        fprintf(f, "ratio_bw_of_spec\t%.3f\t\n", r_bw);
        fprintf(f, "tensor_health_fp16\t%.3f\t(unlocked ~1.0)\n", tensor_health);
        fprintf(f, "tensor_health_int8\t%.3f\t(unlocked ~1.0)\n", int8_health);
        const char* verdict;
        if (r_fp32 > 0.65 && r_bw > 0.75 && (tensor_health < 0.55 || int8_health < 0.55))
            verdict = "LOCK-SUSPECTED (FP32 & BW healthy but tensor path far below silicon-normal)";
        else if (r_fp32 > 0.65 && r_bw > 0.75 && tensor_health >= 0.8 && int8_health >= 0.8)
            verdict = "NO-LOCK-INDICATED (tensor paths at silicon-normal level)";
        else
            verdict = "INCONCLUSIVE (harness or spec ratio off; inspect raw numbers)";
        fprintf(f, "lock_verdict\t%s\t(indicative only, not authoritative)\n", verdict);
        printf("verdict: %s\n", verdict);
    }
    fclose(f);
    CK(cudaDeviceSynchronize());
    return 0;
}
