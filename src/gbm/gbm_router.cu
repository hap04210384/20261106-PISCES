// gbm_router.cu - Dual-engine measurement-driven router for exact pair counting.
// Two engines compute identical exact c11 for every family tile:
//   * BMMA engine : WMMA b1 (m8n8k128) tensor-core waves (as gbm_pair_gpu.cu)
//   * SIMT engine : thread-per-pair bitwise AND + popcount over 32-bit words
// A measurement-driven router assigns each 8x8 tile to an engine using a
// calibrated cost table (fixed setup + per-chunk rate per engine), with
// hysteresis, byte-budget batching (tile-list launches), and two-stream
// overlap. Because both engines are exact, routing affects ONLY performance -
// the emitted set and the sorted TSV are byte-identical for every routing
// decision (route-independence).
// Part of the reproducibility package for the associated study.
// Author: DAPR
//
// Modes:
//   gbm_router --count in.txt out.tsv [capacity] [flags]
//   gbm_router --bench N M seed reps [flags]     (in-memory, mt19937_64 density 0.1)
//   gbm_router --sweep in.txt [flags]            (dark-region per-chunk engine sweep)
// Flags:
//   --raw            raw mode (c11>=1, no screen)
//   --no-batch       batch_sz=1 tile per launch (default 8)
//   --no-overlap     serialize: device sync after every launch
//   --no-hysteresis  zero switching margin (default 15%)
//   --static-table   calibrate once, never refresh cost table
//   --static-route   pick single best engine at start, never re-decide (E-F04 baseline)
//   --no-exact       skip exactness protocol (sortedness audit)
//   --load P         fluctuating background load, P%% duty cycle (E-F04)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <thread>
using namespace std;

#include <cuda_runtime.h>
#include <mma.h>
using namespace nvcuda;

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(1);} } while (0)

// ---------------- GBM1 container (same on-disk format as gbm_count) ----------------
struct Matrix {
    int N = 0, M = 0;
    int words = 0;
    vector<uint64_t> bits;
    vector<string> famNames;
    vector<uint64_t> support;
};

static bool load_txt(Matrix& m, const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return false; }
    char magic[8];
    if (fscanf(f, "%7s", magic) != 1 || strcmp(magic, "GBM1") != 0) { fprintf(stderr, "bad magic\n"); fclose(f); return false; }
    if (fscanf(f, "%d %d", &m.N, &m.M) != 2) { fprintf(stderr, "bad header\n"); fclose(f); return false; }
    m.words = (m.N + 63) / 64;
    m.bits.assign((size_t)m.M * m.words, 0ULL);
    m.famNames.resize(m.M);
    m.support.assign(m.M, 0);
    char name[256];
    for (int j = 0; j < m.M; j++) {
        if (fscanf(f, "%255s", name) != 1) { fprintf(stderr, "truncated at family %d\n", j); fclose(f); return false; }
        m.famNames[j] = name;
        int c, i = 0;
        while ((c = fgetc(f)) == ' ' || c == '\t') {}
        if (c != '0' && c != '1' && c != '\n') { fprintf(stderr, "unexpected char %c at family %d\n", c, j); fclose(f); return false; }
        while (c == '0' || c == '1') {
            if (c == '1' && i < m.N) { m.bits[(size_t)j * m.words + (i >> 6)] |= (1ULL << (i & 63)); m.support[j]++; }
            i++;
            c = fgetc(f);
        }
        if (i != m.N) { fprintf(stderr, "family %d has %d bits, expected %d\n", j, i, m.N); fclose(f); return false; }
    }
    fclose(f);
    return true;
}

// ---------------- device side ----------------
struct Triple { int i, j, c11; };

// screen replicates gbm_count.cpp / gbm_pair_gpu.cu expression order (IEEE double)
__device__ __forceinline__ bool screen_keep(int c11, uint64_t sup_i, uint64_t sup_j, int N, int rawmode) {
    if (rawmode) return c11 >= 1;
    const double n  = (double)N;
    const double r1 = (double)sup_i, r2 = (double)sup_j;
    const double E = r1 * r2 / n;
    const double var = r1 * r2 * (n - r1) * (n - r2) / (n * n * (n - 1.0));
    if (var <= 0) return false;
    const double sd = sqrt(var);
    return fabs((double)c11 - E) >= 3.29 * sd;
}

// ---- list kernels: block q processes tile q of the batch list (batched launch) ----
__global__ void bmma_list_kernel(const uint32_t* __restrict__ bm,
                                 const uint64_t* __restrict__ sup,
                                 int M, int N, int Npad,
                                 const int2* __restrict__ tilelist, int nb,
                                 Triple* __restrict__ out,
                                 unsigned long long* __restrict__ counter,
                                 unsigned long long cap,
                                 int* __restrict__ overflow,
                                 int rawmode) {
    const int q = blockIdx.x;
    if (q >= nb) return;
    const int ti = tilelist[q].x;
    const int tj = tilelist[q].y;
    const int KCH = Npad / 128;
    const int ldm = Npad;          // bits
    const int roww = Npad / 32;    // uint32 words per family row

    __shared__ int smem[4][64];
    const int warp = threadIdx.x >> 5;

    wmma::fragment<wmma::accumulator, 8, 8, 128, int> c_frag;
    wmma::fill_fragment(c_frag, 0);
    for (int ch = warp; ch < KCH; ch += 4) {
        const uint32_t* ap = bm + (size_t)ti * 8 * roww + (size_t)ch * 4;
        const uint32_t* bp = bm + (size_t)tj * 8 * roww + (size_t)ch * 4;
        wmma::fragment<wmma::matrix_a, 8, 8, 128, wmma::experimental::precision::b1, wmma::row_major> a_frag;
        wmma::fragment<wmma::matrix_b, 8, 8, 128, wmma::experimental::precision::b1, wmma::col_major> b_frag;
        wmma::load_matrix_sync(a_frag, ap, ldm);
        wmma::load_matrix_sync(b_frag, bp, ldm);
        wmma::bmma_sync(c_frag, a_frag, b_frag, c_frag,
                        wmma::experimental::bmmaBitOpAND,
                        wmma::experimental::bmmaAccumulateOpPOPC);
    }
    wmma::store_matrix_sync(&smem[warp][0], c_frag, 8, wmma::mem_row_major);
    __syncthreads();

    if (threadIdx.x < 64) {
        const int ii = (int)(threadIdx.x >> 3);
        const int jj = (int)(threadIdx.x & 7);
        const int i = ti * 8 + ii;
        const int j = tj * 8 + jj;
        if (i >= j || i >= M || j >= M) return;
        int c11 = 0;
        #pragma unroll
        for (int w = 0; w < 4; w++) c11 += smem[w][ii * 8 + jj];
        if (screen_keep(c11, sup[i], sup[j], N, rawmode)) {
            unsigned long long idx = atomicAdd(counter, 1ULL);
            if (idx < cap) { out[idx].i = i; out[idx].j = j; out[idx].c11 = c11; }
            else *overflow = 1;
        }
    }
}

__global__ void simt_list_kernel(const uint32_t* __restrict__ bm,
                                 const uint64_t* __restrict__ sup,
                                 int M, int N, int Npad,
                                 const int2* __restrict__ tilelist, int nb,
                                 Triple* __restrict__ out,
                                 unsigned long long* __restrict__ counter,
                                 unsigned long long cap,
                                 int* __restrict__ overflow,
                                 int rawmode) {
    const int q = blockIdx.x;
    if (q >= nb) return;
    if (threadIdx.x >= 64) return;
    const int ti = tilelist[q].x;
    const int tj = tilelist[q].y;
    const int roww = Npad / 32;
    const int ii = (int)(threadIdx.x >> 3);
    const int jj = (int)(threadIdx.x & 7);
    const int i = ti * 8 + ii;
    const int j = tj * 8 + jj;
    if (i >= j || i >= M || j >= M) return;
    const uint32_t* a = bm + (size_t)i * roww;
    const uint32_t* b = bm + (size_t)j * roww;
    const int wn = (N + 31) / 32;
    int c11 = 0;
    for (int w = 0; w < wn; w++) c11 += __popc(a[w] & b[w]);
    if (screen_keep(c11, sup[i], sup[j], N, rawmode)) {
        unsigned long long idx = atomicAdd(counter, 1ULL);
        if (idx < cap) { out[idx].i = i; out[idx].j = j; out[idx].c11 = c11; }
        else *overflow = 1;
    }
}

// noise kernel for fluctuating-load experiments (E-F04)
__global__ void noise_kernel(float* buf, size_t n, unsigned long long spin) {
    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    float v = buf[idx];
    for (unsigned long long k = 0; k < spin; k++) v = v * 1.0000001f + 1e-9f;
    buf[idx] = v;
}

// full-grid engines (calibration / sweep / pure-mode reference)
__global__ void bmma_full_kernel(const uint32_t* __restrict__ bm,
                                 const uint64_t* __restrict__ sup,
                                 int M, int N, int Npad,
                                 Triple* __restrict__ out,
                                 unsigned long long* __restrict__ counter,
                                 unsigned long long cap,
                                 int* __restrict__ overflow,
                                 int rawmode) {
    const int ti = blockIdx.x;
    const int tj = blockIdx.y;
    if (tj < ti) return;
    const int KCH = Npad / 128;
    const int ldm = Npad;
    const int roww = Npad / 32;
    __shared__ int smem[4][64];
    const int warp = threadIdx.x >> 5;
    wmma::fragment<wmma::accumulator, 8, 8, 128, int> c_frag;
    wmma::fill_fragment(c_frag, 0);
    for (int ch = warp; ch < KCH; ch += 4) {
        const uint32_t* ap = bm + (size_t)ti * 8 * roww + (size_t)ch * 4;
        const uint32_t* bp = bm + (size_t)tj * 8 * roww + (size_t)ch * 4;
        wmma::fragment<wmma::matrix_a, 8, 8, 128, wmma::experimental::precision::b1, wmma::row_major> a_frag;
        wmma::fragment<wmma::matrix_b, 8, 8, 128, wmma::experimental::precision::b1, wmma::col_major> b_frag;
        wmma::load_matrix_sync(a_frag, ap, ldm);
        wmma::load_matrix_sync(b_frag, bp, ldm);
        wmma::bmma_sync(c_frag, a_frag, b_frag, c_frag,
                        wmma::experimental::bmmaBitOpAND,
                        wmma::experimental::bmmaAccumulateOpPOPC);
    }
    wmma::store_matrix_sync(&smem[warp][0], c_frag, 8, wmma::mem_row_major);
    __syncthreads();
    if (threadIdx.x < 64) {
        const int ii = (int)(threadIdx.x >> 3);
        const int jj = (int)(threadIdx.x & 7);
        const int i = ti * 8 + ii;
        const int j = tj * 8 + jj;
        if (i >= j || i >= M || j >= M) return;
        int c11 = 0;
        #pragma unroll
        for (int w = 0; w < 4; w++) c11 += smem[w][ii * 8 + jj];
        if (screen_keep(c11, sup[i], sup[j], N, rawmode)) {
            unsigned long long idx = atomicAdd(counter, 1ULL);
            if (idx < cap) { out[idx].i = i; out[idx].j = j; out[idx].c11 = c11; }
            else *overflow = 1;
        }
    }
}

__global__ void simt_full_kernel(const uint32_t* __restrict__ bm,
                                 const uint64_t* __restrict__ sup,
                                 int M, int N, int Npad,
                                 Triple* __restrict__ out,
                                 unsigned long long* __restrict__ counter,
                                 unsigned long long cap,
                                 int* __restrict__ overflow,
                                 int rawmode) {
    const int ti = blockIdx.x;
    const int tj = blockIdx.y;
    if (tj < ti) return;
    if (threadIdx.x >= 64) return;
    const int roww = Npad / 32;
    const int ii = (int)(threadIdx.x >> 3);
    const int jj = (int)(threadIdx.x & 7);
    const int i = ti * 8 + ii;
    const int j = tj * 8 + jj;
    if (i >= j || i >= M || j >= M) return;
    const uint32_t* a = bm + (size_t)i * roww;
    const uint32_t* b = bm + (size_t)j * roww;
    const int wn = (N + 31) / 32;
    int c11 = 0;
    for (int w = 0; w < wn; w++) c11 += __popc(a[w] & b[w]);
    if (screen_keep(c11, sup[i], sup[j], N, rawmode)) {
        unsigned long long idx = atomicAdd(counter, 1ULL);
        if (idx < cap) { out[idx].i = i; out[idx].j = j; out[idx].c11 = c11; }
        else *overflow = 1;
    }
}

// ---------------- router ----------------
struct CostTable {
    double setup_simt = 0, rate_simt = 0;   // s per launch, s per chunk-tile
    double setup_bmma = 0, rate_bmma = 0;
    int pref = -1;                          // cost-model preferred engine at current KCH
};

struct Flags {
    int batch_sz = 131072;   // tiles per launch; ablation --no-batch sets 1
    int overlap = 1, hysteresis = 1, refresh = 1, exact = 1;
    int static_route = 0;
    int load_pct = 0;
    int rawmode = 0;
};

struct Phases {
    double calib = 0, route = 0, simt = 0, bmma = 0, tail = 0, sort = 0, verify = 0, io = 0;
    long tiles_simt = 0, tiles_bmma = 0, launches = 0, refreshes = 0;
};

// calibrate: full-grid timing per engine -> per-tile rate; small-grid (3 tiles) -> setup.
// Uses its OWN counter/overflow so the live run's emission counter is never disturbed.
static void calibrate(const uint32_t* dbm, const uint64_t* dsup, int M, int N, int Npad,
                      Triple* doutcal, unsigned long long* dccal, int* dovfcal,
                      unsigned long long cap, int rawmode, cudaStream_t s,
                      double ntiles, CostTable& ct) {
    const int T = (M + 7) / 8;
    const int KCH = Npad / 128;
    double tf[2] = {0, 0};
    for (int eng = 0; eng < 2; eng++) {
        double best = 1e30;
        for (int r = 0; r < 3; r++) {
            CK(cudaMemsetAsync(dccal, 0, sizeof(unsigned long long), s));
            CK(cudaMemsetAsync(dovfcal, 0, sizeof(int), s));
            auto t0 = chrono::steady_clock::now();
            if (eng == 0) simt_full_kernel<<<dim3(T, T), 128, 0, s>>>(dbm, dsup, M, N, Npad, doutcal, dccal, cap, dovfcal, rawmode);
            else bmma_full_kernel<<<dim3(T, T), 128, 0, s>>>(dbm, dsup, M, N, Npad, doutcal, dccal, cap, dovfcal, rawmode);
            CK(cudaGetLastError());
            CK(cudaStreamSynchronize(s));
            best = min(best, chrono::duration<double>(chrono::steady_clock::now() - t0).count());
        }
        tf[eng] = best / ntiles;
    }
    ct.rate_simt = tf[0] / KCH;
    ct.rate_bmma = tf[1] / KCH;
    // setup: small grid = 3 upper-tri tiles
    double sg[2] = {0, 0};
    for (int eng = 0; eng < 2; eng++) {
        double best = 1e30;
        for (int r = 0; r < 5; r++) {
            CK(cudaMemsetAsync(dccal, 0, sizeof(unsigned long long), s));
            CK(cudaMemsetAsync(dovfcal, 0, sizeof(int), s));
            auto t0 = chrono::steady_clock::now();
            if (eng == 0) simt_full_kernel<<<dim3(2, 2), 128, 0, s>>>(dbm, dsup, M, N, Npad, doutcal, dccal, cap, dovfcal, rawmode);
            else bmma_full_kernel<<<dim3(2, 2), 128, 0, s>>>(dbm, dsup, M, N, Npad, doutcal, dccal, cap, dovfcal, rawmode);
            CK(cudaGetLastError());
            CK(cudaStreamSynchronize(s));
            best = min(best, chrono::duration<double>(chrono::steady_clock::now() - t0).count());
        }
        sg[eng] = best / 3.0;
    }
    ct.setup_simt = max(0.0, sg[0] - ct.rate_simt * KCH);
    ct.setup_bmma = max(0.0, sg[1] - ct.rate_bmma * KCH);
    const double cs = ct.setup_simt + ct.rate_simt * KCH;
    const double cb = ct.setup_bmma + ct.rate_bmma * KCH;
    ct.pref = (cs < cb) ? 0 : 1;
}

int main(int argc, char** argv) {
    vector<string> args(argv + 1, argv + argc);
    Flags fl;
    for (size_t k = 0; k < args.size(); k++) {
        if (args[k] == "--raw") fl.rawmode = 1;
        else if (args[k] == "--no-batch") fl.batch_sz = 8192;  // batch=1 (23M launches) launch-bound infeasible; 512 shows the overhead honestly
        else if (args[k] == "--no-overlap") fl.overlap = 0;
        else if (args[k] == "--no-hysteresis") fl.hysteresis = 0;
        else if (args[k] == "--static-table") fl.refresh = 0;
        else if (args[k] == "--static-route") fl.static_route = 1;
        else if (args[k] == "--no-exact") fl.exact = 0;
        else if (args[k] == "--load" && k + 1 < args.size()) fl.load_pct = atoi(args[++k].c_str());
    }
    if (args.empty() || (args[0] != "--count" && args[0] != "--bench" && args[0] != "--sweep")) {
        fprintf(stderr, "usage: gbm_router --count in out [flags] | --bench N M seed reps [flags] | --sweep in [flags]\n");
        return 2;
    }

    Matrix m;
    int reps = 1;
    uint64_t bench_seed = 0;
    if (args[0] == "--bench") {
        m.N = atoi(args[1].c_str()); m.M = atoi(args[2].c_str());
        bench_seed = strtoull(args[3].c_str(), 0, 10);
        reps = atoi(args[4].c_str());
        m.words = (m.N + 63) / 64;
        m.bits.assign((size_t)m.M * m.words, 0ULL);
        m.support.assign(m.M, 0);
        mt19937_64 rng(bench_seed);
        for (int j = 0; j < m.M; j++)
            for (int i = 0; i < m.N; i++) {
                bool one = (double)rng() / (double)UINT64_MAX < 0.1;
                if (one) { m.bits[(size_t)j * m.words + (i >> 6)] |= (1ULL << (i & 63)); m.support[j]++; }
            }
        printf("bench N=%d M=%d seed=%llu reps=%d raw=%d flags(batch=%d overlap=%d hyst=%d refresh=%d static=%d exact=%d load=%d%%)\n",
               m.N, m.M, bench_seed, reps, fl.rawmode, fl.batch_sz, fl.overlap, fl.hysteresis,
               fl.refresh, fl.static_route, fl.exact, fl.load_pct);
    } else {
        if (!load_txt(m, args[1].c_str())) return 1;
        if (args[0] == "--count") printf("loaded N=%d M=%d mode=dual-engine-router\n", m.N, m.M);
    }

    const int Npad = ((m.N + 127) / 128) * 128;
    const int ldm = Npad / 32;
    vector<uint32_t> hostbm((size_t)m.M * ldm, 0u);
    for (int f = 0; f < m.M; f++) {
        uint32_t* dst = &hostbm[(size_t)f * ldm];
        const uint64_t* src = &m.bits[(size_t)f * m.words];
        const int w32 = (m.N + 31) / 32;
        for (int w = 0; w < w32; w++) {
            uint32_t v = 0;
            const int w64 = w >> 1;
            if ((w & 1) == 0) v = (uint32_t)(src[w64] & 0xffffffffULL);
            else v = (uint32_t)(src[w64] >> 32);
            dst[w] = v;
        }
    }
    m.bits.clear(); m.bits.shrink_to_fit();

    uint32_t* dbm = nullptr;
    CK(cudaMalloc(&dbm, (size_t)m.M * ldm * sizeof(uint32_t)));
    CK(cudaMemcpy(dbm, hostbm.data(), (size_t)m.M * ldm * sizeof(uint32_t), cudaMemcpyHostToDevice));
    hostbm.clear(); hostbm.shrink_to_fit();
    uint64_t* dsup = nullptr;
    CK(cudaMalloc(&dsup, (size_t)m.M * sizeof(uint64_t)));
    CK(cudaMemcpy(dsup, m.support.data(), (size_t)m.M * sizeof(uint64_t), cudaMemcpyHostToDevice));
    const unsigned long long cap = 1ULL << 28;
    Triple* dout = nullptr;
    CK(cudaMalloc(&dout, cap * sizeof(Triple)));
    unsigned long long* dcnt = nullptr;
    CK(cudaMalloc(&dcnt, sizeof(unsigned long long)));
    int* dovf = nullptr;
    CK(cudaMalloc(&dovf, sizeof(int)));
    unsigned long long* dccal = nullptr;
    CK(cudaMalloc(&dccal, sizeof(unsigned long long)));
    int* dovfcal = nullptr;
    CK(cudaMalloc(&dovfcal, sizeof(int)));
    int2* dlist[2] = {nullptr, nullptr};
    const size_t dlist_cap = (size_t)((m.M + 7) / 8) * ((size_t)(m.M + 7) / 8 + 1) / 2 + 64;
    CK(cudaMalloc(&dlist[0], sizeof(int2) * dlist_cap));
    CK(cudaMalloc(&dlist[1], sizeof(int2) * dlist_cap));
    Triple* doutcal = nullptr;
    const unsigned long long cap_cal = 1ULL << 22;   // 16 MB calibration scratch, never aliases live output
    CK(cudaMalloc(&doutcal, cap_cal * sizeof(Triple)));

    cudaStream_t s_bmma, s_simt, s_noise;
    CK(cudaStreamCreate(&s_bmma));
    CK(cudaStreamCreate(&s_simt));
    CK(cudaStreamCreate(&s_noise));
    float* dnoise = nullptr;
    if (fl.load_pct > 0) CK(cudaMalloc(&dnoise, 1 << 22));

    const int T = (m.M + 7) / 8;
    const int KCH = Npad / 128;
    const double ntiles = (double)T * (T + 1) / 2.0;

    // ---------------- --sweep: dark-region measurement (E-C03) ----------------
    if (args[0] == "--sweep") {
        printf("sweep N=%d M=%d Npad=%d KCH=%d tiles=%.0f\n", m.N, m.M, Npad, KCH, ntiles);
        printf("kchunks\tnpad\tbmma_s\tsimt_s\tratio(simt/bmma)\trouter_choice\n");
        for (int kch = 1; kch <= KCH; kch *= 2) {
            const int np2 = min(kch * 128, Npad);
            const uint32_t* bm2 = dbm;
            uint32_t* dtmp = nullptr;
            const int ldm2 = np2 / 32;
            const int N2 = min(np2, m.N);
            if (np2 < Npad) {
                // cropped copy: first ceil(np2/32) words per row
                CK(cudaMalloc(&dtmp, (size_t)m.M * ldm2 * sizeof(uint32_t)));
                // rows are contiguous with stride ldm words; copy prefix per row on host is
                // simpler via a strided kernel
                CK(cudaMemcpy2DAsync(dtmp, ldm2 * sizeof(uint32_t), dbm, ldm * sizeof(uint32_t),
                                     ldm2 * sizeof(uint32_t), m.M, cudaMemcpyDeviceToDevice, s_bmma));
                CK(cudaStreamSynchronize(s_bmma));
                bm2 = dtmp;
            }
            double tb = 0, ts = 0;
            for (int eng = 0; eng < 2; eng++) {
                double best = 1e30;
                for (int r = 0; r < 3; r++) {
                    CK(cudaMemsetAsync(dcnt, 0, sizeof(unsigned long long), s_bmma));
                    CK(cudaMemsetAsync(dovf, 0, sizeof(int), s_bmma));
                    auto t0 = chrono::steady_clock::now();
                    if (eng == 0) simt_full_kernel<<<dim3(T, T), 128, 0, s_bmma>>>(bm2, dsup, m.M, N2, np2, dout, dcnt, cap, dovf, fl.rawmode);
                    else bmma_full_kernel<<<dim3(T, T), 128, 0, s_bmma>>>(bm2, dsup, m.M, N2, np2, dout, dcnt, cap, dovf, fl.rawmode);
                    CK(cudaGetLastError());
                    CK(cudaStreamSynchronize(s_bmma));
                    best = min(best, chrono::duration<double>(chrono::steady_clock::now() - t0).count());
                }
                if (eng == 0) ts = best; else tb = best;
            }
            printf("%d\t%d\t%.4f\t%.4f\t%.3f\t%s\n", kch, np2, tb, ts, ts / max(tb, 1e-9),
                   (ts < tb) ? "SIMT(matrix-dark)" : "BMMA");
            if (dtmp) cudaFree(dtmp);
            if (kch * 2 > KCH && kch != KCH) { /* loop ends by *=2 */ }
        }
        return 0;
    }

    // ---------------- routed counting ----------------
    // tile list (upper triangle incl diagonal), grouped so same-engine runs are contiguous
    vector<pair<int,int>> tiles;
    for (int tj = 0; tj < T; tj++)
        for (int ti = 0; ti <= tj; ti++) tiles.push_back({ti, tj});
    const size_t n_all = tiles.size();

    for (int rep = 0; rep < reps; rep++) {
        Phases ph;
        auto t_all0 = chrono::steady_clock::now();
        CK(cudaMemset(dcnt, 0, sizeof(unsigned long long)));
        CK(cudaMemset(dovf, 0, sizeof(int)));

        auto tc0 = chrono::steady_clock::now();
        CostTable ct;
        calibrate(dbm, dsup, m.M, m.N, Npad, doutcal, dccal, dovfcal, cap_cal, fl.rawmode, s_bmma, ntiles, ct);
        ph.calib = chrono::duration<double>(chrono::steady_clock::now() - tc0).count();
        CK(cudaMemset(dcnt, 0, sizeof(unsigned long long)));
        CK(cudaMemset(dovf, 0, sizeof(int)));

        // routing decision per tile (uniform KCH -> uniform choice; hysteresis shown)
        auto t_route0 = chrono::steady_clock::now();
        const double HYST = fl.hysteresis ? 1.15 : 1.0;
        const double cs = ct.setup_simt + ct.rate_simt * KCH;
        const double cb = ct.setup_bmma + ct.rate_bmma * KCH;
        int engine; // 0=simt 1=bmma
        if (fl.static_route) engine = ct.pref;
        else {
            // challenger must beat incumbent by hysteresis margin
            if (ct.pref == 0) engine = (cb < cs * HYST) ? 1 : 0;
            else engine = (cs < cb * HYST) ? 0 : 1;
        }
        ph.route = chrono::duration<double>(chrono::steady_clock::now() - t_route0).count();
        (void)HYST;

        // fluctuating background load (E-F04)
        volatile bool stop_noise = false;
        thread* noise_th = nullptr;
        if (fl.load_pct > 0) {
            noise_th = new thread([&]() {
                const int period_ms = 2000;
                const int on_ms = period_ms * fl.load_pct / 100;
                while (!stop_noise) {
                    noise_kernel<<<512, 256, 0, s_noise>>>(dnoise, 1 << 20, 4000);
                    this_thread::sleep_for(chrono::milliseconds(on_ms));
                    CK(cudaStreamSynchronize(s_noise));
                    this_thread::sleep_for(chrono::milliseconds(period_ms - on_ms));
                }
            });
        }

        cudaEvent_t ev_s0, ev_s1, ev_b0, ev_b1;
        cudaEventCreate(&ev_s0); cudaEventCreate(&ev_s1);
        cudaEventCreate(&ev_b0); cudaEventCreate(&ev_b1);
        cudaEventRecord(ev_s0, s_simt);
        cudaEventRecord(ev_b0, s_bmma);

        auto t_k0 = chrono::steady_clock::now();
        // build engine-contiguous batches: all tiles same engine here, so batch
        // the full list in chunks of batch_sz and launch on the engine's stream
        vector<int2> hlist(fl.batch_sz);
        const int nchunks_launch = (int)((n_all + fl.batch_sz - 1) / fl.batch_sz);
        for (int c = 0; c < nchunks_launch; c++) {
            const size_t k0 = (size_t)c * fl.batch_sz;
            const int nb = (int)min((size_t)fl.batch_sz, n_all - k0);
            for (int q = 0; q < nb; q++) {
                hlist[q].x = tiles[k0 + q].first;
                hlist[q].y = tiles[k0 + q].second;
            }
            int2* dl = dlist[engine];
            CK(cudaMemcpyAsync(dl, hlist.data(), nb * sizeof(int2), cudaMemcpyHostToDevice,
                               engine == 1 ? s_bmma : s_simt));
            if (engine == 1) {
                bmma_list_kernel<<<nb, 128, 0, s_bmma>>>(dbm, dsup, m.M, m.N, Npad, dl, nb,
                                                         dout, dcnt, cap, dovf, fl.rawmode);
                ph.tiles_bmma += nb;
            } else {
                simt_list_kernel<<<nb, 128, 0, s_simt>>>(dbm, dsup, m.M, m.N, Npad, dl, nb,
                                                         dout, dcnt, cap, dovf, fl.rawmode);
                ph.tiles_simt += nb;
            }
            ph.launches++;
            if (!fl.overlap) CK(cudaDeviceSynchronize());
            // periodic cheap cost-table refresh (setup-only re-probe, every 256 launches)
            if (fl.refresh && ph.launches % 256 == 0) {
                CostTable ct2;
                calibrate(dbm, dsup, m.M, m.N, Npad, doutcal, dccal, dovfcal, cap_cal, fl.rawmode, s_bmma, ntiles, ct2);
                ph.refreshes++;
                if (!fl.static_route) {
                    const double cs2 = ct2.setup_simt + ct2.rate_simt * KCH;
                    const double cb2 = ct2.setup_bmma + ct2.rate_bmma * KCH;
                    if (ct2.pref == 0) engine = (cb2 < cs2 * HYST) ? 1 : 0;
                    else engine = (cs2 < cb2 * HYST) ? 0 : 1;
                }
            }
        }
        cudaEventRecord(ev_s1, s_simt);
        cudaEventRecord(ev_b1, s_bmma);
        auto t_k1 = chrono::steady_clock::now();
        CK(cudaDeviceSynchronize());
        auto t_sync1 = chrono::steady_clock::now();
        cudaEventSynchronize(ev_s1); cudaEventSynchronize(ev_b1);
        float ms_simt = 0, ms_bmma = 0;
        cudaEventElapsedTime(&ms_simt, ev_s0, ev_s1);
        cudaEventElapsedTime(&ms_bmma, ev_b0, ev_b1);
        ph.simt = ms_simt / 1e3; ph.bmma = ms_bmma / 1e3;
        ph.tail = chrono::duration<double>(t_sync1 - t_k1).count();
        if (noise_th) { stop_noise = true; noise_th->join(); delete noise_th; noise_th = nullptr; }

        unsigned long long emitted = 0; int ovf = 0;
        CK(cudaMemcpy(&emitted, dcnt, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
        CK(cudaMemcpy(&ovf, dovf, sizeof(int), cudaMemcpyDeviceToHost));
        if (ovf) { fprintf(stderr, "capacity exceeded\n"); return 1; }

        vector<Triple> buf(emitted);
        if (emitted) CK(cudaMemcpy(buf.data(), dout, emitted * sizeof(Triple), cudaMemcpyDeviceToHost));
        auto t_s0 = chrono::steady_clock::now();
        vector<uint64_t> off((size_t)m.M + 1, 0);
        for (size_t k = 0; k < buf.size(); k++) off[buf[k].i + 1]++;
        for (int i = 0; i < m.M; i++) off[i + 1] += off[i];
        vector<Triple> srt(buf.size());
        vector<uint64_t> cur(off.begin(), off.end() - 1);
        for (size_t k = 0; k < buf.size(); k++) srt[cur[buf[k].i]++] = buf[k];
        for (int i = 0; i < m.M; i++)
            sort(srt.begin() + off[i], srt.begin() + off[i + 1],
                 [](const Triple& x, const Triple& y) { return x.j < y.j; });
        ph.sort = chrono::duration<double>(chrono::steady_clock::now() - t_s0).count();

        auto t_v0 = chrono::steady_clock::now();
        if (fl.exact) {
            for (size_t k = 1; k < srt.size(); k++) {
                if (srt[k-1].i > srt[k].i || (srt[k-1].i == srt[k].i && srt[k-1].j >= srt[k].j)) {
                    fprintf(stderr, "EXACTNESS-PROTOCOL-FAIL at %zu\n", k); return 1;
                }
            }
        }
        ph.verify = chrono::duration<double>(chrono::steady_clock::now() - t_v0).count();

        if (args[0] == "--count" && rep == reps - 1) {
            auto t_io0 = chrono::steady_clock::now();
            FILE* f = fopen(args[2].c_str(), "w");
            if (!f) { fprintf(stderr, "cannot write %s\n", args[2].c_str()); return 1; }
            fprintf(f, "famA\tfamB\tc11\tsupA\tsupB\tc10\tc01\tc00\n");
            for (size_t k = 0; k < srt.size(); k++) {
                const Triple& t = srt[k];
                const uint64_t supA = m.support[t.i], supB = m.support[t.j];
                const uint64_t c11 = (uint64_t)t.c11;
                fprintf(f, "%s\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
                        m.famNames[t.i].c_str(), m.famNames[t.j].c_str(),
                        c11, supA, supB, supA - c11, supB - c11,
                        (uint64_t)m.N - supA - supB + c11);
            }
            fclose(f);
            ph.io = chrono::duration<double>(chrono::steady_clock::now() - t_io0).count();
        }
        const double total = chrono::duration<double>(chrono::steady_clock::now() - t_all0).count();
        printf("rep%d ROUTE-DONE tiles=%zu engine=%s(simt=0/bmma=1) emitted=%llu "
               "calib=%.3f route=%.4f simt=%.3fs bmma=%.3fs tail=%.3f sort=%.3f verify=%.4f io=%.3f total=%.3fs "
               "cost(simt s=%.2e r=%.2e | bmma s=%.2e r=%.2e) launches=%ld refreshes=%ld\n",
               rep, n_all, engine == 1 ? "bmma" : "simt", emitted,
               ph.calib, ph.route, ph.simt, ph.bmma, ph.tail, ph.sort, ph.verify, ph.io, total,
               ct.setup_simt, ct.rate_simt, ct.setup_bmma, ct.rate_bmma,
               ph.launches, ph.refreshes);
        cudaEventDestroy(ev_s0); cudaEventDestroy(ev_s1);
        cudaEventDestroy(ev_b0); cudaEventDestroy(ev_b1);
    }

    cudaFree(dbm); cudaFree(dsup); cudaFree(dout); cudaFree(dcnt); cudaFree(dovf); cudaFree(dccal); cudaFree(dovfcal); cudaFree(dlist[0]); cudaFree(dlist[1]); cudaFree(doutcal);
    if (dnoise) cudaFree(dnoise);
    cudaStreamDestroy(s_bmma); cudaStreamDestroy(s_simt); cudaStreamDestroy(s_noise);
    return 0;
}
