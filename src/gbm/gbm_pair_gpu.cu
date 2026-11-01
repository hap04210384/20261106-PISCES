// gbm_pair_gpu.cu - Gene Bitmap Mining: all-pairs exact co-occurrence counting on Tensor Cores (BMMA wave)
// Genome x gene-family presence/absence bitmap engine (WMMA b1, CUDA, Release).
// Tensor-core path: each mma.m8n8k128.b1 computes an 8-family x 8-family chunk
// co-occurrence count (popcount of pairwise AND over 128 genomes); chunks are
// accumulated in an int32 fragment, giving exact c11 for all N genomes.
// The host-side hypergeometric screen replicates gbm_count.cpp bit-for-bit
// (same expression order, IEEE double, no fast-math) so the emitted TSV is
// byte-identical to the AVX2 CPU reference.
// Part of the Nature Methods submission reproducibility package.
// Author: YalongZhang
//
// Modes:
//   gbm_pair_gpu --count in.txt out.tsv [capacity_entries]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
using namespace std;

#include <cuda_runtime.h>
#include <mma.h>
using namespace nvcuda;

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(1);} } while (0)

// ---------------- GBM1 container (same on-disk format as gbm_count) ----------------
struct Matrix {
    int N = 0, M = 0;
    int words = 0;                 // uint64 words per family row (unpadded)
    vector<uint64_t> bits;         // M x words (unpadded rows)
    vector<string> famNames;
    vector<uint64_t> support;        // per-family popcount over the FIRST N bits
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

__global__ void bmma_count_kernel(const uint32_t* __restrict__ bm,   // families x Npad bits, 256B-aligned rows
                                  const uint64_t* __restrict__ sup,  // per-family support (over real N bits)
                                  int M, int N, int Npad,
                                  Triple* __restrict__ out,
                                  unsigned long long* __restrict__ counter,
                                  unsigned long long cap,
                                  int* __restrict__ overflow,
                                  int rawmode) {
    const int T = (M + 7) / 8;
    const int ti = blockIdx.x;   // family tile (rows of A)
    const int tj = blockIdx.y;   // family tile (cols of B)
    if (tj < ti) return;         // upper triangle only
    const int KCH = Npad / 128;     // 128-genome chunks
    const int ldm = Npad;           // b1 load ldm is in BITS; row addr = p + i*ldm/8 (16B-aligned rows)
    const int roww = Npad / 32;     // uint32 words per family row

    __shared__ int smem[4][64];
    const int warp = threadIdx.x >> 5;
    const int lane = threadIdx.x & 31;

    wmma::fragment<wmma::accumulator, 8, 8, 128, int> c_frag;
    wmma::fill_fragment(c_frag, 0);

    for (int ch = warp; ch < KCH; ch += 4) {
        // A: 8 family rows x 128 genome bits (row-major)
        const uint32_t* ap = bm + (size_t)ti * 8 * roww + (size_t)ch * 4;
        // B: 8 family rows x 128 genome bits, read as col-major 128x8
        //     (column j of B = genome-axis copy of family j's chunk)
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

    // 64 (ii,jj) slots; threads 0..63 each finalize one pair.
    // Screen replicates gbm_count.cpp line-for-line (IEEE double, same order).
    if (threadIdx.x < 64) {
        const int ii = (int)(threadIdx.x >> 3);
        const int jj = (int)(threadIdx.x & 7);
        const int i = ti * 8 + ii;
        const int j = tj * 8 + jj;
        if (i >= j || i >= M || j >= M) return;
        int c11 = 0;
        #pragma unroll
        for (int w = 0; w < 4; w++) c11 += smem[w][ii * 8 + jj];
        const double n  = (double)N;
        const double r1 = (double)sup[i], r2 = (double)sup[j];
        const double E = r1 * r2 / n;
        const double var = r1 * r2 * (n - r1) * (n - r2) / (n * n * (n - 1.0));
        bool keep;
        if (rawmode) {
            keep = (c11 >= 1);
        } else {
            keep = false;
            if (var > 0) {
                const double sd = sqrt(var);
                keep = fabs((double)c11 - E) >= 3.29 * sd;
            }
        }
        if (keep) {
            unsigned long long idx = atomicAdd(counter, 1ULL);
            if (idx < cap) {
                out[idx].i = i; out[idx].j = j; out[idx].c11 = c11;
            } else {
                *overflow = 1;
            }
        }
    }
}

// ---------------- host side ----------------
int main(int argc, char** argv) {
    if (argc >= 2 && strcmp(argv[1], "--count") == 0) {
        const char* in = argv[2];
        const char* outp = argv[3];
        unsigned long long cap = (argc >= 5) ? strtoull(argv[4], 0, 10) : (1ULL << 26);
        const int rawmode = (argc >= 6 && strcmp(argv[5], "raw") == 0) ? 1 : 0;
        Matrix m;
        if (!load_txt(m, in)) return 1;
        const int Npad = ((m.N + 127) / 128) * 128;
        const int ldm = Npad / 32;  // uint32 per family row
        printf("loaded N=%d M=%d Npad=%d mode=bmma-wave\n", m.N, m.M, Npad);

        // build padded device bitmap (zero padding contributes nothing to counts)
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
        uint32_t* dbm = nullptr;
        CK(cudaMalloc(&dbm, (size_t)m.M * ldm * sizeof(uint32_t)));
        CK(cudaMemcpy(dbm, hostbm.data(), (size_t)m.M * ldm * sizeof(uint32_t), cudaMemcpyHostToDevice));
        uint64_t* dsup = nullptr;
        CK(cudaMalloc(&dsup, (size_t)m.M * sizeof(uint64_t)));
        CK(cudaMemcpy(dsup, m.support.data(), (size_t)m.M * sizeof(uint64_t), cudaMemcpyHostToDevice));
        Triple* dout = nullptr;
        CK(cudaMalloc(&dout, cap * sizeof(Triple)));
        unsigned long long* dcnt = nullptr;
        CK(cudaMalloc(&dcnt, sizeof(unsigned long long)));
        CK(cudaMemset(dcnt, 0, sizeof(unsigned long long)));
        int* dovf = nullptr;
        CK(cudaMalloc(&dovf, sizeof(int)));
        CK(cudaMemset(dovf, 0, sizeof(int)));

        const int T = (m.M + 7) / 8;
        dim3 grid(T, T);
        auto t0 = chrono::steady_clock::now();
        bmma_count_kernel<<<grid, 128>>>(dbm, dsup, m.M, m.N, Npad, dout, dcnt, cap, dovf, rawmode);
        CK(cudaGetLastError());
        CK(cudaDeviceSynchronize());
        auto t1 = chrono::steady_clock::now();
        const double kernel_s = chrono::duration<double>(t1 - t0).count();

        unsigned long long emitted = 0;
        int ovf = 0;
        CK(cudaMemcpy(&emitted, dcnt, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
        CK(cudaMemcpy(&ovf, dovf, sizeof(int), cudaMemcpyDeviceToHost));
        printf("kernel=%.3fs  emitted=%llu  overflow=%d\n", kernel_s, emitted, ovf);
        if (ovf) { fprintf(stderr, "capacity exceeded, raise capacity_entries\n"); return 1; }

        vector<Triple> buf(emitted);
        if (emitted) CK(cudaMemcpy(buf.data(), dout, emitted * sizeof(Triple), cudaMemcpyDeviceToHost));

        // deterministic (i,j) order via bucketed counting sort
        vector<uint64_t> off((size_t)m.M + 1, 0);
        for (size_t k = 0; k < buf.size(); k++) off[buf[k].i + 1]++;
        for (int i = 0; i < m.M; i++) off[i + 1] += off[i];
        vector<Triple> srt(buf.size());
        vector<uint64_t> cur(off.begin(), off.end() - 1);
        for (size_t k = 0; k < buf.size(); k++) srt[cur[buf[k].i]++] = buf[k];
        // secondary sort by j inside each i-bucket for deterministic output
        for (int i = 0; i < m.M; i++)
            sort(srt.begin() + off[i], srt.begin() + off[i + 1],
                 [](const Triple& x, const Triple& y) { return x.j < y.j; });
        auto t2 = chrono::steady_clock::now();

        FILE* f = fopen(outp, "w");
        if (!f) { fprintf(stderr, "cannot write %s\n", outp); return 1; }
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
        auto t3 = chrono::steady_clock::now();
        const double total_s = chrono::duration<double>(t3 - t0).count();
        printf("GPU-COUNT-DONE pairs=%llu emitted=%llu kernel=%.3fs sort=%.3fs io=%.3fs total=%.2fs\n",
               (unsigned long long)((uint64_t)m.M * (m.M - 1) / 2), emitted,
               kernel_s,
               chrono::duration<double>(t2 - t1).count(),
               chrono::duration<double>(t3 - t2).count(), total_s);

        cudaFree(dbm); cudaFree(dsup); cudaFree(dout); cudaFree(dcnt); cudaFree(dovf);
        return 0;
    }
    fprintf(stderr, "usage: gbm_pair_gpu --count in.txt out.tsv [capacity_entries]\n");
    return 2;
}
