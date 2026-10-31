// gbm_count.cpp - Gene Bitmap Mining: all-pairs exact co-occurrence counting
// Genome x gene-family presence/absence bitmap engine (AVX2, Release).
// Part of the Nature Methods submission reproducibility package.
// Author: YalongZhang
//
// Modes:
//   gbm_count --gen N M seed out.txt          generate random 0/1 matrix (families x genomes)
//   gbm_count --count in.txt out.tsv [minc]   all-pairs exact count (c11 >= minc), marginals from support index
//   gbm_count --selftest N M seed             E-A01: AVX2 bitmap vs naive 3-loop, all pairs must agree
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <chrono>
using namespace std;

// ---------------- AVX2 popcount (Harley-Seal over 256-bit) ----------------
#include <immintrin.h>

static inline uint64_t popcnt256(__m256i v) {
    // split into two 128-bit halves and use POPCNTQ
    __m128i lo = _mm256_castsi256_si128(v);
    __m128i hi = _mm256_extracti128_si256(v, 1);
    return (uint64_t)_mm_popcnt_u64((uint64_t)_mm_cvtsi128_si64(lo)) +
           (uint64_t)_mm_popcnt_u64((uint64_t)_mm_cvtsi128_si64(hi)) +
           (uint64_t)_mm_popcnt_u64((uint64_t)_mm_extract_epi64(lo, 1)) +
           (uint64_t)_mm_popcnt_u64((uint64_t)_mm_extract_epi64(hi, 1));
}

static uint64_t and_popcnt_rows(const uint64_t* a, const uint64_t* b, int words) {
    uint64_t s = 0;
    int i = 0;
    for (; i + 4 <= words; i += 4) {
        __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
        __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
        s += popcnt256(_mm256_and_si256(va, vb));
    }
    for (; i < words; i++) s += (uint64_t)_mm_popcnt_u64(a[i] & b[i]);
    return s;
}

// ---------------- matrix container ----------------
struct Matrix {
    int N = 0;                 // genomes
    int M = 0;                 // families
    int words = 0;             // uint64 words per family row
    vector<uint64_t> bits;     // M x words
    vector<string> famNames;
    vector<uint64_t> support;  // per-family popcount

    uint64_t* row(int f) { return bits.data() + (size_t)f * words; }
    const uint64_t* row(int f) const { return bits.data() + (size_t)f * words; }
};

static bool save_txt(const Matrix& m, const char* path) {
    FILE* f = fopen(path, "w");
    if (!f) return false;
    fprintf(f, "GBM1 %d %d\n", m.N, m.M);
    for (int j = 0; j < m.M; j++) {
        fprintf(f, "%s ", m.famNames[j].c_str());
        const uint64_t* r = m.row(j);
        for (int i = 0; i < m.N; i++) {
            int bit = (int)((r[i >> 6] >> (i & 63)) & 1ULL);
            fputc(bit ? '1' : '0', f);
        }
        fputc('\n', f);
    }
    fclose(f);
    return true;
}

static bool load_txt(Matrix& m, const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return false; }
    char magic[8];
    if (fscanf(f, "%7s", magic) != 1 || strcmp(magic, "GBM1") != 0) {
        fprintf(stderr, "bad magic\n"); fclose(f); return false;
    }
    if (fscanf(f, "%d %d", &m.N, &m.M) != 2) { fprintf(stderr, "bad header\n"); fclose(f); return false; }
    m.words = (m.N + 63) / 64;
    m.bits.assign((size_t)m.M * m.words, 0ULL);
    m.famNames.resize(m.M);
    m.support.assign(m.M, 0);
    char name[256];
    for (int j = 0; j < m.M; j++) {
        if (fscanf(f, "%255s", name) != 1) { fprintf(stderr, "truncated at family %d\n", j); fclose(f); return false; }
        m.famNames[j] = name;
        // read N bits: consume whitespace-separated? our format has bits immediately after name
        int c, i = 0;
        // skip one separator char (space/tab)
        while ((c = fgetc(f)) == ' ' || c == '\t') {}
        if (c != '0' && c != '1' && c != '\n') { fprintf(stderr, "unexpected char %c at family %d\n", c, j); fclose(f); return false; }
        while (c == '0' || c == '1') {
            if (c == '1' && i < m.N) { m.row(j)[i >> 6] |= (1ULL << (i & 63)); m.support[j]++; }
            i++;
            c = fgetc(f);
        }
        if (i != m.N) { fprintf(stderr, "family %d has %d bits, expected %d\n", m.famNames[j].c_str() ? j : j, i, m.N); fclose(f); return false; }
    }
    fclose(f);
    return true;
}

static Matrix gen_random(int N, int M, uint64_t seed, double density = 0.1) {
    Matrix m;
    m.N = N; m.M = M;
    m.words = (N + 63) / 64;
    m.bits.assign((size_t)M * m.words, 0ULL);
    m.famNames.resize(M);
    m.support.assign(M, 0);
    mt19937_64 rng(seed);
    for (int j = 0; j < M; j++) {
        char buf[64];
        snprintf(buf, sizeof buf, "fam_%06d", j);
        m.famNames[j] = buf;
        for (int i = 0; i < N; i++) {
            bool one = (double)rng() / (double)UINT64_MAX < density;
            if (one) { m.row(j)[i >> 6] |= (1ULL << (i & 63)); m.support[j]++; }
        }
    }
    return m;
}

// naive count for self-test
static uint64_t naive_and(const Matrix& m, int a, int b) {
    uint64_t s = 0;
    const uint64_t* ra = m.row(a);
    const uint64_t* rb = m.row(b);
    for (int w = 0; w < m.words; w++) s += (uint64_t)_mm_popcnt_u64(ra[w] & rb[w]);
    return s;
}

int main(int argc, char** argv) {
    if (argc >= 2 && strcmp(argv[1], "--gen") == 0) {
        int N = atoi(argv[2]), M = atoi(argv[3]);
        uint64_t seed = strtoull(argv[4], 0, 10);
        Matrix m = gen_random(N, M, seed);
        if (!save_txt(m, argv[5])) { fprintf(stderr, "save failed\n"); return 1; }
        printf("generated N=%d M=%d -> %s\n", N, M, argv[5]);
        return 0;
    }
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) {
        int N = atoi(argv[2]), M = atoi(argv[3]);
        uint64_t seed = strtoull(argv[4], 0, 10);
        Matrix m = gen_random(N, M, seed);
        // verify support index
        for (int j = 0; j < M; j++) {
            uint64_t s = 0;
            for (int w = 0; w < m.words; w++) s += (uint64_t)_mm_popcnt_u64(m.row(j)[w]);
            if (s != m.support[j]) { printf("SUPPORT MISMATCH fam %d: %llu vs %llu\n", j, s, m.support[j]); return 1; }
        }
        // all-pairs: AVX2 loop vs naive
        for (int a = 0; a < M; a++) {
            for (int b = a + 1; b < M; b++) {
                uint64_t fast = and_popcnt_rows(m.row(a), m.row(b), m.words);
                uint64_t slow = naive_and(m, a, b);
                if (fast != slow) { printf("MISMATCH pair %d %d: %llu vs %llu\n", a, b, fast, slow); return 1; }
            }
        }
        printf("SELFTEST-OK N=%d M=%d pairs=%lld\n", N, M, (long long)M * (M - 1) / 2);
        return 0;
    }
    if (argc >= 2 && (strcmp(argv[1], "--count") == 0 || strcmp(argv[1], "--countpruned") == 0)) {
        bool pruned = strcmp(argv[1], "--countpruned") == 0;
        const char* in = argv[2];
        const char* out = argv[3];
        uint64_t minc = (argc >= 5) ? strtoull(argv[4], 0, 10) : 2;
        Matrix m;
        if (!load_txt(m, in)) return 1;
        printf("loaded N=%d M=%d words=%d mode=%s\n", m.N, m.M, m.words, pruned ? "pruned" : "full");
        FILE* f = fopen(out, "w");
        if (!f) { fprintf(stderr, "cannot write %s\n", out); return 1; }
        fprintf(f, "famA\tfamB\tc11\tsupA\tsupB\tc10\tc01\tc00\n");
        auto t0 = chrono::steady_clock::now();
        uint64_t emitted = 0, total = 0, pruned_pairs = 0;
        for (int a = 0; a < m.M; a++) {
            if (pruned && m.support[a] < minc) continue;  // anti-monotonic: sup < minc can never reach minc
            const uint64_t* ra = m.row(a);
            for (int b = a + 1; b < m.M; b++) {
                if (pruned && m.support[b] < minc) { pruned_pairs++; continue; }  // upper bound min(supA,supB) < minc
                uint64_t c11 = and_popcnt_rows(ra, m.row(b), m.words);
                total++;
                if (c11 >= minc) {
                    uint64_t supA = m.support[a], supB = m.support[b];
                    fprintf(f, "%s\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
                            m.famNames[a].c_str(), m.famNames[b].c_str(),
                            c11, supA, supB,
                            supA - c11, supB - c11,
                            (uint64_t)m.N - supA - supB + c11);
                    emitted++;
                }
            }
            if ((a & 1023) == 0) {
                double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
                fprintf(stderr, "\rrow %d/%d  pairs=%llu  elapsed=%.1fs", a, m.M, total, el);
            }
        }
        fclose(f);
        double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        printf("\nCOUNT-DONE pairs=%llu emitted=%llu minc=%llu time=%.2fs rate=%.1f Mpair/s pruned_pairs=%llu\n",
               total, emitted, minc, el, total / el / 1e6, pruned_pairs);
        return 0;
    }
    fprintf(stderr, "usage: gbm_count --gen|--count|--selftest ...\n");
    return 2;
}
