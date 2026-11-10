// gbm_approx.cpp - E-B09 approximate baselines: MinHash-LSH prefilter vs exact.
// Author: DAPR
//
// Modes:
//   gbm_approx --lsh in.txt out_pairs.tsv K B R
//       MinHash sketch (K hashes) + LSH banding (B bands x R rows) over family bitvectors.
//       Emits candidate pairs (union of band-bucket pairs), each with exact c11 computed
//       and the MinHash Jaccard estimate.
//   gbm_approx --compare lsh_pairs.tsv exact_pairs.tsv out.tsv
//       Precision/recall of the LSH candidate set vs the exact tested set (E-B05 amended rule).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#if defined(_MSC_VER)
#define POPCNT64(x) __popcnt64(x)
#else
#define POPCNT64(x) __builtin_popcountll(x)
#endif
using namespace std;

static inline uint64_t sm64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

struct Mat { int N, M, words; vector<uint64_t> bits; vector<string> names; vector<uint64_t> sup; };
static const uint64_t* rowp(const Mat& m, int j) { return m.bits.data() + (size_t)j * m.words; }
static bool load(Mat& m, const char* path) {
    FILE* f = fopen(path, "r"); if (!f) return false;
    char magic[8], name[256];
    if (fscanf(f, "%7s %d %d", magic, &m.N, &m.M) != 3 || strcmp(magic, "GBM1")) return false;
    m.words = (m.N + 63) / 64; m.bits.assign((size_t)m.M * m.words, 0ULL); m.names.resize(m.M); m.sup.assign(m.M, 0);
    for (int j = 0; j < m.M; j++) {
        if (fscanf(f, "%255s", name) != 1) return false;
        m.names[j] = name;
        int c, i = 0;
        while ((c = fgetc(f)) == ' ' || c == '\t') {}
        while (c == '0' || c == '1') { if (c == '1' && i < m.N) { m.bits[(size_t)j * m.words + (i >> 6)] |= (1ULL << (i & 63)); m.sup[j]++; } i++; c = fgetc(f); }
        if (i != m.N) { fprintf(stderr, "row %d bits %d != %d\n", j, i, m.N); return false; }
    }
    fclose(f); return true;
}

int main(int argc, char** argv) {
    if (argc >= 2 && strcmp(argv[1], "--lsh") == 0) {
        const char* in = argv[2]; const char* out = argv[3];
        int K = atoi(argv[4]), B = atoi(argv[5]), R = atoi(argv[6]);
        if (K != B * R) { fprintf(stderr, "K must equal B*R\n"); return 1; }
        Mat m; if (!load(m, in)) { fprintf(stderr, "load failed\n"); return 1; }
        printf("loaded N=%d M=%d, building %d-minhash sketches...\n", m.N, m.M, K);
        auto t0 = chrono::steady_clock::now();
        // sketch[f][k] = min over members g of sm64(g ^ seed_k)
        vector<vector<uint64_t>> sketch(m.M, vector<uint64_t>(K, UINT64_MAX));
        for (int f = 0; f < m.M; f++) {
            const uint64_t* r = rowp(m, f);
            for (int i = 0; i < m.N; i++) {
                if ((r[i >> 6] >> (i & 63)) & 1) {
                    for (int k = 0; k < K; k++) {
                        uint64_t h = sm64(((uint64_t)i << 32) ^ (uint64_t)k);
                        if (h < sketch[f][k]) sketch[f][k] = h;
                    }
                }
            }
        printf("\rsketch %d/%d", f, m.M); }
        printf("\nsketches done in %.1fs; banding %dx%d...\n",
               chrono::duration<double>(chrono::steady_clock::now() - t0).count(), B, R);
        // banding: bucket key = hash of band rows; collect family lists
        unordered_map<uint64_t, vector<int>> buckets;
        size_t total_bucket_refs = 0;
        for (int f = 0; f < m.M; f++) {
            for (int b = 0; b < B; b++) {
                uint64_t key = 1469598103934665603ULL ^ (uint64_t)b;
                for (int rr = 0; rr < R; rr++) {
                    uint64_t v = sketch[f][b * R + rr];
                    key = (key ^ v) * 1099511628211ULL;
                }
                buckets[key].push_back(f);
                total_bucket_refs++;
            }
        }
        printf("buckets: %zu, refs: %zu\n", buckets.size(), total_bucket_refs);
        // candidate pairs from buckets
        unordered_map<uint64_t, pair<int,int>> cand;
        for (auto& kv : buckets) {
            auto& v = kv.second;
            for (size_t i = 0; i < v.size(); i++)
                for (size_t j = i + 1; j < v.size(); j++) {
                    int a = min(v[i], v[j]), b = max(v[i], v[j]);
                    uint64_t key = ((uint64_t)a << 32) | (uint64_t)b;
                    cand.emplace(key, make_pair(a, b));
                }
        }
        printf("LSH candidate pairs: %zu\n", cand.size());
        // exact count + Jaccard estimate for each candidate
        FILE* fo = fopen(out, "w");
        fprintf(fo, "famA\tfamB\tc11_exact\tjac_est\tc11_est\n");
        for (auto& kv : cand) {
            int a = kv.second.first, b = kv.second.second;
            uint64_t c11 = 0;
            const uint64_t* ra = rowp(m, a); const uint64_t* rb = rowp(m, b);
            for (int w = 0; w < m.words; w++) c11 += POPCNT64(ra[w] & rb[w]);
            int inter = 0;
            for (int k = 0; k < K; k++) if (sketch[a][k] == sketch[b][k]) inter++;
            double J = (double)inter / K;
            double est = (m.sup[a] + m.sup[b]) > 0 ? J * (double)(m.sup[a] + m.sup[b]) / (1.0 + J) : 0.0;
            fprintf(fo, "%s\t%s\t%llu\t%.4f\t%.2f\n", m.names[a].c_str(), m.names[b].c_str(), c11, J, est);
        }
        fclose(fo);
        printf("LSH-DONE -> %s (%.1fs)\n", out, chrono::duration<double>(chrono::steady_clock::now() - t0).count());
        return 0;
    }
    if (argc >= 2 && strcmp(argv[1], "--compare") == 0) {
        // exact tested set from gbm_count output (first two columns)
        unordered_map<string, pair<int,int>> exact;
        FILE* fe = fopen(argv[3], "r");
        if (!fe) { fprintf(stderr, "cannot open exact %s\n", argv[3]); return 1; }
        char A[512], B[512]; int c11, sa, sb, c10, c01, c00;
        char line[2048];
        fgets(line, sizeof line, fe);
        while (fscanf(fe, "%511s%511s%d%d%d%d%d%d", A, B, &c11, &sa, &sb, &c10, &c01, &c00) == 8) {
            exact[string(A) + "\t" + string(B)] = {c11, 0};
        }
        fclose(fe);
        FILE* fl = fopen(argv[2], "r");
        if (!fl) { fprintf(stderr, "cannot open lsh %s\n", argv[2]); return 1; }
        fgets(line, sizeof line, fl);
        double j; double est;
        long long tp = 0, lsh_n = 0;
        FILE* fo = fopen(argv[4], "w");
        fprintf(fo, "famA\tfamB\tc11_exact\tin_exact_set\n");
        while (fscanf(fl, "%511s%511s%d%lf%lf", A, B, &c11, &j, &est) == 5) {
            lsh_n++;
            bool in = exact.count(string(A) + "\t" + string(B)) > 0;
            if (in) tp++;
            fprintf(fo, "%s\t%s\t%d\t%d\n", A, B, c11, in ? 1 : 0);
        }
        fclose(fl); fclose(fo);
        printf("COMPARE-DONE lsh=%lld exact_tested=%zu tp=%lld recall=%.4f precision=%.4f -> %s\n",
               lsh_n, exact.size(), tp, (double)tp / (double)exact.size(),
               lsh_n ? (double)tp / (double)lsh_n : 0.0, argv[4]);
        return 0;
    }
    fprintf(stderr, "usage: gbm_approx --lsh|--compare ...\n");
    return 2;
}
