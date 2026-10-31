// gbm_sim.cpp - Simulation: E-B08 ground-truth planting + E-A05 permutation-null verification.
// Author: YalongZhang
//
// Modes:
//   gbm_sim --plant in.txt out.txt truth.tsv K_cooc K_mutex K_cond seed
//       Plant interactions into a GBM1 matrix. Ground truth TSV columns:
//       type famA famB support_target genome_group
//       (conditional pairs are planted only inside a random genome subgroup of size N/2)
//   gbm_sim --permnull in.txt famA famB [L]
//       Observed two-sided Fisher p for the pair + empirical permutation null over L random
//       label permutations (uniform shuffle of genome labels). Prints observed p, empirical
//       permutation p, and (if N<=10) exact enumeration over all N! permutations.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
using namespace std;

static vector<double> LF;
static void init_lf(int n) { LF.assign(n + 1, 0.0); for (int i = 1; i <= n; i++) LF[i] = LF[i - 1] + log((double)i); }
static inline double logc(int n, int k) { if (k < 0 || k > n) return -INFINITY; return LF[n] - LF[k] - LF[n - k]; }
static inline double loghyper(int k, int r1, int c1, int n) { return logc(c1, k) + logc(n - c1, r1 - k) - logc(n, r1); }
static double fisher_two_sided(int a, int b, int c, int d) {
    int r1 = a + b, c1 = a + c, n = a + b + c + d;
    if (n == 0) return 1.0;
    int lo = max(0, r1 - (n - c1)), hi = min(r1, c1);
    double obs = loghyper(a, r1, c1, n), tol = 1e-7, p = 0.0;
    for (int k = lo; k <= hi; k++) { double lp = loghyper(k, r1, c1, n); if (lp <= obs + tol) p += exp(lp); }
    return min(1.0, p);
}

struct Mat { int N, M, words; vector<uint64_t> bits; vector<string> names; };
static uint64_t* row(Mat& m, int j) { return m.bits.data() + (size_t)j * m.words; }
static bool load(Mat& m, const char* path) {
    FILE* f = fopen(path, "r"); if (!f) return false;
    char magic[8], name[256];
    if (fscanf(f, "%7s %d %d", magic, &m.N, &m.M) != 3 || strcmp(magic, "GBM1")) return false;
    m.words = (m.N + 63) / 64; m.bits.assign((size_t)m.M * m.words, 0ULL); m.names.resize(m.M);
    for (int j = 0; j < m.M; j++) {
        if (fscanf(f, "%255s", name) != 1) return false;
        m.names[j] = name;
        int c, i = 0;
        while ((c = fgetc(f)) == ' ' || c == '\t') {}
        while (c == '0' || c == '1') { if (c == '1' && i < m.N) row(m, j)[i >> 6] |= (1ULL << (i & 63)); i++; c = fgetc(f); }
        if (i != m.N) { fprintf(stderr, "row %d bits %d != %d\n", j, i, m.N); return false; }
    }
    fclose(f); return true;
}
static bool save(Mat& m, const char* path) {
    FILE* f = fopen(path, "w"); if (!f) return false;
    fprintf(f, "GBM1 %d %d\n", m.N, m.M);
    for (int j = 0; j < m.M; j++) {
        fprintf(f, "%s ", m.names[j].c_str());
        for (int i = 0; i < m.N; i++) fputc(((row(m, j)[i >> 6] >> (i & 63)) & 1) ? '1' : '0', f);
        fputc('\n', f);
    }
    fclose(f); return true;
}
static uint64_t support(Mat& m, int j) {
    uint64_t s = 0; for (int w = 0; w < m.words; w++) s += (uint64_t)__popcnt64(row(m, j)[w]);
    return s;
}

int main(int argc, char** argv) {
    if (argc >= 2 && strcmp(argv[1], "--plant") == 0) {
        const char* in = argv[2]; const char* out = argv[3]; const char* truth = argv[4];
        int Kcooc = atoi(argv[5]), Kmutex = atoi(argv[6]), Kcond = atoi(argv[7]);
        uint64_t seed = strtoull(argv[8], 0, 10);
        Mat m; if (!load(m, in)) { fprintf(stderr, "load failed\n"); return 1; }
        mt19937_64 rng(seed);
        FILE* ft = fopen(truth, "w");
        fprintf(ft, "type\tfamA\tfamB\ttarget_support\tgroup_size\n");
        // unique family pool: each planted pair uses two families never used before
        vector<int> pool(m.M);
        for (int i = 0; i < m.M; i++) pool[i] = i;
        shuffle(pool.begin(), pool.end(), rng);
        size_t pool_pos = 0;
        auto pick_two = [&](int& a, int& b) {
            a = pool[pool_pos++]; b = pool[pool_pos++];
        };
        auto clear_row = [&](int j) { fill(row(m, j), row(m, j) + m.words, 0ULL); };
        auto plant_group = [&](int j, vector<int>& gs) { for (int g : gs) row(m, j)[g >> 6] |= (1ULL << (g & 63)); };
        // supports log-uniform 1e-3..1e-1 of N
        auto rand_support = [&]() {
            double lo = log10(max(2.0, m.N * 1e-3)), hi = log10(m.N * 1e-1);
            double v = lo + (double)rng() / (double)UINT64_MAX * (hi - lo);
            int s = (int)llround(pow(10.0, v));
            return max(2, min(s, m.N / 2));
        };
        vector<int> all(m.N); for (int i = 0; i < m.N; i++) all[i] = i;
        int planted = 0;
        for (int t = 0; t < Kcooc && planted < m.M / 2; t++, planted++) {
            int a, b; pick_two(a, b); int s = rand_support();
            shuffle(all.begin(), all.end(), rng);
            vector<int> gs(all.begin(), all.begin() + s);
            clear_row(a); clear_row(b);
            plant_group(a, gs); plant_group(b, gs);
            fprintf(ft, "cooc\t%s\t%s\t%d\t%d\n", m.names[a].c_str(), m.names[b].c_str(), s, m.N);
        }
        for (int t = 0; t < Kmutex && planted < m.M / 2; t++, planted++) {
            int a, b; pick_two(a, b); int s = rand_support();
            shuffle(all.begin(), all.end(), rng);
            vector<int> ga(all.begin(), all.begin() + s);
            vector<int> gb(all.begin() + s, all.begin() + 2 * s > all.end() ? all.end() : all.begin() + 2 * s);
            clear_row(a); clear_row(b);
            plant_group(a, ga); plant_group(b, gb);
            fprintf(ft, "mutex\t%s\t%s\t%d\t%d\n", m.names[a].c_str(), m.names[b].c_str(), s, m.N);
        }
        for (int t = 0; t < Kcond && planted < m.M / 2; t++, planted++) {
            int a, b; pick_two(a, b); int s = rand_support();
            int gsize = m.N / 2;
            shuffle(all.begin(), all.end(), rng);
            vector<int> group(all.begin(), all.begin() + gsize);
            shuffle(group.begin(), group.end(), rng);
            vector<int> gs(group.begin(), group.begin() + min(s, gsize));
            clear_row(a); clear_row(b);
            plant_group(a, group);            // a: present in subgroup
            plant_group(b, gs);             // b: co-occurs with a only inside subgroup
            fprintf(ft, "cond\t%s\t%s\t%d\t%d\n", m.names[a].c_str(), m.names[b].c_str(), s, gsize);
        }
        fclose(ft);
        if (!save(m, out)) { fprintf(stderr, "save failed\n"); return 1; }
        printf("PLANT-OK planted=%d -> %s (truth %s)\n", planted, out, truth);
        return 0;
    }
    if (argc >= 2 && strcmp(argv[1], "--permnull") == 0) {
        const char* in = argv[2]; string fa = argv[3], fb = argv[4];
        int L = (argc >= 6) ? atoi(argv[5]) : 10000;
        Mat m; if (!load(m, in)) { fprintf(stderr, "load failed\n"); return 1; }
        int ia = -1, ib = -1;
        for (int j = 0; j < m.M; j++) { if (m.names[j] == fa) ia = j; if (m.names[j] == fb) ib = j; }
        if (ia < 0 || ib < 0) { fprintf(stderr, "families not found\n"); return 1; }
        init_lf(m.N + 8);
        auto stat_of = [&](vector<int>& perm) {
            // permuted labels: c11 = overlap of row ia with row ib under permuted genome indices
            int c11 = 0; int sa = 0, sb = 0;
            for (int i = 0; i < m.N; i++) {
                int pi = perm[i];
                int a = (int)((row(m, ia)[pi >> 6] >> (pi & 63)) & 1);
                int b = (int)((row(m, ib)[i >> 6] >> (i & 63)) & 1);
                c11 += a & b; sa += a; sb += b;
            }
            int c10 = sa - c11, c01 = sb - c11, c00 = m.N - c11 - c10 - c01;
            return fisher_two_sided(c11, c10, c01, c00);
        };
        vector<int> id(m.N); for (int i = 0; i < m.N; i++) id[i] = i;
        double p_obs = stat_of(id);
        mt19937_64 rng(20261031);
        int exceed = 0;
        for (int l = 0; l < L; l++) {
            vector<int> p = id; shuffle(p.begin(), p.end(), rng);
            if (stat_of(p) <= p_obs) exceed++;
        }
        double p_emp = (exceed + 1) / (double)(L + 1);
        printf("PERMNULL obs_p=%.6g emp_p(L=%d)=%.6g\n", p_obs, L, p_emp);
        if (m.N <= 10) {
            // exhaustive enumeration of all permutations
            vector<int> p = id;
            long long tot = 0, exc = 0;
            sort(p.begin(), p.end());
            do {
                tot++;
                if (stat_of(p) <= p_obs) exc++;
            } while (next_permutation(p.begin(), p.end()));
            printf("PERMNULL-EXACT total_perm=%lld exceed=%lld exact_p=%.6g\n", tot, exc, (exc + 1) / (double)(tot + 1));
        }
        return 0;
    }
    fprintf(stderr, "usage: gbm_sim --plant|--permnull ...\n");
    return 2;
}
