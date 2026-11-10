// gbm_stats.cpp - Gene Bitmap Mining: pair-level statistics (E-B05) and E-A04 reference checks.
// Author: DAPR
//
// Modes:
//   gbm_stats --randtables K N seed out.csv        generate K random 2x2 tables (a b c d) with margins<=N
//   gbm_stats --fisherfile in.csv out.csv         two-sided Fisher exact p (R fisher.test definition) per row a,b,c,d
//   gbm_stats --pairs pairs.tsv out.tsv [fdr]      candidate table: OR + Fisher p + BH-FDR (over tested pairs)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
using namespace std;

// ---------- log-factorial table ----------
static vector<double> LF;
static void init_lf(int n) {
    LF.assign(n + 1, 0.0);
    for (int i = 1; i <= n; i++) LF[i] = LF[i - 1] + log((double)i);
}
static inline double logc(int n, int k) {
    if (k < 0 || k > n) return -INFINITY;
    return LF[n] - LF[k] - LF[n - k];
}
// hypergeometric log prob of table with fixed margins: P(X=k), k=c11
// margins: r1=c11+c10 (supA), c1=c11+c01 (supB), n=N
static inline double loghyper(int k, int r1, int c1, int n) {
    return logc(c1, k) + logc(n - c1, r1 - k) - logc(n, r1);
}
// two-sided Fisher p (R fisher.test): sum of hypergeom probs <= observed prob (with tolerance)
static double fisher_two_sided(int a, int b, int c, int d) {
    int r1 = a + b, r2 = c + d, c1 = a + c, n = a + b + c + d;
    if (n == 0) return 1.0;
    int lo = max(0, r1 - (n - c1));
    int hi = min(r1, c1);
    double obs = loghyper(a, r1, c1, n);
    double tol = 1e-7;
    double p = 0.0;
    for (int k = lo; k <= hi; k++) {
        double lp = loghyper(k, r1, c1, n);
        if (lp <= obs + tol) p += exp(lp);
    }
    return min(1.0, p);
}

// ---------- BH-FDR ----------
static vector<double> bh_fdr(vector<double> p) {
    size_t m = p.size();
    vector<size_t> idx(m);
    for (size_t i = 0; i < m; i++) idx[i] = i;
    sort(idx.begin(), idx.end(), [&](size_t x, size_t y) { return p[x] < p[y]; });
    vector<double> q(m);
    double prev = 1.0;
    for (size_t rank = m; rank >= 1; rank--) {
        size_t i = idx[rank - 1];
        double val = min(prev, p[i] * (double)m / (double)rank);
        q[i] = val;
        prev = val;
    }
    return q;
}

int main(int argc, char** argv) {
    if (argc >= 2 && strcmp(argv[1], "--randtables") == 0) {
        int K = atoi(argv[2]), N = atoi(argv[3]);
        uint64_t seed = strtoull(argv[4], 0, 10);
        FILE* f = fopen(argv[5], "w");
        fprintf(f, "a,b,c,d\n");
        mt19937_64 rng(seed);
        for (int i = 0; i < K; i++) {
            int r1 = (int)(rng() % (N + 1));
            int c1 = (int)(rng() % (N + 1));
            int lo = max(0, r1 + c1 - N), hi = min(r1, c1);
            int a = lo + (int)(rng() % (uint64_t)(hi - lo + 1));
            int b = r1 - a, c = c1 - a, d = N - r1 - c1 + a;
            fprintf(f, "%d,%d,%d,%d\n", a, b, c, d);
        }
        fclose(f);
        printf("RANDTABLES-OK K=%d -> %s\n", K, argv[5]);
        return 0;
    }
    if (argc >= 2 && strcmp(argv[1], "--fisherfile") == 0) {
        FILE* fi = fopen(argv[2], "r");
        FILE* fo = fopen(argv[3], "w");
        if (!fi || !fo) { fprintf(stderr, "io error\n"); return 1; }
        int a, b, c, d;
        init_lf(1 << 22);
        char line[256];
        fgets(line, sizeof line, fi); // header
        fprintf(fo, "a,b,c,d,p_fisher,odds_ratio\n");
        int nrows = 0;
        while (fscanf(fi, "%d,%d,%d,%d", &a, &b, &c, &d) == 4) {
            double p = fisher_two_sided(a, b, c, d);
            double or_ = (b * c == 0) ? (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))
                                      : (double)a * d / ((double)b * c);
            fprintf(fo, "%d,%d,%d,%d,%.12g,%.12g\n", a, b, c, d, p, or_);
            nrows++;
        }
        fclose(fi); fclose(fo);
        printf("FISHERFILE-OK rows=%d -> %s\n", nrows, argv[3]);
        return 0;
    }
    if (argc >= 2 && strcmp(argv[1], "--pairs") == 0) {
        const char* in = argv[2];
        const char* out = argv[3];
        FILE* fi = fopen(in, "r");
        if (!fi) { fprintf(stderr, "cannot open %s\n", in); return 1; }
        char line[4096];
        fgets(line, sizeof line, fi); // header
        struct Row { string A, B; int c11, c10, c01, c00; };
        vector<Row> rows;
        char A[512], B[512];
        int c11, supA, supB, c10, c01, c00;
        while (fscanf(fi, "%511s%511s%d%d%d%d%d%d", A, B, &c11, &supA, &supB, &c10, &c01, &c00) == 8) {
            rows.push_back({A, B, c11, c10, c01, c00});
        }
        fclose(fi);
        int maxn = 1;
        for (auto& r : rows) maxn = max(maxn, r.c11 + r.c10 + r.c01 + r.c00);
        init_lf(maxn + 8);
        auto t0 = chrono::steady_clock::now();
        vector<double> pvals(rows.size());
        FILE* fo = fopen(out, "w");
        fprintf(fo, "famA\tfamB\tc11\tc10\tc01\tc00\todds_ratio\tp_fisher\n");
        for (size_t i = 0; i < rows.size(); i++) {
            auto& r = rows[i];
            pvals[i] = fisher_two_sided(r.c11, r.c10, r.c01, r.c00);
            double or_ = (r.c10 * r.c01 == 0) ? (r.c11 + 0.5) * (r.c00 + 0.5) / ((r.c10 + 0.5) * (r.c01 + 0.5))
                                              : (double)r.c11 * r.c00 / ((double)r.c10 * r.c01);
            fprintf(fo, "%s\t%s\t%d\t%d\t%d\t%d\t%.6g\t%.6g\n",
                    r.A.c_str(), r.B.c_str(), r.c11, r.c10, r.c01, r.c00, or_, pvals[i]);
        }
        fclose(fo);
        double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
        auto padj = bh_fdr(pvals);  // screened-set BH (m = number of screened pairs)
        // conservative all-pairs BH (m = C(M,2)), optional M via argv[4]
        vector<double> padj_all = padj;
        bool has_all = argc >= 5;
        if (has_all) {
            long long Mfam = atoll(argv[4]);
            long double mall = (long double)Mfam * (Mfam - 1) / 2.0L;
            vector<size_t> idx2(rows.size());
            for (size_t i = 0; i < rows.size(); i++) idx2[i] = i;
            sort(idx2.begin(), idx2.end(), [&](size_t x, size_t y) { return pvals[x] < pvals[y]; });
            double prev = 1.0;
            for (size_t rank = rows.size(); rank >= 1; rank--) {
                size_t i = idx2[rank - 1];
                double val = min(prev, (double)((long double)pvals[i] * mall / (long double)rank));
                padj_all[i] = val;
                prev = val;
            }
        }
        vector<size_t> idx(rows.size());
        for (size_t i = 0; i < rows.size(); i++) idx[i] = i;
        sort(idx.begin(), idx.end(), [&](size_t x, size_t y) { return padj[x] < padj[y]; });
        FILE* ff = fopen(out, "w");
        if (has_all)
            fprintf(ff, "famA\tfamB\tc11\tc10\tc01\tc00\todds_ratio\tp_fisher\tpadj_BH_screened\tpadj_BH_allpairs\trank\n");
        else
            fprintf(ff, "famA\tfamB\tc11\tc10\tc01\tc00\todds_ratio\tp_fisher\tpadj_BH\trank\n");
        int rank = 0;
        for (size_t i : idx) {
            auto& r = rows[i];
            double or_ = (r.c10 * r.c01 == 0) ? (r.c11 + 0.5) * (r.c00 + 0.5) / ((r.c10 + 0.5) * (r.c01 + 0.5))
                                              : (double)r.c11 * r.c00 / ((double)r.c10 * r.c01);
            if (has_all)
                fprintf(ff, "%s\t%s\t%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%d\n",
                        r.A.c_str(), r.B.c_str(), r.c11, r.c10, r.c01, r.c00, or_, pvals[i], padj[i], padj_all[i], ++rank);
            else
                fprintf(ff, "%s\t%s\t%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\t%d\n",
                        r.A.c_str(), r.B.c_str(), r.c11, r.c10, r.c01, r.c00, or_, pvals[i], padj[i], ++rank);
        }
        fclose(ff);
        printf("PAIRS-STATS-OK tested=%zu time=%.2fs -> %s\n", rows.size(), el, out);
        return 0;
    }
    fprintf(stderr, "usage: gbm_stats --randtables|--fisherfile|--pairs ...\n");
    return 2;
}
