#!/usr/bin/env python
# ee04_overlap.py - E-E04 external-cohort replication: main KPN cohort vs pang774.
# Both sides significant at BH-FDR 0.05 over all pairs C(M,2) (the manuscript's
# headline multiplicity口径). The 774-side adjusted p is recomputed with the
# hybrid step-up formula p*C(M,2)/rank over the screened set (the file's
# padj_BH column is screened-set口径 and is intentionally not used).
# Cross-cohort family matching by identical locus-tag family names: a family is
# "matched" iff it occurs in a main-cohort significant pair AND its name exists
# in the 774 pangenome family list.
import time, os

MW = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PD = MW + r"\data\processed"
MAIN_SIG = MW + r"\tables\tb_eb05_candidates_sig_allpairs.tsv"
S774 = PD + r"\pang774_stats.tsv"
GBM_774 = PD + r"\pang774_gbm1.txt"
OUT = os.path.dirname(os.path.abspath(__file__)) + os.sep + ".." + os.sep + "results" + os.sep + "tb_ee04_external_replication.txt"

M1, M2 = 54590, 55185
C774 = M2 * (M2 - 1) // 2
print(f"main M = {M1} | 774 M = {M2} | C(774,2) = {C774:,}")

t0 = time.time()
fams774 = set()
with open(GBM_774, encoding="utf-8", errors="replace") as f:
    f.readline()
    for line in f:
        if line.strip():
            fams774.add(line.split(" ", 1)[0].strip())
print(f"774 families: {len(fams774):,} ({time.time()-t0:.0f}s)")

def canon(a, b):
    return (a, b) if a <= b else (b, a)

main_sig = {}
fams_in_sig = set()
with open(MAIN_SIG, encoding="utf-8", errors="replace") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 5:
            continue
        try:
            pv = float(p[2]) if p[2] else 0.0
        except ValueError:
            pv = 0.0
        a, b = p[0], p[1]
        main_sig[canon(a, b)] = pv
        if a in fams774:
            fams_in_sig.add(a)
        if b in fams774:
            fams_in_sig.add(b)
print(f"main BH 0.05 significant pairs: {len(main_sig):,} ({time.time()-t0:.0f}s)")

matched = fams_in_sig   # families in main sig pairs whose names exist in pang774
print(f"matched families: {len(matched):,}")

sig774 = set()
n774 = 0
with open(S774, encoding="utf-8", errors="replace") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 10:
            continue
        n774 += 1
        try:
            pv = float(p[7])
            rk = float(p[9])
        except ValueError:
            continue
        if rk > 0 and pv * C774 / rk < 0.05:
            sig774.add(canon(p[0], p[1]))
print(f"774 screened pairs: {n774:,} | 774 BH 0.05 significant: {len(sig774):,} ({time.time()-t0:.0f}s)")

main_matched = {pr for pr in main_sig if pr[0] in matched and pr[1] in matched}
sig774_matched = {pr for pr in sig774 if pr[0] in matched and pr[1] in matched}
rep = len(main_matched & sig774_matched)
union = len(main_matched | sig774_matched)
top10k = sorted(main_matched, key=lambda pr: main_sig[pr])[:10000]
rep10k = sum(1 for pr in top10k if pr in sig774)

res = [
    ("basis", "BH-FDR 0.05 all-pairs C(M,2) both sides; 774 padj hybrid p*C(M,2)/rank"),
    ("main_M", M1), ("main_sig_BH005", len(main_sig)),
    ("774_M", M2), ("774_sig_BH005", len(sig774)),
    ("matched_families", len(matched)),
    ("main_sig_matched", len(main_matched)),
    ("774_sig_matched", len(sig774_matched)),
    ("recall_main_sig_in_774", f"{rep}/{len(main_matched)} = {rep/len(main_matched):.4f}"),
    ("jaccard_sig_matched", f"{rep}/{union} = {rep/union:.4f}"),
    ("recall_top10k_by_p", f"{rep10k}/10000 = {rep10k/10000:.4f}"),
]
with open(OUT, "w", encoding="utf-8") as f:
    f.write("metric\tvalue\n")
    for k, v in res:
        f.write(f"{k}\t{v}\n")
for k, v in res:
    print(f"{k} = {v}")
print(f"done in {time.time()-t0:.0f}s -> {OUT}")
