#!/usr/bin/env python
# ee02_fpr.py - E-E02 negative control: genome-label permutation null.
# Independently permute the genome axis of every family column (margins preserved,
# within-genome co-occurrence destroyed), rerun the exact GPU count + Fisher stats,
# and check empirical FPR against nominal levels.
import subprocess, time, os
import numpy as np
import os

MW = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PD = MW + r"\data\processed"
EXE = MW + r"\code\src\gbm\gbm_pair_gpu.exe"
STATS = MW + r"\code\src\gbm\gbm_stats.exe"
SRC = PD + r"\eco615_gbm1.txt"
M = 58774
BONF = 0.05 / (M * (M - 1) / 2)

def permute(seed):
    out = PD + rf"\eco615_perm{seed}.txt"
    rng = np.random.RandomState(seed)
    t0 = time.time()
    with open(SRC) as f, open(out, "w") as g:
        g.write(f.readline())
        for line in f:
            sp = line.index(" ")
            nm, bits = line[:sp], line[sp + 1:].strip()
            g.write(nm + " " + "".join(np.array(list(bits))[rng.permutation(len(bits))]) + "\n")
    print(f"permuted seed={seed} -> {out} ({time.time()-t0:.0f}s)", flush=True)
    return out

def run_null(tag, gbm1):
    t0 = time.time()
    pairs = PD + rf"\eco615_{tag}_pairs.tsv"
    r = subprocess.run([EXE, "--count", gbm1, pairs, str(1 << 26)],
                       capture_output=True, text=True, timeout=280)
    print(tag, "count:", r.stdout.strip().splitlines()[-1] if r.stdout else r.stderr, flush=True)
    st = PD + rf"\eco615_{tag}_stats.tsv"
    r = subprocess.run([STATS, "--pairs", pairs, st, str(M)],
                       capture_output=True, text=True, timeout=280)
    print(tag, "stats:", r.stdout.strip()[-120:], flush=True)
    n = 0; bonf = 0; bh = 0; ball = 0; minp = 1.0
    with open(st) as f:
        head = f.readline().rstrip("\n").split("\t")
        ip = head.index("p_fisher")
        ibh = head.index("padj_BH_screened") if "padj_BH_screened" in head else head.index("padj_BH")
        iba = head.index("padj_BH_allpairs") if "padj_BH_allpairs" in head else ibh
        for line in f:
            p = line.rstrip("\n").split("\t")
            pv = float(p[ip]); n += 1
            if pv < minp: minp = pv
            if pv <= BONF: bonf += 1
            if float(p[ibh]) <= 0.05: bh += 1
            if float(p[iba]) <= 0.05: ball += 1
    print(f"{tag}: tested={n} min_p={minp:.3g} bonf={bonf} ({bonf/n:.2e}) bh05={bh} ({bh/n:.2e}) all05={ball} ({ball/n:.2e}) ({time.time()-t0:.0f}s)", flush=True)
    return (tag, n, minp, bonf, bonf / n, bh, bh / n, ball, ball / n)

rows = []
for seed in [20261102, 777]:
    g = permute(seed)
    rows.append(run_null(f"perm{seed}", g))
    os.remove(g)  # permuted matrix regenerable; keep pairs/stats small? pairs too
    os.remove(PD + rf"\eco615_perm{seed}_pairs.tsv")

# real-data reference counts
n = 0; bonf = 0; bh = 0; ball = 0
with open(PD + r"\eco615_stats.tsv") as f:
    head = f.readline().rstrip("\n").split("\t")
    ip, ibh, iba = head.index("p_fisher"), head.index("padj_BH"), head.index("padj_all")
    for line in f:
        p = line.rstrip("\n").split("\t")
        n += 1
        if float(p[ip]) <= BONF: bonf += 1
        if float(p[ibh]) <= 0.05: bh += 1
        if float(p[iba]) <= 0.05: ball += 1
rows.append(("real", n, "", bonf, bonf / n, bh, bh / n, ball, ball / n))
print(f"real: tested={n} bonf={bonf} bh05={bh} all05={ball}", flush=True)

res = MW + r"\code\results\tb_ee02_fpr_calibration.tsv"
with open(res, "w") as f:
    f.write("dataset\ttested\tmin_p\tbonf_count\tbonf_fpr\tbh05_count\tbh05_fdr\tallpairs05_count\tallpairs05_fdr\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
print("done ->", res, flush=True)
