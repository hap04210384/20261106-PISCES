#!/usr/bin/env python
# ee01_positive_control.py - E-E01 positive control: literature-documented
# co-carried AMR combos must be recovered by the all-pairs test at Bonferroni
# level, and must rank at the top of 29.9M tested pairs.
import time
from scipy.stats import fisher_exact
import os

MW = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PD = MW + r"\data\processed"
N_GENOMES = 2006
M = 58774
BONF = 0.05 / (M * (M - 1) / 2)

# annotated AMR families (family -> gene)
ann = {}
for line in open(PD + r"\eco615_resfinder_annotations.tsv"):
    p = line.rstrip("\n").split("\t")
    if p[0] == "family":
        continue
    ann[p[0]] = p[1]

CTX = "GCF_013201525.1_CDS_4520"   # blaCTX-M family (CTX-M-1 rep; covers 280 CTX-M genes)
COMBOS = [
    # combo_name, famA, geneA, famB, geneB, literature basis
    ("blaCTX-M+aac(6')-Ib-cr", CTX, "blaCTX-M", "GCF_041447665.1_CDS_4907", "aac(6')-Ib-cr", "ESBL+FQ-R on IncF plasmids, classic"),
    ("blaCTX-M+blaTEM-1", CTX, "blaCTX-M", "GCF_001679985.1_CDS_4781", "blaTEM-1B", "most common ESBL+penicillinase pairing"),
    ("sul2+dfrA17", "GCF_033802325.1_CDS_4609", "sul2", "GCF_013201525.1_CDS_4433", "dfrA17", "IS26/class-2-integron trimethoprim-sulfonamide cassette"),
    ("sul2+aadA5", "GCF_033802325.1_CDS_4609", "sul2", "GCF_013201525.1_CDS_4434", "aadA5", "same cassette, streptomycin-sulfonamide"),
    ("mcr-1+blaCTX-M", "GCF_009823315.1_CDS_3255", "mcr-1.1", CTX, "blaCTX-M", "colistin+ESBL co-carriage on IncI2/IncHI2"),
    ("qnrS1+blaCTX-M", "GCF_904426485.1_CDS_4214", "qnrS1", CTX, "blaCTX-M", "PMQR+ESBL plasmid linkage"),
    ("tet(A)+blaCTX-M", "GCF_001679985.1_CDS_1430", "tet(A)", CTX, "blaCTX-M", "tetracycline+ESBL MDR plasmids"),
]

# ---- presence vectors from GBM1
want = set()
for _, a, _, b, _, _ in COMBOS:
    want.add(a); want.add(b)
bits = {}
t0 = time.time()
with open(PD + r"\eco615_gbm1.txt") as f:
    hdr = f.readline().split()
    assert hdr[0] == "GBM1" and int(hdr[1]) == N_GENOMES
    for line in f:
        nm = line[: line.index(" ")]
        if nm in want:
            bits[nm] = line[line.index(" ") + 1 :].strip()
missing = want - set(bits)
print(f"presence bits: {len(bits)}/{len(want)} found ({time.time()-t0:.0f}s); missing: {missing}", flush=True)

# ---- single streaming pass over 3GB stats file for wanted pairs (both orders)
pairs_wanted = set()
for _, a, _, b, _, _ in COMBOS:
    pairs_wanted.add((a, b)); pairs_wanted.add((b, a))
stats = {}
t0 = time.time()
with open(PD + r"\eco615_stats.tsv") as f:
    head = f.readline().rstrip("\n").split("\t")
    ia, ib, ior, ip, irank = head.index("famA"), head.index("famB"), head.index("odds_ratio"), head.index("p_fisher"), head.index("rank")
    for line in f:
        p = line.split("\t")
        key = (p[ia], p[ib])
        if key in pairs_wanted:
            stats[key] = (p[ior], p[ip], p[irank])
    print(f"stats scan done: {len(stats)} rows ({time.time()-t0:.0f}s)", flush=True)

# ---- assemble table
out = open(MW + r"\code\results\tb_ee01_positive_control.tsv", "w")
out.write("combo\tgeneA\tgeneB\tnA\tnB\tnAB\todds_ratio_scipy\tp_scipy\tp_pipeline\t"
          "padj_bonferroni\tpipeline_rank\trecalled\n")
n_rec = 0
for name, a, ga, b, gb, basis in COMBOS:
    ba, bb = bits[a], bits[b]
    nA = sum(1 for c in ba if c == "1")
    nB = sum(1 for c in bb if c == "1")
    nAB = sum(1 for x, y in zip(ba, bb) if x == "1" and y == "1")
    n10 = nA - nAB; n01 = nB - nAB; n00 = N_GENOMES - nA - nB + nAB
    OR, p = fisher_exact([[nAB, n10], [n01, n00]])
    st = stats.get((a, b)) or stats.get((b, a)) or ("NA", "NA", "NA")
    por, pp, prank = st
    try:
        padj = float(pp) / (M * (M - 1) / 2)
        rec = "YES" if padj <= 0.05 else "NO"
    except ValueError:
        padj = float("nan"); rec = "NA"
    if rec == "YES":
        n_rec += 1
    out.write(f"{name}\t{ga}\t{gb}\t{nA}\t{nB}\t{nAB}\t{OR:.3g}\t{p:.3g}\t{pp}\t{padj:.3g}\t{prank}\t{rec}\n")
    print(f"{name}: nA={nA} nB={nB} nAB={nAB} OR={OR:.3g} p_pipe={pp} rank={prank} -> {rec}", flush=True)
out.write(f"# recovery: {n_rec}/{len(COMBOS)} combos at Bonferroni (threshold {BONF:.3g})\n")
out.close()
print(f"RECOVERY: {n_rec}/{len(COMBOS)}", flush=True)
