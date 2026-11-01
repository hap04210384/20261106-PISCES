#!/usr/bin/env python
# ee01_pyseer_baseline.py - reference GWAS baseline: pyseer FE single-variant,
# phenotype = presence of AMR gene A, variant set = all 110 annotated AMR families;
# record pyseer's p for the literature-documented partner B.
import os, subprocess, time

MW = r"E:\20261031-paper-Tensor core solves big problems\Manuscript_work"
PD = MW + r"\data\processed"
WD = PD + r"\ee01"
PY = MW + r"\code\envs\pyseer-env\Scripts\python.exe"
os.makedirs(WD, exist_ok=True)
N = 2006

CTX = "GCF_013201525.1_CDS_4520"
COMBOS = [
    ("blaCTX-M+aac(6')-Ib-cr", CTX, "aac(6')-Ib-cr", "GCF_041447665.1_CDS_4907"),
    ("blaCTX-M+blaTEM-1", CTX, "blaTEM-1B", "GCF_001679985.1_CDS_4781"),
    ("sul2+dfrA17", "GCF_033802325.1_CDS_4609", "dfrA17", "GCF_013201525.1_CDS_4433"),
    ("sul2+aadA5", "GCF_033802325.1_CDS_4609", "aadA5", "GCF_013201525.1_CDS_4434"),
    ("mcr-1+blaCTX-M", "GCF_009823315.1_CDS_3255", "blaCTX-M", CTX),
    ("qnrS1+blaCTX-M", "GCF_904426485.1_CDS_4214", "blaCTX-M", CTX),
    ("tet(A)+blaCTX-M", "GCF_001679985.1_CDS_1430", "blaCTX-M", CTX),
]

ann = []  # (fam, gene) preserving file order
for line in open(PD + r"\eco615_resfinder_annotations.tsv"):
    p = line.rstrip("\n").split("\t")
    if p[0] == "family":
        continue
    ann.append((p[0], p[1]))
want = {f for f, _ in ann}
bits = {}
with open(PD + r"\eco615_gbm1.txt") as f:
    f.readline()
    for line in f:
        nm = line[: line.index(" ")]
        if nm in want:
            bits[nm] = line[line.index(" ") + 1 :].strip()
assert len(bits) == len(ann), (len(bits), len(ann))

# Rtab: all 110 annotated AMR families as variants
rtab = WD + r"\amr110.Rtab"
strains = [f"g{j}" for j in range(N)]
with open(rtab, "w") as f:
    f.write("variant\t" + "\t".join(strains) + "\n")
    for fam, gene in ann:
        f.write(f"{fam}\t" + "\t".join(bits[fam]) + "\n")
print("rtab written", flush=True)

rows = []
for name, a, gB, b in COMBOS:
    for phen_fam, phen_gene, partner_fam, partner_gene, direction in [
        (a, name.split("+")[0], b, gB, "A->B"),
        (b, gB, a, name.split("+")[0], "B->A"),
    ]:
        t0 = time.time()
        trait = WD + rf"\trait_{len(rows):03d}.tsv"
        with open(trait, "w") as f:
            f.write("strain\tgroup\n")
            f.write("\n".join(f"g{j}\t{bits[phen_fam][j]}" for j in range(N)))
        out = subprocess.run(
            [PY, "-m", "pyseer", "--phenotypes", trait, "--pres", rtab,
             "--no-distances", "--min-af", "0.0001", "--cpu", "8"],
            capture_output=True, text=True, timeout=280)
        lrt = af = fp = beta = ""
        for line2 in out.stdout.splitlines():
            parts = line2.split("\t")
            if len(parts) >= 8 and parts[0] == partner_fam:
                af, fp, lrt, beta = parts[1], parts[2], parts[3], parts[4]
        rows.append((name, direction, phen_gene, partner_gene,
                     bits[phen_fam].count("1"), lrt, fp, beta, f"{time.time()-t0:.1f}s"))
        print(f"{name} {direction}: partner={partner_gene} lrt_p={lrt} ({rows[-1][-1]})", flush=True)

res = MW + r"\code\results\tb_ee01_pyseer_baseline.tsv"
with open(res, "w") as f:
    f.write("combo\tdirection\tphenotype_gene\tpartner_gene\tn_pheno\tpyseer_lrt_p\tpyseer_filter_p\tpyseer_beta\truntime\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
print("done ->", res, flush=True)
