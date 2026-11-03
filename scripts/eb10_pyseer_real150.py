#!/usr/bin/env python
# eb10_pyseer_real150.py - E-B10: re-baseline pyseer (FE + LMM) on the REAL
# KPN planted-150 benchmark (kpn_planted150, seed42) - same instance as the
# headline 84/150 (editor review A1). Protocol mirrors E-B07 dry-run + E-G03:
# per planted pair, both directions; variant = presence(famX) single-variant
# Rtab; phenotype = presence(famY); significance = lrt_p < Bonferroni
# 0.05/54590 = 9.16e-7. Resumable; concurrent.
import os, subprocess, time, threading
from concurrent.futures import ThreadPoolExecutor

MW = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PD = MW + r"\data\processed"
WD = PD + r"\eb10"
PY = MW + r"\code\envs\pyseer-env\Scripts\python.exe"
SIM = PD + r"\eg03\similarity.tsv"
OUT = MW + r"\code\results\tb_eb10_pyseer_real150.tsv"
BONF = 0.05 / 54590.0
os.makedirs(WD, exist_ok=True)
N = 1992

# ---- load only needed family bitstrings ----
truth = []
for line in open(PD + r"\kpn_planted150_truth.tsv"):
    p = line.rstrip("\n").split("\t")
    if p[0] == "type":
        continue
    truth.append((p[0], p[1], p[2]))
need = {f for _, a, b in truth for f in (a, b)}
bits = {}
with open(PD + r"\kpn_planted150.txt") as f:
    hdr = f.readline()
    for line in f:
        sp = line.index(" ")
        nm = line[:sp]
        if nm in need:
            bits[nm] = line[sp + 1:].strip()
missing = need - set(bits)
assert not missing, ("families missing from bitmap", list(missing)[:5])
print("loaded", len(bits), "family bitstrings", flush=True)

# ---- resumability ----
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        p = line.rstrip("\n").split("\t")
        if p and p[0] not in ("pair", "") and len(p) >= 4:
            done.add((int(p[0]), p[2], int(p[3])))
tasks = []
for i, (typ, a, b) in enumerate(truth):
    for mode in ("FE", "LMM"):
        for d, (var, phen) in enumerate(((a, b), (b, a))):
            if (i, mode, d) not in done:
                tasks.append((i, typ, mode, d, var, phen))
print("tasks to run:", len(tasks), flush=True)

lock = threading.Lock()
fh = open(OUT, "a")
if fh.tell() == 0:
    fh.write("pair\ttype\tmode\tdirection\tvariant\tphenotype\tlrt_p\tfilter_p\tcall\truntime_s\n")

def run_one(task):
    i, typ, mode, d, var, phen = task
    t0 = time.time()
    tag = WD + rf"\p{i:03d}_{mode}_{d}"
    with open(tag + ".Rtab", "w") as f:
        f.write("variant\t" + "\t".join(f"g{j}" for j in range(N)) + "\n")
        f.write(var + "\t" + "\t".join(bits[var]) + "\n")
    with open(tag + ".tsv", "w") as f:
        f.write("strain\tgroup\n")
        f.write("\n".join(f"g{j}\t{bits[phen][j]}" for j in range(N)))
    cmd = [PY, "-m", "pyseer", "--phenotypes", tag + ".tsv",
           "--pres", tag + ".Rtab", "--min-af", "0.0001", "--cpu", "1"]
    if mode == "FE":
        cmd += ["--no-distances"]
    else:
        cmd += ["--distances", SIM]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=200)
        lrt = filt = ""
        for ln in out.stdout.splitlines():
            parts = ln.split("\t")
            if len(parts) >= 8 and parts[0] == var:
                filt, lrt = parts[2], parts[3]
        try:
            lp = float(lrt)
            call = 1 if lp < BONF else 0
            lrt_s = "%.3e" % lp
        except (ValueError, TypeError):
            lrt_s, call = "", 0
        rc = out.returncode
    except subprocess.TimeoutExpired:
        lrt_s, filt, call, rc = "", "", 0, "TMO"
    row = (f"{i}\t{typ}\t{mode}\t{d}\t{var}\t{phen}\t{lrt_s}\t{filt}\t{call}"
           f"\t{time.time()-t0:.1f}\n")
    with lock:
        fh.write(row)
        fh.flush()
    return i, mode, d, call, lrt_s

t_start = time.time()
with ThreadPoolExecutor(max_workers=8) as ex:
    n_done = 0
    for res in ex.map(run_one, tasks):
        n_done += 1
        if n_done % 25 == 0:
            el = time.time() - t_start
            print(f"{n_done}/{len(tasks)} elapsed {el:.0f}s", flush=True)
print("ALL DONE", time.time() - t_start, flush=True)
