#!/usr/bin/env python
# eb10_scoary_real150.py - E-B10 Scoary arm: re-baseline Scoary (structure-naive
# Fisher, --no_pairwise) on the REAL KPN planted-150 benchmark, same instance
# as the headline 84/150 and the pyseer FE/LMM rerun. Protocol mirrors the
# E-B07 dry-run: per planted pair, BOTH directions; trait = presence(partner
# family); gene of interest = the other family; criterion per direction =
# Bonferroni_p < 0.05 (Scoary's own Bonferroni column); a pair is called when
# BOTH directions are significant. Resumable + concurrent; time-boxed via
# --budget seconds so it can be re-invoked until ALL DONE.
import os, sys, csv, time, threading, subprocess
from concurrent.futures import ThreadPoolExecutor, FIRST_COMPLETED, wait

MW = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PD = MW + r"\data\processed"
WD = PD + r"\eb10"
SCOARY = MW + r"\code\envs\pyseer-env\Scripts\scoary.exe"
GENES = WD + r"\gene_presence_absence.csv"
OUT = MW + r"\code\results\tb_eb10_scoary_real150.tsv"
BUDGET = float(sys.argv[sys.argv.index("--budget") + 1]) if "--budget" in sys.argv else 130.0
WORKERS = 12

truth = []
for line in open(PD + r"\kpn_planted150_truth.tsv"):
    p = line.rstrip("\n").split("\t")
    if p[0] == "type":
        continue
    truth.append((p[0], p[1], p[2]))

done = set()
if os.path.exists(OUT):
    with open(OUT, newline='', encoding='utf-8') as f:
        for p in csv.reader(f, delimiter='\t'):
            if p and p[0] not in ("pair", "") and len(p) >= 3:
                done.add((int(p[0]), int(p[1])))
tasks = [(i, d) for i in range(len(truth)) for d in (0, 1)
         if (i, d) not in done]
print("tasks to run:", len(tasks), flush=True)
if not tasks:
    print("ALL DONE", flush=True)
    sys.exit(0)

lock = threading.Lock()
fh = open(OUT, "a", newline="")
if fh.tell() == 0:
    fh.write("pair\tdirection\ttype\tvariant\ttrait\tbonferroni_p\tcalled\truntime_s\n")

def run_one(task):
    i, d = task
    typ, a, b = truth[i]
    # trait files were generated with phen=(a,b)[d]; gene of interest is the
    # partner family (mirrors eb10_pyseer_real150.py direction convention)
    phen, var = (a, b) if d == 0 else (b, a)
    od = WD + rf"\out_p{i:03d}_{d}"
    os.makedirs(od, exist_ok=True)
    t0 = time.time()
    logf = open(od + r"\run.log", "w")
    rc = subprocess.run(
        [SCOARY, "-t", WD + rf"\trait_p{i:03d}_{d}.csv",
         "-g", GENES, "-o", od, "--no-time", "--no_pairwise",
         "--p_value_cutoff", "1", "--max_hits", "60000"],
        stdout=logf, stderr=subprocess.STDOUT, timeout=900).returncode
    logf.close()
    bonf, called = "", 0
    res = os.path.join(od, phen + ".results.csv")
    if os.path.exists(res):
        with open(res, newline='', encoding='utf-8') as f:
            for row in csv.reader(f):
                if row and row[0] == var:
                    try:
                        bonf = float(row[11])
                        called = 1 if bonf < 0.05 else 0
                    except (ValueError, IndexError):
                        pass
                    break
    row = (f"{i}\t{d}\t{typ}\t{var}\t{phen}\t{bonf}\t{called}"
           f"\t{time.time()-t0:.1f}\n")
    with lock:
        fh.write(row)
        fh.flush()
    return i, d, called, rc

t_start = time.time()
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    pending = {}
    ti = iter(tasks)
    for _ in range(min(WORKERS, len(tasks))):
        try:
            t = next(ti)
            pending[ex.submit(run_one, t)] = t
        except StopIteration:
            break
    n_done = 0
    while pending:
        if time.time() - t_start > BUDGET:
            print("BUDGET EXHAUSTED - re-invoke to resume", flush=True)
            for fut in pending:
                fut.cancel()
            break
        finished, _ = wait(pending, timeout=5,
                           return_when=FIRST_COMPLETED)
        for fut in finished:
            pending.pop(fut)
            n_done += 1
            try:
                i, d, called, rc = fut.result()
                print(f"done pair={i} dir={d} called={called} rc={rc} "
                      f"({n_done} this session)", flush=True)
            except Exception as e:
                print("task error:", e, flush=True)
            try:
                t = next(ti)
                pending[ex.submit(run_one, t)] = t
            except StopIteration:
                pass
print("session complete", time.time() - t_start, flush=True)
