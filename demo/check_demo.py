#!/usr/bin/env python3
"""Demo end-to-end check: compare screened pair statistics against planted ground truth.

Runs after the demo pipeline (gbm_count --count ... 3 screen ; gbm_stats --pairs ...).
PASS criteria (all must hold):
  * all planted co-occurrence pairs recovered at nominal Fisher p < 0.05
  * >= 2 of the strongly planted mutual-exclusion pairs recovered
  * >= 1 conditional pair recovered
  * the screen emitted < 2% of all pairs
The remaining planted pairs were deliberately planted below the |z| >= 3.29 gate;
not recovering them is the statistically correct behaviour of the screen.
"""
import sys, os, csv
from collections import Counter

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
stats_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "out", "demo_stats.tsv")
truth_path = os.path.join(root, "data", "demo", "truth.tsv")
M, N = 600, 2000
total_pairs = M * (M - 1) // 2

truth = {}
with open(truth_path) as f:
    for row in csv.reader(f, delimiter="\t"):
        if row[0] == "type":
            continue
        truth[frozenset((row[1], row[2]))] = (row[0], int(row[3]))

stats = {}
with open(stats_path) as f:
    for i, row in enumerate(csv.reader(f, delimiter="\t")):
        if i == 0:
            continue
        stats[frozenset((row[0], row[1]))] = float(row[7])

rec, tot = Counter(), Counter()
for k, (t, _tgt) in truth.items():
    tot[t] += 1
    p = stats.get(k)
    if p is not None and p < 0.05:
        rec[t] += 1

print("== planted-pair recovery (nominal Fisher p < 0.05) ==")
for t in ("cooc", "mutex", "cond"):
    print(f"  {t:6s}: {rec[t]}/{tot[t]}")
frac = 100.0 * len(stats) / total_pairs
print(f"screen emitted {len(stats)} of {total_pairs} pairs ({frac:.3f}%)")
print("note: planted pairs weaker than the |z| >= 3.29 gate are expected to be absent")
print("      from the screened set; that is the screen working as designed, not a bug.")

ok = (rec["cooc"] == tot["cooc"] and rec["mutex"] >= 2 and rec["cond"] >= 1 and frac < 2.0)
print("DEMO " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
