#!/usr/bin/env python
# ec05_anytime.py - E-C05 anytime behaviour: (1) killed-run partial output must be
# a byte-exact prefix of the full deterministic output; (2) slice concat must equal
# the full output (kill-resume compositionality).
import subprocess, time, hashlib, os, sys

MW = r"E:\20261031-paper-Tensor core solves big problems\Manuscript_work"
EXE = MW + r"\code\src\gbm\gbm_count.exe"
IN = MW + r"\data\processed\eco615_gbm1.txt"
PD = MW + r"\data\processed"
FULL = PD + r"\ec05_full.tsv"

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 24), b""):
            h.update(ch)
    return h.hexdigest()

rows = []

# 0. reference full output md5 (produced in previous step) + GPU pairs file
full_md5 = md5(FULL)
ref_pairs_md5 = md5(PD + r"\eco615_pairs.tsv")
rows.append(("ref_full_vs_gpu_pairs", full_md5, ref_pairs_md5, full_md5 == ref_pairs_md5))
print("full==gpu pairs file:", full_md5 == ref_pairs_md5, flush=True)

# 1. kill test: interrupt a fresh run at ~20 s
killed = PD + r"\ec05_killed.tsv"
t0 = time.time()
p = subprocess.Popen([EXE, "--count", IN, killed, "2", "screen", "bound"],
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
time.sleep(20)
p.kill()
p.wait()
part = open(killed, "rb").read()
full_head = open(FULL, "rb").read(len(part))
prefix_ok = part == full_head
rows.append(("kill_at_20s_prefix", len(part), "bytes", prefix_ok))
print(f"killed run: {len(part)} bytes, prefix of full output: {prefix_ok}", flush=True)

# 2. slices [0,29387) + [29387,58774), concat == full
s1, s2 = PD + r"\ec05_slice0.tsv", PD + r"\ec05_slice1.tsv"
for a0, a1, out in [(0, 29387, s1), (29387, 58774, s2)]:
    t1 = time.time()
    r = subprocess.run([EXE, "--count", IN, out, "2", "screen", "bound", str(a0), str(a1)],
                       capture_output=True, text=True, timeout=290)
    tail = [l for l in r.stdout.splitlines() if "COUNT-DONE" in l]
    print(f"slice [{a0},{a1}): {tail[-1] if tail else r.stdout[-200:]} ({time.time()-t1:.0f}s)", flush=True)

cat = PD + r"\ec05_slicecat.tsv"
with open(cat, "wb") as f:
    f.write(open(s1, "rb").read())
    f.write(open(s2, "rb").read())
# NOTE: slice outputs each carry the header line; full output has exactly one.
# Compare after stripping per-slice headers (line 1).
def strip_header(p):
    d = open(p, "rb").read()
    return d[d.index(b"\n") + 1:]
s12 = strip_header(s1) + strip_header(s2)
full_body = strip_header(FULL)
concat_ok = s12 == full_body
rows.append(("slice_concat_vs_full_body", len(s12), len(full_body), concat_ok))
print("slice concat == full body:", concat_ok, flush=True)

# 3. idempotent resume: rerunning slice 0 after 'interruption' gives identical bytes
s1b = PD + r"\ec05_slice0_retry.tsv"
subprocess.run([EXE, "--count", IN, s1b, "2", "screen", "bound", "0", "29387"],
               capture_output=True, text=True, timeout=290)
resume_ok = open(s1, "rb").read() == open(s1b, "rb").read()
rows.append(("resume_idempotent", md5(s1), md5(s1b), resume_ok))
print("resume idempotent:", resume_ok, flush=True)

# write results
res = MW + r"\code\results\tb_ec05_anytime.tsv"
with open(res, "w") as f:
    f.write("check\ta\tb\tpassed\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
print("done ->", res, flush=True)
sys.stdout.flush()
