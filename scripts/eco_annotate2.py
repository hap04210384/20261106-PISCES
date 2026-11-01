#!/usr/bin/env python
# eco_annotate2.py - chunked resumable ResFinder annotation of eco615 families.
# usage: python eco_annotate2.py <start> <end>
import sys, time
import numpy as np
import h5py, hdf5plugin

K = 11
START, END = int(sys.argv[1]), int(sys.argv[2])
MW = r"E:\20261031-paper-Tensor core solves big problems\Manuscript_work"
H5 = MW + r"\data\raw\pangbank\eco615.h5"
RF = MW + r"\data\raw\resfinder\rfdb\all.fsa"
MAP = MW + r"\data\processed\eco615_fam_rep.tsv"
OUT = MW + r"\data\processed\eco615_resfinder_annotations.tsv"

t0 = time.time()
RC = str.maketrans('ACGT', 'TGCA')
def canon(k):
    rc = k.translate(RC)[::-1]
    return k if k <= rc else rc

# ---- resfinder inverted index (rebuilt per run; all.fsa is small)
idx = {}
rf_names = []
rf_kmers = []
name = None; buf = []
def flush():
    global name, buf
    if name is None: return
    seq = ''.join(buf).upper()
    km = set()
    for i in range(len(seq) - K + 1):
        k = seq[i:i+K]
        if 'N' not in k:
            km.add(canon(k))
    gid = len(rf_names)
    rf_names.append(name)
    rf_kmers.append(km)
    for k in km:
        idx.setdefault(k, []).append(gid)
    name = None; buf = []
with open(RF) as fh:
    for line in fh:
        if line.startswith('>'):
            flush()
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
flush()
print(f"index: {len(rf_names)} genes, {len(idx)} kmers ({time.time()-t0:.0f}s)", flush=True)

# ---- family map
fams = []
with open(MAP) as fh:
    next(fh)
    for line in fh:
        p = line.rstrip('\n').split('\t')
        fams.append((p[0], int(p[2])))
print(f"families: {len(fams)}, chunk [{START},{END})", flush=True)

# ---- seqid -> row
f = h5py.File(H5, 'r')
sq = f['annotations/sequences']
sids = sq['seqid'][:]
row_of = {int(s): i for i, s in enumerate(sids)}

# ---- resume: skip families already annotated in OUT
done = set()
try:
    with open(OUT) as fh:
        next(fh)
        for line in fh:
            done.add(line.split('\t', 1)[0])
except FileNotFoundError:
    pass
new_file = len(done) == 0 and START == 0
out = open(OUT, 'a' if not new_file else 'w')
if new_file:
    out.write("family\trf_gene\tjaccard\tshared_kmers\n")
out.flush()

n_hit = 0; n_done = 0
dna_ds = sq['dna']
BATCH = 512
i = START
while i < END:
    j = min(i + BATCH, END)
    rows = []
    for k in range(i, j):
        sid = fams[k][1]
        r = row_of.get(sid)
        rows.append((k, r))
    sel = [r for _, r in rows if r is not None]
    dnas = dna_ds[sel] if sel else []
    di = 0
    for k, r in rows:
        fam = fams[k][0]
        seq = None
        if r is not None:
            seq = dnas[di].decode().rstrip('\x00').upper(); di += 1
        if fam in done:
            continue
        if not seq:
            continue
        fk = set()
        for x in range(len(seq) - K + 1):
            kk = seq[x:x+K]
            if 'N' not in kk:
                fk.add(canon(kk))
        if not fk:
            continue
        cand = {}
        for kk in fk:
            for gid in idx.get(kk, ()):
                cand[gid] = cand.get(gid, 0) + 1
        best, bestj, bestsh = -1, 0.0, 0
        for gid, sh in cand.items():
            uni = len(fk) + len(rf_kmers[gid]) - sh
            jj = sh / uni if uni else 0.0
            if jj > bestj:
                bestj, best, bestsh = jj, gid, sh
        if best >= 0 and bestj >= 0.5 and bestsh >= 30:
            out.write(f"{fam}\t{rf_names[best]}\t{bestj:.4f}\t{bestsh}\n")
            n_hit += 1
        n_done += 1
    out.flush()
    print(f"  [{i}:{j}] done, {n_hit} new hits ({time.time()-t0:.0f}s)", flush=True)
    i = j
out.close()
print(f"chunk complete: {n_done} processed, {n_hit} hits ({time.time()-t0:.0f}s)", flush=True)
