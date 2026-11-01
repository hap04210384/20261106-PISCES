#!/usr/bin/env python
# eco_prep2.py - family -> rep gene map, preferring non-fragment members.
import h5py, hdf5plugin, time
import numpy as np

MW = r"E:\20261031-paper-Tensor core solves big problems\Manuscript_work"
H5 = MW + r"\data\raw\pangbank\eco615.h5"
OUT = MW + r"\data\processed\eco615_fam_rep.tsv"

t0 = time.time()
f = h5py.File(H5, 'r')
gf = f['geneFamilies']
g = gf['gene'][:]
fm = gf['geneFam'][:]
print(f"gf loaded {len(g)} ({time.time()-t0:.0f}s)", flush=True)

gn = f['annotations/genes']
gID = gn['ID'][:]
gfrag = gn['is_fragment'][:]
print(f"genes loaded ({time.time()-t0:.0f}s)", flush=True)

if bool(np.all(gID[1:] >= gID[:-1])):
    sID, sfrag = gID, gfrag
    print("genes.ID already sorted", flush=True)
else:
    perm = np.argsort(gID, kind='stable')
    sID, sfrag = gID[perm], gfrag[perm]
    print(f"genes.ID sorted via argsort ({time.time()-t0:.0f}s)", flush=True)

pos = np.searchsorted(sID, g)
assert bool(np.all(sID[pos] == g)), "gene IDs not found in genes table"
frag = sfrag[pos]
print(f"frag flags mapped ({time.time()-t0:.0f}s)", flush=True)

order = np.lexsort((np.arange(len(g)), fm, frag))
fs = fm[order]
uniq_fam, first_idx = np.unique(fs, return_index=True)
sel = order[first_idx]
rep_genes = g[sel]
print(f"families: {len(uniq_fam)}, non-fragment reps: {int((frag[sel]==0).sum())} ({time.time()-t0:.0f}s)", flush=True)

gs = f['annotations/geneSequences']
gs_gene = gs['gene'][:]
gs_seqid = gs['seqid'][:]
pos2 = np.searchsorted(gs_gene, rep_genes)
assert bool(np.all(gs_gene[pos2] == rep_genes))
seqids = gs_seqid[pos2]

with open(OUT, 'w') as out:
    out.write("family\trep_gene\tseqid\n")
    for fam, rg, sid in zip(uniq_fam, rep_genes, seqids):
        out.write(f"{fam.decode()}\t{rg.decode()}\t{int(sid)}\n")
print(f"done -> {OUT} ({time.time()-t0:.0f}s)", flush=True)
