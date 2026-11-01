#!/usr/bin/env python
# pangbank_to_gbm1.py - Convert PanGBank pangenome HDF5 (flat layout, blosc2)
# to GBM1 bitmap text matrix. Data preparation only; experiments run in C++.
# Requires: pip install hdf5plugin (provides the blosc2 HDF5 filter).
# Usage: python pangbank_to_gbm1.py in.h5 out.txt
import sys
import numpy as np
import h5py
import hdf5plugin  # noqa: F401  (registers blosc2 filter on import)

def main():
    inp, out = sys.argv[1], sys.argv[2]
    f = h5py.File(inp, 'r')

    genomes = [n.decode() for n in f['annotations/genomes']['name']]
    gidx = {g: i for i, g in enumerate(genomes)}
    N = len(genomes)
    print(f'N={N} genomes')

    # contig id -> genome name
    ctg = f['annotations/contigs']
    ctg_genome = {int(cid): gen.decode() for cid, gen in zip(ctg['ID'], ctg['genome'])}

    genes = f['annotations/genes']
    g_ids = genes['ID']          # S21 gene IDs
    g_ctg = genes['contig']      # u4 contig ids
    print(f'{len(g_ids)} genes loaded')

    gf = f['geneFamilies']
    gf_gene = gf['gene']         # S21
    gf_fam = gf['geneFam']       # S21
    print(f'{len(gf_gene)} gene-family rows loaded')

    # family names -> codes
    fams_all = f['geneFamiliesInfo']['name']
    fam_list = [n.decode() for n in fams_all]
    fam_index = {fn: i for i, fn in enumerate(fam_list)}
    M = len(fam_list)
    print(f'M={M} families')

    # sort geneFamilies by gene id for join
    order = np.argsort(gf_gene, kind='stable')
    gf_gene_s = gf_gene[order]
    gf_fam_s = gf_fam[order]

    pos = np.searchsorted(gf_gene_s, g_ids)
    pos_clip = np.minimum(pos, len(gf_gene_s) - 1)
    hit = gf_gene_s[pos_clip] == g_ids
    fam_of_gene = gf_fam_s[pos_clip]
    print(f'genes with family: {int(hit.sum())}/{len(g_ids)}')

    # map gene -> (genome, family codes)
    gi = np.empty(len(g_ids), dtype=np.int64)
    missing_genome = 0
    for i, cid in enumerate(g_ctg):
        gname = ctg_genome.get(int(cid))
        if gname is None:
            gi[i] = -1
            missing_genome += 1
        else:
            gi[i] = gidx[gname]
    if missing_genome:
        print(f'WARN genes with unknown contig: {missing_genome}')

    # family codes via dict lookup on bytes (vectorized-ish)
    fam_bytes = {fn.encode(): i for i, fn in enumerate(fam_list)}
    fi = np.fromiter((fam_bytes.get(x, -1) for x in fam_of_gene), dtype=np.int64, count=len(g_ids))

    mat = np.zeros((M, (N + 63) // 64), dtype=np.uint64)
    ok = (gi >= 0) & (fi >= 0) & hit
    g_ok = gi[ok]
    f_ok = fi[ok]
    words = g_ok >> 6
    bits = g_ok & 63
    # scatter OR: use flat index
    flat = f_ok * mat.shape[1] + words
    uflat, inv = np.unique(flat, return_inverse=True)
    vals = np.zeros(len(uflat), dtype=np.uint64)
    np.bitwise_or.at(vals, inv, np.uint64(1) << bits.astype(np.uint64))
    mat.ravel()[uflat] = vals
    presences = int(ok.sum())
    print(f'total presences={presences}')

    # support distribution sanity
    bitsmat = np.unpackbits(mat.view(np.uint8), axis=1, bitorder='little')[:, :N]
    sup = bitsmat.sum(axis=1)
    print(f'support: min={sup.min()} max={sup.max()} mean={sup.mean():.1f} '
          f'median={np.median(sup):.0f} (expect heavy-tailed)')

    with open(out, 'w') as fo:
        fo.write(f'GBM1 {N} {M}\n')
        for fi_i, fname in enumerate(fam_list):
            row = bitsmat[fi_i]
            fo.write(fname + ' ')
            fo.write(''.join('1' if b else '0' for b in row))
            fo.write('\n')
    print('WROTE', out)

if __name__ == '__main__':
    main()
