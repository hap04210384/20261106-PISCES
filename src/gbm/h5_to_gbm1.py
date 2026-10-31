#!/usr/bin/env python
# h5_to_gbm1.py - Convert PanGBank/PPanggolin pangenome HDF5 to GBM1 bitmap text matrix.
# Data preparation step (not an experiment); experiments run in C++ (gbm_count/gbm_stats).
# Usage: python h5_to_gbm1.py in.h5 out.txt [--inspect]
import sys
import numpy as np
import h5py

def inspect(f):
    def walk(name, obj):
        kind = 'Group' if isinstance(obj, h5py.Group) else 'Dataset'
        shape = '' if isinstance(obj, h5py.Group) else str(obj.shape)
        print(f'{kind:8s} {name} {shape}')
    f.visititems(walk)

def main():
    inp, out = sys.argv[1], sys.argv[2]
    f = h5py.File(inp, 'r')
    if '--inspect' in sys.argv:
        inspect(f)
        # also print root-level attrs
        for k in f.attrs:
            print('ATTR', k, f.attrs[k])
        print('root keys:', list(f.keys()))
        return

    # PPanggolin layout: geneFamilies/<family_id> contains datasets 'genes' (gene IDs)
    # and annotations/genomes/<genome_name> contains per-genome gene info.
    # Strategy: read gene->genome mapping from annotations, family->genes from geneFamilies.
    gf = f['geneFamilies']
    fam_names = list(gf.keys())
    # annotations: genome -> genes dataset; figure out columns
    ann = f['annotations']
    genomes = list(ann.keys())
    # inspect one genome's genes dataset
    g0 = ann[genomes[0]]
    gkey = 'genes' if 'genes' in g0 else list(g0.keys())[0]
    dset = g0[gkey]
    print('genome genes dataset:', dset.shape, dset.dtype, 'keys of genome grp:', list(g0.keys()))
    # find column indices
    cols = {name: i for i, name in enumerate(dset.dtype.names)} if dset.dtype.names else {}
    print('columns:', cols)
    # gene id column and family column guesses
    # PPanggolin 'genes' structured array typically has fields: dna, name, start, stop, strand, type, family, ...
    fam_field = None
    for cand in ('family', 'gene_family', 'geneFamily'):
        if cand in cols:
            fam_field = cand
            break
    name_field = 'name' if 'name' in cols else ('dna' if 'dna' in cols else None)
    if fam_field is None:
        print('ERROR: cannot find family column; available:', cols)
        sys.exit(1)

    fam_index = {fn: i for i, fn in enumerate(fam_names)}
    M = len(fam_names)
    N = len(genomes)
    print(f'N={N} genomes, M={M} families')

    # write GBM1: header + one line per family would need transposed accumulation.
    # families x genomes matrix in memory: M * N bits = N/8 bytes per family
    mat = np.zeros((M, (N + 63) // 64), dtype=np.uint64)
    word_bits = 64
    for gi, gname in enumerate(genomes):
        genes = f['annotations'][gname][gkey][:]
        fams = genes[fam_field]
        # family values are bytes like b'G_01234' or family object reference; normalize
        for fam in fams:
            if isinstance(fam, bytes):
                fam = fam.decode()
            idx = fam_index.get(fam)
            if idx is not None:
                mat[idx, gi >> 6] |= np.uint64(1) << np.uint64(gi & 63)
        if gi % 200 == 0:
            print(f'  genome {gi}/{N}', flush=True)

    # write text GBM1
    with open(out, 'w') as fo:
        fo.write(f'GBM1 {N} {M}\n')
        bits = np.unpackbits(mat.view(np.uint8), axis=1, bitorder='little')
        # bits shape (M, words*64) -> slice N
        bits = bits[:, :N]
        for fi, fname in enumerate(fam_names):
            row = bits[fi]
            fo.write(f'{fname} ')
            fo.write(''.join('1' if b else '0' for b in row))
            fo.write('\n')
    print('WROTE', out)

if __name__ == '__main__':
    main()
