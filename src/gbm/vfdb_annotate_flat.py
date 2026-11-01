#!/usr/bin/env python
# vfdb_annotate_flat.py - Map PanGBank gene families (flat geneFamiliesInfo layout)
# to VFDB virulence factors via 7-mer Jaccard, using an inverted k-mer index.
# Usage: python vfdb_annotate_flat.py pangenome.h5 VFDB_setA_pro.fas.gz out.tsv [K] [threshold]
import sys, gzip
from collections import Counter, defaultdict
import hdf5plugin  # noqa: F401
import h5py

def read_fasta_gz(path):
    seqs = {}
    name = None; buf = []
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if name: seqs[name] = ''.join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name: seqs[name] = ''.join(buf)
    return seqs

def kmerset(seq, K):
    return set(seq[i:i+K] for i in range(len(seq) - K + 1))

def main():
    h5path, vfdb_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
    K = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    THR = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
    print('reading VFDB...')
    vfdb = read_fasta_gz(vfdb_path)
    vf_names = []
    vf_kmers = []
    index = defaultdict(list)  # kmer -> [vf idx]
    for n, s in vfdb.items():
        ks = kmerset(s, K)
        if not ks: continue
        idx = len(vf_names)
        vf_names.append(n)
        vf_kmers.append(ks)
        for k in ks:
            index[k].append(idx)
    print(f'VFDB proteins indexed: {len(vf_names)}, unique {K}-mers: {len(index)}')
    del vfdb

    f = h5py.File(h5path, 'r')
    info = f['geneFamiliesInfo']
    names = info['name']
    proteins = info['protein']
    fo = open(out, 'w')
    fo.write('family\tvf_protein\tjaccard\tshared_kmers\n')
    n_total = 0; n_seq = 0; n_hit = 0
    for i in range(len(names)):
        n_total += 1
        raw = proteins[i]
        if isinstance(raw, bytes):
            raw = raw.split(b'\x00')[0].decode('ascii', 'ignore')
        else:
            raw = str(raw).split('\x00')[0]
        seq = raw.strip()
        if len(seq) < K:
            continue
        n_seq += 1
        fk = kmerset(seq, K)
        cand = Counter()
        for k in fk:
            for vi in index.get(k, ()):
                cand[vi] += 1
        best = None; bestj = 0.0
        for vi, inter in cand.items():
            j = inter / (len(fk) + len(vf_kmers[vi]) - inter)
            if j > bestj:
                bestj = j; best = (vi, inter)
        if best and bestj >= THR:
            n_hit += 1
            fo.write(f'{names[i].decode()}\t{vf_names[best[0]]}\t{bestj:.3f}\t{best[1]}\n')
        if n_seq % 10000 == 0:
            print(f'  processed {n_seq} families, hits {n_hit}', flush=True)
    fo.close()
    print(f'total={n_total} with_seq={n_seq} hits={n_hit} -> {out}')

if __name__ == '__main__':
    main()
