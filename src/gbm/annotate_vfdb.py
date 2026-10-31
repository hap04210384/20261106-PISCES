#!/usr/bin/env python
# annotate_vfdb.py - Map gene families to VFDB virulence factors via k-mer Jaccard.
# Data-prep/analysis step. Input: pangenome h5 (family protein seqs), VFDB_setA_pro.fas.gz
# Usage: python annotate_vfdb.py pangenome.h5 VFDB_setA_pro.fas.gz out.tsv [K]
import sys, gzip, re
import h5py

def kmers(seq, K):
    return set(seq[i:i+K] for i in range(len(seq)-K+1))

def read_fasta_gz(path):
    seqs={}
    name=None; buf=[]
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('>'):
                if name: seqs[name]=''.join(buf)
                name=line[1:].split()[0]; buf=[]
            else:
                buf.append(line.strip())
    if name: seqs[name]=''.join(buf)
    return seqs

def main():
    h5path, vfdb_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
    K = int(sys.argv[4]) if len(sys.argv)>4 else 7
    print('reading VFDB...')
    vfdb = read_fasta_gz(vfdb_path)
    vf_kmers = {n: kmers(s, K) for n,s in vfdb.items() if len(s)>=K}
    print(f'VFDB proteins: {len(vf_kmers)}')
    f = h5py.File(h5path,'r')
    gf = f['geneFamilies']
    fo = open(out,'w')
    fo.write('family\tvf_protein\tvf_header\tjaccard\tshared_kmers\n')
    n_with_seq=0; n_total=0
    # try to find protein sequences per family
    sample = list(gf.keys())[:5]
    print('sample family keys:', {k: list(gf[k].keys()) for k in sample[:2]})
    for fam in gf.keys():
        n_total+=1
        g = gf[fam]
        seq = None
        for cand in ('protein','sequence','seq','representative_protein'):
            if cand in g:
                try:
                    d = g[cand]
                    if d.shape and d.shape[0]>0:
                        v = d[0]
                        seq = v.decode() if isinstance(v,bytes) else str(v)
                except Exception:
                    pass
                if seq: break
        if not seq:
            # some layouts: 'sequences' group
            if 'sequences' in g:
                try:
                    s=g['sequences']
                    key=list(s.keys())[0]
                    v=s[key][0] if hasattr(s[key],'__len__') and s[key].shape else s[key][()]
                    seq = v.decode() if isinstance(v,bytes) else None
                except Exception:
                    seq=None
        if not seq:
            continue
        n_with_seq+=1
        fk = kmers(seq, K)
        if not fk: continue
        best=None; bestj=0.0
        for name,vk in vf_kmers.items():
            inter = len(fk & vk)
            if inter==0: continue
            j = inter/len(fk | vk)
            if j>bestj: bestj=j; best=(name,inter)
        if best and bestj>=0.3:
            fo.write(f'{fam}\t{best[0]}\t{bestj:.3f}\t{best[1]}\n')
        if n_with_seq % 5000 == 0:
            print(f'  annotated {n_with_seq} families', flush=True)
    fo.close()
    print(f'families total={n_total} with_seq={n_with_seq} -> {out}')

if __name__=='__main__':
    main()
