#!/usr/bin/env python
# merge_chunks.py - merge ranged chunks into the final h5 and verify size.
import os, sys
wd = sys.argv[1] if len(sys.argv)>1 else '.'
total = 1120068900
n = 16
csize = (total + n - 1)//n
parts = []
for i in range(n):
    if i == 8:
        for j in range(4):
            parts.append(os.path.join(wd, f'chunk_08_{j}.bin'))
        continue
    parts.append(os.path.join(wd, f'chunk_{i:02d}.bin'))
ok = True
for p in parts:
    sz = os.path.getsize(p)
    exp = None
    base = os.path.basename(p)
    if base.startswith('chunk_08_'):
        j = int(base.split('_')[2].split('.')[0])
        s0 = 8*csize; sc = 17501077
        exp = min(sc, 9*csize-1 - (s0+j*sc) + 1)
    else:
        i = int(base.split('_')[1].split('.')[0])
        exp = min(csize, total - i*csize)
    status = 'OK' if sz == exp else f'BAD (have {sz}, want {exp})'
    if sz != exp: ok = False
    print(f'{base}: {sz} {status}')
if not ok:
    print('CHUNKS-INCOMPLETE'); sys.exit(1)
out = os.path.join(wd, 'klebsiella_pneumoniae_refseq.h5')
with open(out, 'wb') as fo:
    for p in parts:
        with open(p, 'rb') as fi:
            while True:
                b = fi.read(1 << 24)
                if not b: break
                fo.write(b)
sz = os.path.getsize(out)
print(f'MERGED {out} size={sz} expected={total} {"MATCH" if sz==total else "MISMATCH"}')
sys.exit(0 if sz == total else 1)
