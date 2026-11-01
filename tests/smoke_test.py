#!/usr/bin/env python3
# smoke_test.py - end-to-end smoke test for the GBM engines.
# Usage: python tests/smoke_test.py [path-to-exe-dir]   (default: <repo>/src/gbm)
# Verifies, on a small synthetic matrix generated on the fly:
#   1. gbm_count --gen  -> deterministic matrix (regenerated hash matches)
#   2. gbm_count --count -> ALL pairs' c11/c10/c01/c00 match an independent
#      Python recount bit-for-bit; two repeated runs are byte-identical (md5)
#   3. gbm_stats --fisherfile -> two-sided Fisher p matches brute-force
#      hypergeometric enumeration (R-compatible definition)
#   4. gbm_sim --plant -> planted co-occurrence pairs recover c11 == target_support
# Exit code 0 = all green. No third-party dependencies (stdlib only).
import hashlib, math, os, subprocess, sys, tempfile

EXEDIR = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else \
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'gbm')
COUNT = os.path.join(EXEDIR, 'gbm_count.exe' if os.name == 'nt' else 'gbm_count')
STATS = os.path.join(EXEDIR, 'gbm_stats.exe' if os.name == 'nt' else 'gbm_stats')
SIM = os.path.join(EXEDIR, 'gbm_sim.exe' if os.name == 'nt' else 'gbm_sim')
N, M, SEED = 300, 800, 42


def run(exe, args):
    r = subprocess.run([exe] + args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{os.path.basename(exe)} {' '.join(args)} failed:\n{r.stdout}\n{r.stderr}")
    return r


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load_gbm1(path):
    with open(path) as f:
        magic, n, m = f.readline().split()
        assert magic == 'GBM1', magic
        names, bits = [], []
        for line in f:
            nm, b = line.split(' ', 1)
            names.append(nm)
            bits.append(int(b.strip(), 2))
        assert len(names) == int(m), (len(names), m)
    return names, bits


def fisher_two_sided_r(a, b, c, d):
    """R fisher.test two-sided: sum probabilities of tables with hypergeom
    probability <= probability of the observed table."""
    n = a + b + c + d
    r1, c1 = a + b, a + c
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)

    def logpmf(x):
        # log [ C(c1,x) * C(n-c1, r1-x) / C(n,r1) ]
        return (math.lgamma(c1 + 1) - math.lgamma(x + 1) - math.lgamma(c1 - x + 1)
                + math.lgamma(n - c1 + 1) - math.lgamma(r1 - x + 1)
                - math.lgamma(n - c1 - r1 + x + 1)
                - math.lgamma(n + 1) + math.lgamma(r1 + 1) + math.lgamma(n - r1 + 1))
    p_obs = logpmf(a)
    return min(1.0, sum(math.exp(logpmf(x)) for x in range(lo, hi + 1)
                        if logpmf(x) <= p_obs + 1e-12))


def main():
    for exe in (COUNT, STATS, SIM):
        assert os.path.isfile(exe), f"missing executable: {exe} (build the engines first)"
    tmp = tempfile.mkdtemp(prefix='gbm_smoke_')
    print(f"[smoke] exe dir: {EXEDIR}")
    print(f"[smoke] work dir: {tmp}")

    # 1. --gen determinism
    m1, m2 = os.path.join(tmp, 'm1.txt'), os.path.join(tmp, 'm2.txt')
    run(COUNT, ['--gen', str(N), str(M), str(SEED), m1])
    run(COUNT, ['--gen', str(N), str(M), str(SEED), m2])
    h1, h2 = md5(m1), md5(m2)
    assert h1 == h2, f"--gen not deterministic: {h1} vs {h2}"
    print(f"[smoke] 1/4 --gen deterministic (md5 {h1[:16]}...)")

    # 2. --count vs independent recount + repeat determinism
    names, bits = load_gbm1(m1)
    sup = [b.bit_count() for b in bits]
    truth = {}
    for i in range(M):
        bi = bits[i]
        for j in range(i + 1, M):
            c11 = (bi & bits[j]).bit_count()
            if c11 >= 3:
                truth[(names[i], names[j])] = (c11, sup[i], sup[j],
                                               sup[i] - c11, sup[j] - c11,
                                               N - sup[i] - sup[j] + c11)
    out1, out2 = os.path.join(tmp, 'p1.tsv'), os.path.join(tmp, 'p2.tsv')
    run(COUNT, ['--count', m1, out1, '3'])
    run(COUNT, ['--count', m1, out2, '3'])
    assert md5(out1) == md5(out2), "repeated --count runs differ"
    got = {}
    with open(out1) as f:
        f.readline()
        for line in f:
            t = line.split('\t')
            got[(t[0], t[1])] = tuple(map(int, t[2:8]))
    assert got == truth, (f"count mismatch: {len(got)} rows vs {len(truth)} expected; "
                          f"diff sample: {list(set(got.items()) ^ set(truth.items()))[:3]}")
    print(f"[smoke] 2/4 --count exact vs independent recount "
          f"({len(got)} pairs, minc=3), repeated runs byte-identical")

    # 3. --fisherfile vs brute-force enumeration
    tables = os.path.join(tmp, 'tables.csv')
    with open(tables, 'w') as f:
        f.write('a,b,c,d\n')
        for t in [(1, 9, 11, 3), (0, 10, 10, 0), (7, 3, 2, 8), (50, 20, 10, 120)]:
            f.write(','.join(map(str, t)) + '\n')
    fp = os.path.join(tmp, 'fisher.csv')
    run(STATS, ['--fisherfile', tables, fp])
    with open(fp) as f:
        f.readline()
        for line in f:
            a, b, c, d, p_eng, _ = line.strip().split(',')
            p_ref = fisher_two_sided_r(int(a), int(b), int(c), int(d))
            assert abs(float(p_eng) - p_ref) <= 1e-9 * max(1.0, p_ref), \
                f"fisher mismatch table=({a},{b},{c},{d}): engine {p_eng} vs ref {p_ref}"
    print("[smoke] 3/4 --fisherfile matches brute-force two-sided Fisher (4 tables)")

    # 4. --plant recoverability
    planted = os.path.join(tmp, 'planted.txt')
    truthf = os.path.join(tmp, 'truth.tsv')
    run(SIM, ['--plant', m1, planted, truthf, '5', '5', '5', '7'])
    outp = os.path.join(tmp, 'planted_pairs.tsv')
    run(COUNT, ['--count', planted, outp, '2', 'dep'])  # dep adds c11==0 depletion pairs
    gotp = {}
    with open(outp) as f:
        f.readline()
        for line in f:
            t = line.split('\t')
            gotp[(t[0], t[1])] = int(t[2])
    # direct matrix recount of all planted pairs (independent of engine output)
    pnames, pbits = load_gbm1(planted)
    pidx = {nm: i for i, nm in enumerate(pnames)}
    supp = [b.bit_count() for b in pbits]
    checked_cooc = checked_mutex = 0
    with open(truthf) as f:
        f.readline()
        for line in f:
            typ, a, b, ts, _gs = line.rstrip('\n').split('\t')
            ia, ib = pidx[a], pidx[b]
            c11_direct = (pbits[ia] & pbits[ib]).bit_count()
            if typ == 'mutex':
                assert c11_direct == 0, f"planted mutex {a},{b}: direct c11={c11_direct} != 0"
            else:
                assert c11_direct == int(ts), \
                    f"planted {typ} {a},{b}: direct c11={c11_direct} != target {ts}"
            if typ == 'cooc':
                assert gotp.get((a, b), gotp.get((b, a))) == int(ts), \
                    f"planted cooc {a},{b} missing from count output"
                checked_cooc += 1
            elif typ == 'mutex':
                # eligible depletion pairs (supA*supB/N >= minc) must be reported
                # with c11==0; with N=300 log-uniform supports are often too rare
                if min(supp[ia], supp[ib]) >= 2 and supp[ia] * supp[ib] // N >= 2:
                    assert gotp.get((a, b), gotp.get((b, a))) == 0, \
                        f"eligible planted mutex {a},{b} not recovered"
                    checked_mutex += 1
    assert checked_cooc == 5, f"cooc checked={checked_cooc}"
    # hand-crafted depletion case: N=100, two disjoint families with support 50
    # each (expected c11 = 25) MUST appear under the dep rule with c11 == 0;
    # an ultra-rare disjoint pair (support 2 each) must NOT appear (expected < minc)
    hm = os.path.join(tmp, 'hand.txt')
    with open(hm, 'w') as f:
        f.write('GBM1 100 2\n')
        f.write('famX\t' + '1' * 50 + '0' * 50 + '\n')
        f.write('famY\t' + '0' * 50 + '1' * 50 + '\n')
        f.write('famZ\t' + '1' * 1 + '0' * 98 + '1' * 1 + '\n')
        f.write('famW\t' + '0' * 1 + '1' * 98 + '0' * 1 + '\n')
    ho = os.path.join(tmp, 'hand_pairs.tsv')
    run(COUNT, ['--count', hm, ho, '5', 'dep'])
    with open(ho) as f:
        f.readline()
        hrows = {tuple(l.split('\t')[:2]): int(l.split('\t')[2]) for l in f}
    assert hrows.get(('famX', 'famY')) == 0, "eligible depletion pair famX/famY missing"
    assert ('famZ', 'famW') not in hrows, "ultra-rare depletion pair must be excluded"
    print(f"[smoke] 4/4 --plant: {checked_cooc} cooc recovered (c11==target), "
          f"all 15 direct c11 verified; dep rule verified on hand-crafted matrix")

    print("[smoke] ALL GREEN")


if __name__ == '__main__':
    main()
