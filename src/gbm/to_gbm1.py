#!/usr/bin/env python
# to_gbm1.py - Generic presence/absence -> GBM1 bitmap encoder.
# Data preparation step (not an experiment); experiments run in C++ (gbm_count/gbm_stats).
#
# Accepts the de-facto standard outputs of pangenome tools, or any plain
# presence/absence matrix, and emits the GBM1 bitmap text consumed by the
# counting engines. Encoding is a deterministic one-pass linear scan; the same
# input always yields byte-identical GBM1 output.
#
# Formats (--style):
#   matrix     header row: <row-label>,genome1,genome2,...  body: one row per
#              feature; cell value truthy = present ('1','true','yes','present',
#              any non-zero number, any non-empty free text). CSV or TSV
#              (auto-detected from the extension / first line).
#   roary      roary gene_presence_absence.csv: the 15 fixed leading columns are
#              skipped; every non-empty cell under a genome column = present.
#   panaroo    Panaroo gene_presence_absence.Rtab: 'Gene' + one 0/1 column per
#              genome.
#   long       two columns, no header: <feature>,<genome> per line (presence
#              list). Duplicates collapse to a single 1-bit.
#   auto       default: panaroo if the first header cell is 'Gene' and the file
#              is TSV; roary if the header contains 'No. isolates'; long if the
#              file has exactly two columns and no header; otherwise matrix.
#
# Usage:
#   python to_gbm1.py IN OUT [--style auto] [--label-col 0]
#
# The feature (row) order and the genome (column) order of the input are
# preserved in the output; no reordering, no filtering, no re-coding.
import sys
import csv
import os

TRUTHY = {'1', 'true', 't', 'yes', 'y', 'present', 'x'}


def is_truthy(cell):
    c = cell.strip()
    if not c:
        return False
    if c.lower() in TRUTHY:
        return True
    try:
        return float(c) != 0.0
    except ValueError:
        return True  # free text (e.g. a gene name) counts as presence


def open_rows(path):
    """Yield raw rows, auto-detecting the delimiter (tab vs comma)."""
    with open(path, newline='', encoding='utf-8-sig') as f:
        sample = f.read(8192)
        f.seek(0)
        first = sample.splitlines()[0] if sample else ''
        delim = '\t' if ('\t' in first or path.lower().endswith('.rtab')) else ','
        for row in csv.reader(f, delimiter=delim):
            if row:
                yield row


def parse_matrix(path, label_col=0):
    rows = open_rows(path)
    header = next(rows)
    genomes = [h.strip() for i, h in enumerate(header) if i != label_col]
    fams, bits = [], []
    for row in rows:
        if label_col >= len(row):
            continue
        fams.append(row[label_col].strip())
        bits.append([is_truthy(c) for i, c in enumerate(row) if i != label_col])
    return genomes, fams, bits


def parse_panaroo(path):
    return parse_matrix(path, label_col=0)  # same layout as generic matrix


# roary gene_presence_absence.csv fixed leading columns (everything before the
# per-genome columns). Counted by matching header names, so extra/fewer columns
# in future roary versions are tolerated.
ROARY_LEAD_NAMES = {
    'gene', 'non-unique gene name', 'annotation', 'no. isolates',
    'no. sequences', 'avg sequences per isolate', 'genome fragment',
    'order within fragment', 'accessory fragment',
    'accessory order with fragment', 'qc', 'min group size nuc',
    'max group size nuc', 'avg group size nuc',
}


def roary_lead_cols(header):
    n = 0
    for h in header:
        if h.strip().lower() in ROARY_LEAD_NAMES:
            n += 1
        else:
            break
    if n == 0:
        n = 14  # classic roary layout fallback
    return n


def parse_roary(path):
    rows = open_rows(path)
    header = next(rows)
    lead = roary_lead_cols(header)
    genomes = [h.strip() for h in header[lead:]]
    fams, bits = [], []
    for row in rows:
        if not row:
            continue
        fams.append(row[0].strip())
        cells = row[lead:]
        cells += [''] * (len(genomes) - len(cells))
        bits.append([is_truthy(c) for c in cells])
    return genomes, fams, bits


def parse_long(path):
    genomes, fams = [], []
    seen = set()
    for row in open_rows(path):
        if len(row) < 2:
            continue
        fam, gen = row[0].strip(), row[1].strip()
        if (fam, gen) in seen:
            continue
        seen.add((fam, gen))
        fams.append(fam)
        genomes.append(gen)
    # long format carries no header, so input row order is meaningless; sort
    # both axes lexicographically for deterministic output
    g_order = sorted(dict.fromkeys(genomes))
    f_order = sorted(dict.fromkeys(fams))
    gidx = {g: i for i, g in enumerate(g_order)}
    bits = [[False] * len(g_order) for _ in f_order]
    fidx = {f: i for i, f in enumerate(f_order)}
    for (fam, gen) in seen:
        bits[fidx[fam]][gidx[gen]] = True
    return g_order, f_order, bits


def detect_style(path):
    rows = open_rows(path)
    header = next(rows, None)
    if header is None:
        sys.exit('ERROR: empty input')
    low = [h.strip().lower() for h in header]
    n_cols = len(header)
    if low and low[0] == 'gene' and ('\t' in open(path, encoding='utf-8-sig').readline()):
        return 'panaroo'
    if 'no. isolates' in low:
        return 'roary'
    if n_cols == 2 and low[0] not in ('feature', 'gene', 'family', 'genome'):
        return 'long'
    return 'matrix'


PARSERS = {'matrix': parse_matrix, 'panaroo': parse_panaroo,
           'roary': parse_roary, 'long': parse_long}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    style = 'auto'
    label_col = 0
    for i, a in enumerate(sys.argv[1:]):
        if a == '--style':
            style = sys.argv[i + 2]
        elif a == '--label-col':
            label_col = int(sys.argv[i + 2])
    inp, out = args[0], args[1]
    if style == 'auto':
        style = detect_style(inp)
    print(f'style={style}')
    if style == 'long':
        genomes, fams, bits = parse_long(inp)
    else:
        genomes, fams, bits = PARSERS[style](inp, label_col) if style == 'matrix' \
            else PARSERS[style](inp)
    N, M = len(genomes), len(fams)
    if N == 0 or M == 0:
        sys.exit('ERROR: parsed an empty matrix')
    n_pres = sum(sum(r) for r in bits)
    print(f'N={N} genomes, M={M} families, presences={n_pres}')
    # GBM1 files use CRLF line endings (matching the study matrices)
    with open(out, 'w', newline='\r\n') as fo:
        fo.write(f'GBM1 {N} {M}\n')
        for fam, row in zip(fams, bits):
            fo.write(fam + ' ')
            fo.write(''.join('1' if b else '0' for b in row))
            fo.write('\n')
    print('WROTE', out)


if __name__ == '__main__':
    main()
