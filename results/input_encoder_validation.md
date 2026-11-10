# Input-encoder validation (2026-11-09)

Validates `src/gbm/to_gbm1.py`, the generic presence/absence -> GBM1 encoder.

## Round-trip at study scale

Path: the E. coli study matrix (`data/processed/eco615_gbm1.txt`, GBM1, 2,006 genomes
x 58,774 gene families, 9,055,117 presences, 119 MB)
-> exported as a plain CSV matrix (features x genomes, 0/1 cells)
-> re-encoded with `to_gbm1.py --style matrix`
-> GBM1 output is **byte-identical** to the original
(SHA-256 32fb10138c249c130babbda32ffbba4c61a2bc2508664cf657b1ea98b663daf6 for both).

## Round-trip at demo scale

Path: `data/demo/demo_matrix.txt` (GBM1, 2,000 genomes x 600 families)
-> exported as a plain CSV matrix (features x genomes, 0/1 cells)
-> re-encoded with `to_gbm1.py --style matrix`
-> GBM1 output is **byte-identical** to the original
(SHA-256 8cf738b388d35f17b385993ad29ad8fb12dac57d4161b3d85e90b22589f9b272 for both).

Running the reference pipeline (`gbm_count --count ... 3 screen` + `gbm_stats --pairs`)
on the re-encoded bitmap reproduces `demo/out/demo_stats.tsv` **byte-for-byte**;
`demo/check_demo.py` reports DEMO PASS (6/6 co-occurrence, 3/8 mutual-exclusion,
1/4 conditional planted pairs recovered; weak planted pairs below the |z| >= 3.29
gate correctly absent).

## Format adapters (synthetic fixtures, parsed bits verified against ground truth)

- Panaroo `gene_presence_absence.Rtab`: auto-detected, PASS
- roary `gene_presence_absence.csv`: auto-detected (lead columns matched by header
  name), PASS
- two-column long list: auto-detected, PASS (axes sorted lexicographically for
  deterministic output)
- plain CSV/TSV matrix: auto-detected, PASS

Encoding is a deterministic one-pass linear scan: the same input table always
produces byte-identical GBM1, independent of the toolchain that produced the table.
