# Using PISCES on your own data

PISCES is not tied to the study cohorts: the engines consume a
**genome-by-feature presence/absence bitmap (GBM1)**, and one generic encoder
turns standard pangenome tool outputs — or any plain table — into that bitmap.
If you can produce a presence/absence matrix, you can screen every pair of
features exactly.

## One-command pipeline

From the repository root, after building the engines (see README, Quick start):

```
# Windows
pisces_run.bat <your_table.csv> my_run

# Linux / macOS
./pisces_run.sh <your_table.csv> my_run
```

This runs three deterministic stages and writes into `my_run/`:

| File | Content |
|---|---|
| `input.gbm1.txt` | the GBM1 bitmap (byte-identical for identical input) |
| `pairs.tsv` | exact co-occurrence / exclusion / conditional counts for all surviving pairs |
| `candidates.tsv` | odds ratio, two-sided Fisher exact p, BH-FDR (screened-set and conservative all-pairs), rank |

`candidates.tsv` is the multiplicity-controlled result: filter
`padj_BH_allpairs < 0.05` for the conservative call set used in the study.

Accepted inputs (auto-detected; force with `--style`, see
`src/gbm/to_gbm1.py --help`):

- **roary** `gene_presence_absence.csv`
- **Panaroo** `gene_presence_absence.Rtab`
- **PPanGGOLiN / PanGBank** pangenome HDF5 (via `h5_to_gbm1.py` /
  `pangbank_to_gbm1.py`)
- any **CSV/TSV matrix** (features x genomes, 0/1 cells)
- a **two-column long list** (`feature,genome` per line)

## Worked example: roary output to candidate pairs

```
# 1. run roary on your genome annotation GFFs (example: 50 E. coli genomes)
roary -f roary_out -e -n -v *.gff

# 2. encode (roary CSV auto-detected)
python src/gbm/to_gbm1.py roary_out/gene_presence_absence.csv my_run/input.gbm1.txt

# 3. count all pairs exactly (min support 10 here; pairs below it are pruned
#    by the anti-monotone bound and never materialized)
src/gbm/gbm_count --count my_run/input.gbm1.txt my_run/pairs.tsv 10

# 4. statistics — pass M = number of feature families (column 3 of the GBM1 header)
src/gbm/gbm_stats --pairs my_run/pairs.tsv my_run/candidates.tsv 5279
```

The same pipeline with `minc = 5` and a 2,000-genome cohort is the study's
release configuration; the bundled synthetic demo (`demo/run_demo.sh` /
`demo\run_demo.bat`) runs the identical chain on a 2,000 x 600 matrix with 18
planted interactions and self-checks recovery — a good first sanity run before
pointing the pipeline at your own data.

## Scope

Nothing in the engines is gene-specific. Any collection of binary features
profiled across a set of samples works: microbial gene families (this study),
GWAS presence/absence variants, protein domain architectures, k-mer content
binarized at a threshold, or survey checklists. The only requirement is a
feature-by-sample presence/absence table.

## Notes

- Encoding is a one-pass linear scan and deterministic: the same input table
  always yields a byte-identical bitmap, so results never depend on how the
  bitmap was produced.
- The count stage is exact by construction (no sketching, no sampling). The
  optional hypergeometric deviation screen (`gbm_count --count ... <minc> screen`)
  emits exactly the pairs passing |z| >= 3.29 and is byte-identical to the
  unscreened stream on those pairs.
- GPU engines (`gbm_router`, `gbm_pair_gpu`) are drop-in faster routes for the
  count stage; both emit byte-identical output to the AVX2 CPU reference.
