# Gene Bitmap Mining (GBM) reproducibility code

Exact pair-level co-occurrence mining on genome x gene-family presence/absence bitmaps.
Author: YalongZhang (ylzhang@qzc.edu.cn)

## Contents

- `src/gbm/gbm_count.cpp` — all-pairs exact co-occurrence counting (AVX2, Release).
  Modes: `--gen`, `--count` (with optional `dep` depletion rule), `--countpruned`, `--selftest`.
- `src/gbm/gbm_stats.cpp` — pair-level statistics: two-sided Fisher exact (R-compatible
  definition), odds ratio, BH-FDR. Modes: `--randtables`, `--fisherfile`, `--pairs`.
- `src/gbm/gbm_sim.cpp` — ground-truth interaction planting (co-occurrence / mutual
  exclusion / conditional) and permutation-null verification. Modes: `--plant`, `--permnull`.
- `src/gbm/gbm_approx.cpp` — MinHash-LSH approximate prefilter baseline (E-group comparison).
  Modes: `--lsh`, `--compare`.
- `src/gbm/h5_to_gbm1.py` — converts a PPanGGOLiN/PanGBank pangenome HDF5 file into the
  GBM1 text bitmap matrix (data preparation only; all experiments run in the C++ engines).
- `src/gbm/annotate_vfdb.py` — maps gene families to VFDB virulence factors (k-mer Jaccard).
- `src/gbm/*.bat` — MSVC Release build scripts (`/O2 /arch:AVX2 /DNDEBUG`) and the
  real-data end-to-end pipeline (`run_pipeline_real.bat`).

## Data availability

This repository contains code only. The genome x gene-family presence/absence matrix is
derived from a public, precomputed microbial pangenome resource (PanGBank, built with
PPanGGOLiN). Because of data-provider terms, the derived matrices and the source pangenome
file are NOT redistributed here. To reproduce:

1. Query the PanGBank REST API for the species-level pangenome of interest
   (`GET /pangenomes/?pangenome_name=<species>`), then download the pangenome file
   (`GET /pangenomes/{id}/file`). The study used the GTDB_refseq collection release 2.0.0
   pangenome of Klebsiella pneumoniae (2001 genomes, 54590 gene families).
2. Convert with `python src/gbm/h5_to_gbm1.py <pangenome.h5> <matrix.txt>`.
3. Run the pipeline (`run_pipeline_real.bat` shows the exact steps and parameters).

VFDB virulence factor database: http://www.mgc.ac.cn/VFs/ (VFDB_setA_pro).
ResFinder resistance database: https://bitbucket.org/genomicepidemiology/resfinder_db

## Baselines

pyseer (>=1.3.11, FE and LMM modes) and Scoary (1.1.4) are used as gene-level association
baselines. They are standard open-source tools installed from their official repositories;
please validate baseline runs on Linux/macOS (the Windows port used in early development
required two compatibility patches to Scoary and stubs for pyseer's VCF/enet paths; see
`../plans/03_实验fail记录.md`).

## License and attribution

MIT (see LICENSE). If you use this code, cite the associated study and credit
YalongZhang as the author.
