# Exact Bitmap Association Mining Framework

GPU-accelerated exact association mining for population-scale genomic feature analysis. Companion code repository for reproducibility.

> **Status**: under active development. Interfaces may change before the first stable release.

## What it does

- Represents genome × feature matrices as bitmaps and reduces support counting, co-occurrence, exclusion and conditional association to exact bulk Boolean/counting operations.
- Uses a measurement-driven intra-GPU scheduler that routes each operation between SIMT and matrix engines with double-buffered overlap, producing route-independent, bit-identical outputs.
- Provides exact (no sampling, no sketching, no approximation) candidate enumeration with statistical testing, including phylogeny-aware permutation nulls.

## Repository layout

```
src/         core library (counting engine, router, statistics)
tests/       unit and exactness tests (bit-identity checks, MD5 outputs)
tutorials/   minimal end-to-end examples
notebooks/   reproducible analysis notebooks
envs/        conda/docker environment definitions
```

## Installation

(Installation instructions to be added with the first release. Development environment: see `envs/`.)

## Data availability

**This repository does not distribute datasets.** All analyses use publicly available resources that users must download themselves from the original sources, e.g.:

- GTDB (Genome Taxonomy Database): https://gtdb.ecogenomic.org
- NCBI Assembly / GenBank: https://www.ncbi.nlm.nih.gov/assembly
- CARD (Comprehensive Antibiotic Resistance Database): https://card.mcmaster.ca
- VFDB (Virulence Factor Database): http://www.mgc.ac.cn/VFs/

Download manifests (accessions, URLs, SHA256 checksums) and preprocessing scripts are provided so that every derived dataset can be reproduced exactly. Due to licensing terms of the source databases, derived data cannot be redistributed here.

## Reproducing the analyses

Step-by-step reproduction instructions, pinned software versions, and expected checksums of intermediate outputs will be published with the preprint.

## Citation

Citation information will be added upon publication.

## License

MIT (see `LICENSE`).
