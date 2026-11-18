#!/usr/bin/env bash
# ============================================================================
# PISCES one-command pipeline: any feature x genome presence/absence table
# -> exact all-pairs co-occurrence counts -> Fisher exact + BH-FDR candidates
#
# Usage:  ./pisces_run.sh <input> <outdir> [minc]
#   input   presence/absence table: roary gene_presence_absence.csv,
#           Panaroo gene_presence_absence.Rtab, PPanGGOLiN/PanGBank HDF5,
#           plain CSV/TSV matrix, or two-column long list
#           (format auto-detected; see src/gbm/to_gbm1.py --help)
#   outdir  output directory (created if missing)
#   minc    minimum support per feature (default 5); pairs whose joint support
#           cannot reach minc are never materialized
#
# Engines must be built first:
#   g++ -O3 -mavx2 -DNDEBUG -std=c++17 -o src/gbm/gbm_count src/gbm/gbm_count.cpp
#   g++ -O3 -mavx2 -DNDEBUG -std=c++17 -o src/gbm/gbm_stats src/gbm/gbm_stats.cpp
# ============================================================================
set -euo pipefail
if [ $# -lt 2 ]; then echo "Usage: pisces_run.sh <input> <outdir> [minc]"; exit 1; fi
IN="$1"; OUT="$2"; MINC="${3:-5}"
HERE="$(cd "$(dirname "$0")" && pwd)"
GBM1="$OUT/input.gbm1.txt"
PAIRS="$OUT/pairs.tsv"
CAND="$OUT/candidates.tsv"
mkdir -p "$OUT"

command -v python3 >/dev/null || { echo "ERROR: python3 not found on PATH"; exit 1; }
[ -x "$HERE/src/gbm/gbm_count" ] || { echo "ERROR: build gbm_count first (see README, Quick start)"; exit 1; }
[ -x "$HERE/src/gbm/gbm_stats" ] || { echo "ERROR: build gbm_stats first (see README, Quick start)"; exit 1; }

echo "[1/3] encode input to GBM1 bitmap (deterministic, one pass)"
python3 "$HERE/src/gbm/to_gbm1.py" "$IN" "$GBM1"

# GBM1 header line is: GBM1 <genomes N> <families M>
read -r _ N M < "$GBM1"

echo "[2/3] exact all-pairs count (minc=$MINC)"
"$HERE/src/gbm/gbm_count" --count "$GBM1" "$PAIRS" "$MINC"

echo "[3/3] Fisher two-sided p + odds ratio + BH-FDR (screened and all-pairs)"
"$HERE/src/gbm/gbm_stats" --pairs "$PAIRS" "$CAND" "$M"

echo
echo "PISCES pipeline done. Candidate table: $CAND"
echo "Columns: famA famB c11 c10 c01 c00 odds_ratio p_fisher padj_BH_screened padj_BH_allpairs rank"
