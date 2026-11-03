#!/usr/bin/env bash
# End-to-end demo on the bundled synthetic dataset (Linux/macOS).
# Usage: bash demo/run_demo.sh [path-to-exe-dir]
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
BIN="${1:-$root/dist/cpu-windows-x64}"
[ -x "$BIN/gbm_count" ] || BIN="$root/src/gbm"   # fallback: in-tree build dir
mkdir -p "$here/out"

echo "[1/3] exact all-pairs count + hypergeometric screen (minc=3)"
"$BIN/gbm_count" --count "$root/data/demo/demo_matrix.txt" "$here/out/demo_pairs.tsv" 3 screen

echo "[2/3] Fisher exact test per screened pair"
"$BIN/gbm_stats" --pairs "$here/out/demo_pairs.tsv" "$here/out/demo_stats.tsv"

echo "[3/3] check against planted ground truth"
python3 "$here/check_demo.py" "$here/out/demo_stats.tsv"
