#!/usr/bin/env bash
# Re-classify all 6,677 GO CAVs after the retrain.  Expect own=0.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
/groups/clairemcwhite/envs/core_pkgs4/bin/python "$HERE/classify_scalers.py"
echo "--- frame counts after retrain ---"
awk -F'\t' 'NR>1{c[$3]++} END{for(k in c) print k, c[k]}' "$HERE/go_scaler_frames.tsv"
