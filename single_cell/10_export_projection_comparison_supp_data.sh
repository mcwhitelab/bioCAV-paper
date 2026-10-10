source "$(dirname "$0")/../config/paths.sh"

PY=/groups/clairemcwhite/envs/core_pkgs4/bin/python

# Candidate supplemental figure: does orthogonalizing the condition CAV
# against tissue/cell-type structure change the DE-vs-CAV relationship,
# compared to using the raw condition CAV directly? See
# scripts/export_projection_comparison_supp_data.py for the three stages
# computed (raw / L0_removed / L2) and why there are three, not four.

# Since 2026-10-08: official-tokenizer embeddings (cls_penult) and the Fig 6 cases, with the
# cancer - normal row added (METHODS §7.1). Run from single_cell/ (relative paths).
(cd "$(dirname "$0")" && $PY scripts/export_projection_comparison_supp_data.py \
    --figure-data-dir "$BIOCAV_PAPER_ROOT/figures/figure_data" \
    --run-dir /xdisk/clairemcwhite/clairemcwhite/sc_geneformer_official/runs/cls_penult)
