source "$(dirname "$0")/../config/paths.sh"

PY=/groups/clairemcwhite/envs/core_pkgs4/bin/python

# Exports flat CSVs for Figure 5's three case-study pairs (neutrophil/breast,
# epithelial/lung10x, T_cell/colorectum) into ../figures/figure_data/, so
# figures.R can render native ggplot panels instead of embedding matplotlib
# PNGs. All three are restricted to donors with both conditions present
# (see scripts/export_fig5_case_studies.py and scripts/
# paired_donor_gene_corr.py for why -- pooled, donor-blind correlation is
# vulnerable to Simpson's-paradox confounds).

$PY "$(dirname "$0")/scripts/export_fig5_case_studies.py" "$BIOCAV_PAPER_ROOT/figures/figure_data"
