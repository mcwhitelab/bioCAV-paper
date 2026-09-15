source "$(dirname "$0")/../config/paths.sh"

PY=/groups/clairemcwhite/envs/core_pkgs4/bin/python

# Exports flat CSVs for Figure 6's three case-study pairs (neutrophil/breast,
# epithelial/lung10x, fibroblast/colorectum) into ../figures/figure_data/, so
# figures.R can render native ggplot panels instead of embedding matplotlib
# PNGs. All three are restricted to donors with both conditions present --
# pooled, donor-blind statistics are vulnerable to Simpson's-paradox
# confounds (see scripts/paired_donor_gene_corr.py). Both axes of the
# DE-vs-CAV scatter come from that same paired-donor cell set: the CAV
# correlations from scripts/paired_donor_gene_corr.py, the DE from
# scripts/paired_donor_de_mixedlm.py.

$PY "$(dirname "$0")/scripts/export_fig5_case_studies.py" "$BIOCAV_PAPER_ROOT/figures/figure_data"
