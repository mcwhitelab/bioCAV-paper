source "$(dirname "$0")/../config/paths.sh"

PY=/groups/clairemcwhite/envs/core_pkgs4/bin/python

# Exports flat CSVs for Figure 6's three case-study pairs (neutrophil/breast,
# epithelial/lung10x, fibroblast/colorectum) into ../figures/figure_data/, so
# figures.R can render native ggplot panels instead of embedding matplotlib
# PNGs. All three are restricted to donors with both conditions present --
# pooled, donor-blind statistics are vulnerable to Simpson's-paradox
# confounds (see scripts/paired_donor_gene_corr.py). Both axes of the
# DE-vs-CAV scatter come from that same paired-donor cell set: the CAV
# correlations from scripts/paired_donor_gene_corr.py, the DE from pseudobulk
# DESeq2 (scripts/paired_donor_pseudobulk.py, then
# scripts/paired_donor_pseudobulk_de.R 10 in the pseudobulk_de R env), which
# additionally drops pseudobulk samples of < 10 cells.

# Since 2026-10-07 Fig 6 uses the converged C=1 PCA-128 refit (convergence_check arm conv1:
# 01_refit_cavs.py -> 02_build_coords.py -> 03_gene_corr.sh), not the library's stored CAVs.
CC="$(dirname "$0")/convergence_check"
$PY "$(dirname "$0")/scripts/export_fig5_case_studies.py" "$BIOCAV_PAPER_ROOT/figures/figure_data" \
  --coords "$CC/coords_conv1.tsv" --gene-corr-dir "$CC/gene_corr_conv1"

# Lower half of Fig 6D (fibroblast/colorectum fibroblast-program panel; replaced
# that column's continuum strips on 2026-10-05). Run from single_cell/ (relative
# cav_library path).
(cd "$(dirname "$0")" && $PY scripts/export_fig6_fib_modules.py "$BIOCAV_PAPER_ROOT/figures/figure_data" --coords=convergence_check/coords_conv1.tsv)
