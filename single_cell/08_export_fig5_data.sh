source "$(dirname "$0")/../config/paths.sh"

PY=/groups/clairemcwhite/envs/core_pkgs4/bin/python

# Exports flat CSVs for Figure 6's three case-study pairs (neutrophil/breast,
# epithelial/lung10x, neutrophil/lung) into ../figures/figure_data/, so
# figures.R can render native ggplot panels instead of embedding matplotlib
# PNGs. All three are restricted to donors with both conditions present --
# pooled, donor-blind statistics are vulnerable to Simpson's-paradox
# confounds (see scripts/paired_donor_gene_corr.py). Both axes of the
# DE-vs-CAV scatter come from that same paired-donor cell set: the CAV
# correlations from scripts/paired_donor_gene_corr.py, the DE from pseudobulk
# DESeq2 (convergence_check/13a_pseudobulk_allpairs.py, then
# 13b_pseudobulk_de_allpairs.R 10 in the pseudobulk_de R env), which
# additionally drops pseudobulk samples of < 10 cells.

# Since 2026-10-08 Fig 6 uses the official Geneformer V2 tokenization (cls_penult readout)
# and the cancer - normal CAV difference axis; see official_tokenizer/README.md and
# METHODS §6. Inputs (on xdisk):
#   coords_diff.tsv          official_tokenizer/16a_build_diff_coords.py  (cell scores)
#   gene_corr_allpairs_diff/ scripts/paired_donor_gene_corr.py --coords coords_diff.tsv
# The fibroblast/colorectum module panel (export_fig6_fib_modules.py) is no longer
# exported; previous data in figure_data/_backup_pre_fig6_official_diff_2026-10-08/.
RUN=/xdisk/clairemcwhite/clairemcwhite/sc_geneformer_official/runs/cls_penult
(cd "$(dirname "$0")" && $PY scripts/export_fig5_case_studies.py "$BIOCAV_PAPER_ROOT/figures/figure_data" \
  --coords "$RUN/coords_diff.tsv" --gene-corr-dir "$RUN/gene_corr_allpairs_diff")
