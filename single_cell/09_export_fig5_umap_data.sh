source "$(dirname "$0")/../config/paths.sh"

LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d
PY=/groups/clairemcwhite/envs/core_pkgs4/bin/python

# Exports the overview-row t-SNE data for the top of Figure 5: raw
# (unmodified individual baseline CAVs, no averaging/orthogonalization) +
# the CAV hierarchy's three levels using the swapped context-group-condition
# level order (L0=tissue axis/base, L1=cell_type axis with tissue projected
# away, L2=condition axis with tissue+cell_type projected away). L2 spans
# the same subspace either way (tissue+cell_type together), so this is
# consistent with the L2 axes used in results/hierarchy/ elsewhere in
# Figure 5 -- only L0/L1 change which variable is projected away first.
#
# --clean-pairs-csv restricts the cell pool to (cell_type, tissue) combos
# where normal and disease cells were sequenced with the same dominant
# assay (>=80% on both sides) -- so any condition separation visible in
# the t-SNE reflects biology, not a sequencing-technology batch effect
# (see static/clean_assay_matched_pairs.csv for how this was derived).
#
# --paired-donors-only additionally drops cells from donors that only have
# one condition present within their (cell_type, tissue) group, matching
# the paired-donor restriction already used for panels E-G's gene
# correlations (see export_fig5_case_studies.py / paired_donor_gene_corr.py)
# -- keeps the overview t-SNE population consistent with the rest of
# Figure 5 instead of drawing from a broader any-donor pool.

$PY "$(dirname "$0")/scripts/export_fig5_umap_data.py" \
    --coords        $LIB/results/hierarchy_tissue_celltype/cell_coordinates.tsv \
    --h5ad          $LIB/data/cells.h5ad \
    --lib-dir       $LIB \
    --raw-pkl       $LIB/embeddings/cells.pkl \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --baseline-value normal \
    --max-cells     20000 \
    --clean-pairs-csv "$(dirname "$0")/static/clean_assay_matched_pairs.csv" \
    --paired-donors-only \
    --figure-data-dir "$BIOCAV_PAPER_ROOT/figures/figure_data"
