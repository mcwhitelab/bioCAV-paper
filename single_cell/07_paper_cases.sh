export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
source "$(dirname "$0")/../config/paths.sh"

LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d
PY=/groups/clairemcwhite/envs/core_pkgs4/bin/python

mkdir -p $LIB/results/figures

# Paper case studies:
#   - neutrophil, breast, normal vs. cancer   -- clear case (HLA-DR+/CD74 rise)
#   - epithelial_cell, skin_epidermis, melanoma -- clear case
#   - T_cell, skin_epidermis, melanoma          -- less clear case (contrast)
#
# Uses the mixedlm DE results (results/de_mixedlm/) rather than the DESeq2
# ones (results/de_vs_cav/) because DESeq2 pseudobulk failed for
# neutrophil__breast (too few donors); mixedlm has all three pairs, so this
# keeps the DE method consistent across the panel.

$PY $BIOCAV_REPO/specific_scripts/de_pseudobulk.py \
    --lib-dir       $LIB/ \
    --out-dir       $LIB/results/de_mixedlm/ \
    --cav-dir       $LIB/results/gene_correlation_celltype_tissue/ \
    --compare-only \
    --only-pairs \
        neutrophil__breast__normal_vs_breast_cancer \
        epithelial_cell__skin_epidermis__normal_vs_melanoma \
        T_cell__skin_epidermis__normal_vs_melanoma \
    --scatter-out   $LIB/results/figures/de_vs_cav_scatter_paper_cases.png \
    --lfc-threshold 1.5 \
    --label-top-n 15 \
    --show-continuum-genes \
    --continuum-top-n-each 4 \
    --continuum-top-n-each-override skin_epidermis=5 \
    --highlight-genes HLA-DRA HLA-DRB1 CD74

$PY "$(dirname "$0")/scripts/cav_continuum_viz.py" \
    --coords        $LIB/results/hierarchy/cell_coordinates.tsv \
    --h5ad          $LIB/data/cells.h5ad \
    --gene-corr-dir $LIB/results/gene_correlation_celltype_tissue/ \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --level         L2 \
    --pairs \
        neutrophil__breast__normal_vs_breast_cancer \
        epithelial_cell__skin_epidermis__normal_vs_melanoma \
        T_cell__skin_epidermis__normal_vs_melanoma \
    --title "Selected case studies — transcriptional continuum (L2)" \
    --out   $LIB/results/figures/cav_continuum_paper_cases.png \
    --top-n-each 5 \
    --style strips
