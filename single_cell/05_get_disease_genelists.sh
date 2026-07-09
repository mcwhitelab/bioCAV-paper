export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
source "$(dirname "$0")/../config/paths.sh"
LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d

python $BIOCAV_REPO/specific_scripts/cav_gene_correlation.py \
    --coords   $LIB/results/hierarchy/cell_coordinates.tsv \
    --h5ad     $LIB/data/cells.h5ad \
    --lib-dir  $LIB/ \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --level    L2 \
    --out-dir  $LIB/results/gene_correlation_celltype_tissue/ #\
    #--pairs \
    #    endothelial_cell__skin_epidermis__normal \
    #    epithelial_cell__skin_epidermis__normal \
    #    fibroblast__skin_epidermis__normal \
    #    T_cell__skin_epidermis__normal

