export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
source "$(dirname "$0")/../config/paths.sh"
LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d

# Three variants, differing only in which hierarchy run's cell_coordinates.tsv
# (and hence L2 score columns) they correlate against, and group/context
# column order. All three now library-size-normalize + log1p expression
# before correlating (cav_gene_correlation.py's load_expression default) --
# without this, a gene's correlation with the CAV axis can be dominated by
# per-cell total library size (sequencing depth / RNA content) rather than
# that gene's own biology. See MALAT1/fibroblast/melanoma for a worked
# example: raw-count r=+0.46, normalized r=-0.62 (agreeing with DE direction).

python $BIOCAV_REPO/specific_scripts/cav_gene_correlation.py \
    --coords   $LIB/results/hierarchy/cell_coordinates.tsv \
    --h5ad     $LIB/data/cells.h5ad \
    --lib-dir  $LIB/ \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --level    L2 \
    --out-dir  $LIB/results/gene_correlation/

python $BIOCAV_REPO/specific_scripts/cav_gene_correlation.py \
    --coords   $LIB/results/hierarchy_celltype_tissue/cell_coordinates.tsv \
    --h5ad     $LIB/data/cells.h5ad \
    --lib-dir  $LIB/ \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --level    L2 \
    --out-dir  $LIB/results/gene_correlation_celltype_tissue/

python $BIOCAV_REPO/specific_scripts/cav_gene_correlation.py \
    --coords   $LIB/results/hierarchy_tissue_celltype/cell_coordinates.tsv \
    --h5ad     $LIB/data/cells.h5ad \
    --lib-dir  $LIB/ \
    --group-col     tissue \
    --context-col   cell_type \
    --condition-col disease \
    --level    L2 \
    --out-dir  $LIB/results/gene_correlation_tissue_celltype/
