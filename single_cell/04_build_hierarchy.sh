export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
source "$(dirname "$0")/../config/paths.sh"
LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d
PCA=reference_population/global_pca_v1.pkl

mkdir $LIB/results/figures
mkdir $LIB/results/hierarchy


# Step 2: direction map (runs now — no embeddings needed)
#python "$(dirname "$0")/scripts/cav_viz.py" \
#    --plot    direction-map \
#    --lib-dir $LIB \
#    --pca-pkl $PCA \
#    --out     $LIB/results/figures/direction_map.png

## Step 3: build hierarchy coordinates (needs embeddings)
#python $BIOCAV_REPO/specific_scripts/cav_hierarchy.py \
#    --lib-dir $LIB \
#    --pkl     $LIB/embeddings/cells.pkl \
#    --pca-pkl $PCA \
#    --out     $LIB/results/hierarchy/
#
## Step 4: condition scatter
#python "$(dirname "$0")/scripts/cav_viz.py" \
#    --plot    condition-scatter \
#    --lib-dir $LIB \
#    --pca-pkl $PCA \
#    --pkl     $LIB/embeddings/cells.pkl \
#    --obs     $LIB/data/cells.h5ad \
#    --out     $LIB/results/figures/scatter/
#
## Step 5: CAV-space UMAP (needs hierarchy output from step 3)
#python "$(dirname "$0")/scripts/cav_viz.py" \
#    --plot      cav-umap \
#    --coords    $LIB/results/hierarchy/cell_coordinates.tsv \
#    --obs       $LIB/data/cells.h5ad \
#    --shape-col tissue \
#    --condition-col disease \
#    --level     L2 \
#    --out       $LIB/results/figures/cav_tsne_L2.png \
#    --reducer tsne

LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d
PCA=reference_population/global_pca_v1.pkl




# Default: L0=cell_type, L1=tissue (current behaviour)
python $BIOCAV_REPO/specific_scripts/cav_hierarchy.py \
    --lib-dir $LIB/ \
    --pkl     $LIB/embeddings/cells.pkl \
    --out     $LIB/results/hierarchy_celltype_tissue/ \
    --level-order group-context-condition

# Swapped: L0=tissue, L1=cell_type
python $BIOCAV_REPO/specific_scripts/cav_hierarchy.py \
    --lib-dir $LIB/ \
    --pkl     $LIB/embeddings/cells.pkl \
    --out     $LIB/results/hierarchy_tissue_celltype/ \
    --level-order context-group-condition



for level in L0 L1 L2; do

    python "$(dirname "$0")/scripts/cav_viz.py" \
        --plot cav-umap \
        --coords $LIB/results/hierarchy_tissue_celltype/cell_coordinates.tsv \
        --obs $LIB/data/cells.h5ad \
        --color-col tissue \
        --shape-col cell_type \
	--condition-col disease \
	--baseline-value normal \
        --level $level \
        --out $LIB/results/figures/cav_tsne_${level}_tcd_tissue.png \
        --reducer tsne &


    python "$(dirname "$0")/scripts/cav_viz.py" \
        --plot cav-umap \
        --coords $LIB/results/hierarchy/cell_coordinates.tsv \
        --obs $LIB/data/cells.h5ad \
        --color-col cell_type \
        --shape-col tissue \
	--condition-col disease \
	--baseline-value normal \
        --level $level \
        --out $LIB/results/figures/cav_tsne_${level}_tcd_celltype.png \
        --reducer tsne &

done
wait




