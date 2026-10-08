source "$(dirname "$0")/../config/paths.sh"

VIZ_SCRIPT="$(dirname "$0")/../motif_clustering/cav_viz_pattern.py"
EC_PATTERN="$PROFAB_EC_SETS/ec_dataset_part*/ecNo_*/random_positive_train_max1000_cav/concept_v1.npy"

OUT="$(dirname "$0")/results/figures"
mkdir -p "$OUT"

echo "=== ec : UMAP 2D interactive ==="
$CONDA_ENV_DIR/bin/python "$VIZ_SCRIPT" \
    --cav-pattern "$EC_PATTERN" \
    --reducer     umap \
    --interactive \
    --dims        2 \
    --out         "$OUT/ec_umap.html"

echo "=== ec : UMAP 3D interactive ==="
$CONDA_ENV_DIR/bin/python "$VIZ_SCRIPT" \
    --cav-pattern "$EC_PATTERN" \
    --reducer     umap \
    --interactive \
    --dims        3 \
    --out         "$OUT/ec_umap_3d.html"

echo "=== ec : t-SNE static PNG ==="
$CONDA_ENV_DIR/bin/python "$VIZ_SCRIPT" \
    --cav-pattern "$EC_PATTERN" \
    --reducer     tsne \
    --out         "$OUT/ec_direction_map_tsne.png"
