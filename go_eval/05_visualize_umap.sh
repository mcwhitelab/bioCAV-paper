source "$(dirname "$0")/../config/paths.sh"

VIZ_SCRIPT="$(dirname "$0")/../motif_clustering/cav_viz_pattern.py"

# GO terms aren't split into MF/BP/CC subdirectories on disk, so build plain
# include-list text files (one concept name per line) to restrict the glob
# to one namespace via --include-list. No symlinks or directory restructuring.
LISTS="$(dirname "$0")/results/go_namespace_lists"
if [ ! -f "$LISTS/mf_ids.txt" ]; then
    $CONDA_ENV_DIR/bin/python "$(dirname "$0")/scripts/build_go_namespace_lists.py" --out-dir "$LISTS"
fi

GO_PATTERN="$PROFAB_GO_SETS/go_dataset_part*/GO_*/random_positive_train_max1000_cav/concept_v1.npy"

OUT="$(dirname "$0")/results/figures"
mkdir -p "$OUT"

run_umap() {
    NAME="$1"; INCLUDE="$2"

    echo "=== $NAME : UMAP 2D interactive ==="
    $CONDA_ENV_DIR/bin/python "$VIZ_SCRIPT" \
        --cav-pattern   "$GO_PATTERN" \
        --reducer       umap \
        --interactive \
        --dims          2 \
        --include-list  "$INCLUDE" \
        --out           "$OUT/${NAME}_umap.html"

    echo "=== $NAME : UMAP 3D interactive ==="
    $CONDA_ENV_DIR/bin/python "$VIZ_SCRIPT" \
        --cav-pattern   "$GO_PATTERN" \
        --reducer       umap \
        --interactive \
        --dims          3 \
        --include-list  "$INCLUDE" \
        --out           "$OUT/${NAME}_umap_3d.html"

    echo "=== $NAME : t-SNE static PNG ==="
    $CONDA_ENV_DIR/bin/python "$VIZ_SCRIPT" \
        --cav-pattern   "$GO_PATTERN" \
        --reducer       tsne \
        --include-list  "$INCLUDE" \
        --out           "$OUT/${NAME}_direction_map_tsne.png"
}

run_umap mf "$LISTS/mf_ids.txt"
run_umap cc "$LISTS/cc_ids.txt"
run_umap bp "$LISTS/bp_ids.txt"

# BP has ~4700 concepts -- labels overlap into an unreadable mass, so also
# save an unlabeled version.
echo "=== bp : t-SNE static PNG (no labels) ==="
$CONDA_ENV_DIR/bin/python "$VIZ_SCRIPT" \
    --cav-pattern   "$GO_PATTERN" \
    --reducer       tsne \
    --include-list  "$LISTS/bp_ids.txt" \
    --no-labels \
    --out           "$OUT/bp_direction_map_tsne_nolabels.png"

EC_PATTERN="$PROFAB_EC_SETS/ec_dataset_part*/ecNo_*/random_positive_train_max1000_cav/concept_v1.npy"
COMBINED_ARGS=(
    --ec-pattern "$EC_PATTERN"
    --go-pattern "$GO_PATTERN"
    --mf-list    "$LISTS/mf_ids.txt"
    --bp-list    "$LISTS/bp_ids.txt"
    --cc-list    "$LISTS/cc_ids.txt"
    --ec-enzyme-dat "$(dirname "$0")/../ec_eval/static/ec_annotations/enzyme.dat"
    --ec-enzclass   "$(dirname "$0")/../ec_eval/static/ec_annotations/enzclass.txt"
    --go-obo        "$(dirname "$0")/static/go_annotations/go.obo"
)

echo "=== combined EC + GO(MF/BP/CC) : UMAP 3D interactive ==="
$CONDA_ENV_DIR/bin/python "$(dirname "$0")/../motif_clustering/combined_cav_umap_3d.py" \
    "${COMBINED_ARGS[@]}" --dims 3 --out "$OUT/combined_umap_3d.html"

echo "=== combined EC + GO(MF/BP/CC) : UMAP 2D interactive ==="
$CONDA_ENV_DIR/bin/python "$(dirname "$0")/../motif_clustering/combined_cav_umap_3d.py" \
    "${COMBINED_ARGS[@]}" --dims 2 --out "$OUT/combined_umap_2d.html"
