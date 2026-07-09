source "$(dirname "$0")/../config/paths.sh"

python "$(dirname "$0")/scripts/cav_similarity_heatmap.py" \
    --library cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d/ \
        --pkl     cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d/embeddings/cells.pkl \
    --out     cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d/cav_similarity.png \
    --width 20 \
    --height 20
#python "$(dirname "$0")/scripts/cav_similarity_heatmap.py" \
#    --library cav_library/9d8e5dca-03a3-457d-b7fb-844c75735c83/ \
#        --pkl     cav_library/9d8e5dca-03a3-457d-b7fb-844c75735c83/embeddings/cells.pkl \
#    --out     cav_library/9d8e5dca-03a3-457d-b7fb-844c75735c83/cav_similarity.png
