source "$(dirname "$0")/../config/paths.sh"

python "$(dirname "$0")/cav_cluster.py" \
    --cav-pattern      "$AHMAD_CAV_OUTPUTS_DIR/cavs/*/L25_concept_v1.npy" \
    --pfam-annotations results/figures/pfamA.txt \
    --clan-annotations results/figures/Pfam-A.clans.tsv \
    --out              results/figures/pfam_clusters.csv

