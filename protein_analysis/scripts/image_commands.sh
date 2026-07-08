source "$(dirname "$0")/../../config/paths.sh"


python "$(dirname "$0")/cav_viz_pattern.py" \
    --cav-pattern   "$AHMAD_CAV_OUTPUTS_DIR/cavs/*/L25_concept_v1.npy" \
    --reducer       umap \
    --clan-annotations results/figures/Pfam-A.clans.tsv \
    --gif-out       results/figures/pfam_umap_sim_3d.gif \
    --gif-frames    72 \
    --gif-fps       10 \
    --interactive \
    --distance-mode cophenetic \
    --n-neighbors 10 \
    --min-dist 0.3 \
    --dims          3 \
    --pfam-annotations results/figures/pfamA.txt \
    --out           results/figures/pfam_umap_sim_3d.html

python "$(dirname "$0")/cav_viz_pattern.py"     --cav-pattern "$AHMAD_CAV_OUTPUTS_DIR/cavs/*/L25_concept_v1.npy"     --reducer umap     --pfam-annotations results/figures/pfamA.txt     --interactive     --out results/figures/pfam_umap_sim.html --dims 2 --clan-annotations results/figures/Pfam-A.clans.tsv      --distance-mode cophenetic 




python "$(dirname "$0")/cav_viz_pattern.py" \
    --cav-pattern   "$AHMAD_CAV_OUTPUTS_DIR/cavs/*/L25_concept_v1.npy" \
    --reducer       umap \
    --clan-annotations results/figures/Pfam-A.clans.tsv \
    --gif-out       results/figures/pfam_umap_3d.gif \
    --gif-frames    72 \
    --gif-fps       72 \
    --interactive \
    --dims          3 \
    --pfam-annotations results/figures/pfamA.txt \
    --out           results/figures/pfam_umap_3d.html


python "$(dirname "$0")/cav_viz_pattern.py"     --cav-pattern "$AHMAD_CAV_OUTPUTS_DIR/cavs/*/L25_concept_v1.npy"     --reducer umap     --pfam-annotations results/figures/pfamA.txt     --interactive     --out results/figures/pfam_umap.html --dims 2 --clan-annotations results/figures/Pfam-A.clans.tsv 
