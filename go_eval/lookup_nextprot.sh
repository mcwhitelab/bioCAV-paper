source "$(dirname "$0")/../config/paths.sh"

f=$NEXTPROT_MF_FASTA

python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav $PROFAB_GO_SETS/go_dataset_part*/*/random_positive_train_max1000_cav \
    --cav-name-depth 2 --cav-name-strip "_random_positive_train_max1000_cav" \
    --embed-fasta $f \
    --scaler-pkl reference_population/scaler_v1.pkl \
    --out results/nextprot_mf_scores.tsv --k -1 \
    --model $ESM_MODEL_DIR \
    --embed-script $HF_EMBED_SCRIPT

