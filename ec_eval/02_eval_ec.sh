source "$(dirname "$0")/../config/paths.sh"

python $BIOCAV_REPO/specific_scripts/eval_ec.py \
    --gold-standard  static/ec_annotations/EC-Bench/test_ec_long.tsv \
    --ec-base-dirs   $PROFAB_EC_SETS/ec_dataset_part* \
    --out-dir        results/ec_eval/ecbenchtest/ \
    --scaler-pkl     "$BIOCAV_PAPER_ROOT/shared_protein/reference_population/scaler_v1.pkl" \
    --min-n          1

python $BIOCAV_REPO/specific_scripts/eval_ec.py \
    --gold-standard  static/ec_annotations/EC-Bench/price-149_long.tsv \
    --ec-base-dirs   $PROFAB_EC_SETS/ec_dataset_part* \
    --out-dir        results/ec_eval/ecbenchprice/ \
    --scaler-pkl     "$BIOCAV_PAPER_ROOT/shared_protein/reference_population/scaler_v1.pkl" \
    --val-fasta      static/ec_annotations/EC-Bench/price-149.fasta \
    --min-n          1
