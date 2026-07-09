source "$(dirname "$0")/../config/paths.sh"

python $BIOCAV_REPO/specific_scripts/summarize_ec_eval.py \
    --results  results/ec_eval/ecbenchtest/eval_ec_results.tsv \
    --out-dir  results/ec_eval/ecbenchtest/ \
    --exclude  static/ec_annotations/EC-Bench/ec_train_overlap_pairs.tsv

python $BIOCAV_REPO/specific_scripts/summarize_ec_eval.py \
    --results  results/ec_eval/ecbenchprice/eval_ec_results.tsv \
    --out-dir  results/ec_eval/ecbenchprice/
