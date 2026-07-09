source "$(dirname "$0")/../config/paths.sh"

# Not using Price, too many gold standard labels seem wrong -- kept for
# reference / reproducibility only.
python $BIOCAV_REPO/specific_scripts/reformat_ec_goldstandard.py \
    --input  static/ec_annotations/EC-Bench/test_ec.csv \
    --output static/ec_annotations/EC-Bench/test_ec_long.tsv

python $BIOCAV_REPO/specific_scripts/reformat_ec_goldstandard.py \
    --input  static/ec_annotations/EC-Bench/price-149.csv \
    --output static/ec_annotations/EC-Bench/price-149_long.tsv
