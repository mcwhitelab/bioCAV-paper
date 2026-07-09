source "$(dirname "$0")/../config/paths.sh"

python $BIOCAV_REPO/specific_scripts/compare_ec_tools.py \
    --tool-predictions  static/ec_annotations/EC-Bench/results_test_100_go_3.csv \
    --out-dir           results/ec_eval/ecbenchtest \
    --cav-results       results/ec_eval/ecbenchtest/eval_ec_results.tsv \
    --exclude           static/ec_annotations/EC-Bench/ec_train_overlap_pairs.tsv \
    --llr-threshold     2.3 \
    --val-pkl           results/ec_eval/ecbenchtest/eval_embeddings/val_proteins.span.fasta.pkl \
    --ec-base-dirs      $PROFAB_EC_SETS/ec_dataset_part* \
    --scaler-pkl        "$BIOCAV_PAPER_ROOT/shared_protein/reference_population/scaler_v1.pkl" \
    --figure-data-dir   "$BIOCAV_PAPER_ROOT/figures/figure_data"

# ecbenchprice comparison intentionally not run: the price-149 gold standard
# has label-quality issues that don't match the structural matches (see
# ec_eval/01_reformat_gold_standards.sh comment).
