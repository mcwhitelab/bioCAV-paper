source "$(dirname "$0")/../config/paths.sh"

# Single-term debug example:
#python $BIOCAV_REPO/specific_scripts/inspect_go_predictions.py \
#    --go-term          GO:0070403 \
#    --results          results/temporal_eval_mf/eval_temporal_results.tsv \
#    --out-dir          results/temporal_eval_mf/ \
#    --tool-predictions outside_tools/deepgose_preds_mf.tsv \
#    --tool-format      long_tsv

for go in mf bp cc;
do
python $BIOCAV_REPO/specific_scripts/inspect_go_predictions.py \
    --all-terms \
    --results results/temporal_eval_${go}/eval_temporal_results.tsv \
    --out-dir results/temporal_eval_${go}/ \
    --tool-predictions outside_tools/deepgose_preds_${go}.tsv \
    --output results/temporal_eval_${go}/all_term_predictions.tsv
done
