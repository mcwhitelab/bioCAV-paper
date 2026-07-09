source "$(dirname "$0")/../config/paths.sh"

for go in mf bp cc;
do
python $BIOCAV_REPO/specific_scripts/compare_tool_temporal.py \
    --tool-predictions  outside_tools/deepgose_preds_${go}.tsv \
    --results           results/temporal_eval_${go}/eval_temporal_results.tsv \
    --per-term-summary  results/temporal_eval_${go}/eval_temporal_per_term_summary.tsv \
    --out-dir           results/temporal_eval_${go}/ \
    --output            results/temporal_eval_${go}/tool_comparison_deepgose.tsv \
    --tool-format long_tsv \
    --val-pkl results/temporal_eval_${go}/eval_embeddings/val_proteins.span.fasta.pkl \
    --go-base-dirs $PROFAB_GO_SETS/go_dataset_part* \
    --scaler-pkl "$BIOCAV_PAPER_ROOT/shared_protein/reference_population/scaler_v1.pkl" \
    --figure-data-dir "$BIOCAV_PAPER_ROOT/figures/figure_data" \
    --label $go \
    --llr-threshold 2.3

done

# For deepgoweb (unused)
#for go in bp cc mf;
#do
#python $BIOCAV_REPO/specific_scripts/compare_tool_temporal.py \
#    --tool-predictions  outside_tools/val_proteins_${go}.deepgo.csv \
#    --results           results/temporal_eval_${go}/eval_temporal_results.tsv \
#    --per-term-summary  results/temporal_eval_${go}/eval_temporal_per_term_summary.tsv \
#    --out-dir           results/temporal_eval_${go}/ \
#    --output            results/temporal_eval_${go}/tool_comparison.tsv
#done
