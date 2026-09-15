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

# ---------------------------------------------------------------------------
# DeepGoWeb — DISABLED (dropped from Figure 2, 2026-07-12).
#
# DeepGoWeb runs DeepGOPlus = CNN + DIAMOND sequence-similarity. The DIAMOND
# component transfers GO annotations from homologous training-set proteins, so
# it is NOT a sequence-representation-only method and is not comparable to CAV /
# DeepGo-SE for the "what embeddings capture" story. (See DeepGOWeb, NAR 2021;
# DeepGOPlus, Bioinformatics 2020.) Additional problems found before dropping:
# our batch scores (API version 1.0.27, threshold 0) sit on a different/lower
# scale than the website and could not be reconciled, and the batch covers only
# validation proteins (no negatives), so AUC/AUPR were zero-fill inflated.
#
# Recipe kept for reference. To re-enable, uncomment. It writes to a SEPARATE
# figure_data/deepgoweb/ subdir (so it never clobbers the DeepGoSE data above)
# and applies DeepGoWeb's suggested 0.1 minimum confidence via
# --tool-score-threshold. figures.R would also need its 3-way version restored.
#
# mkdir -p "$BIOCAV_PAPER_ROOT/figures/figure_data/deepgoweb"
# for go in mf bp cc;
# do
# python $BIOCAV_REPO/specific_scripts/compare_tool_temporal.py \
#     --tool-predictions  outside_tools/val_proteins_${go}.deepgo.csv \
#     --results           results/temporal_eval_${go}/eval_temporal_results.tsv \
#     --per-term-summary  results/temporal_eval_${go}/eval_temporal_per_term_summary.tsv \
#     --out-dir           results/temporal_eval_${go}/ \
#     --output            results/temporal_eval_${go}/tool_comparison_deepgoweb.tsv \
#     --tool-format wide_csv \
#     --tool-score-threshold 0.1 \
#     --val-pkl results/temporal_eval_${go}/eval_embeddings/val_proteins.span.fasta.pkl \
#     --go-base-dirs $PROFAB_GO_SETS/go_dataset_part* \
#     --scaler-pkl "$BIOCAV_PAPER_ROOT/shared_protein/reference_population/scaler_v1.pkl" \
#     --figure-data-dir "$BIOCAV_PAPER_ROOT/figures/figure_data/deepgoweb" \
#     --label $go \
#     --llr-threshold 2.3
# done
