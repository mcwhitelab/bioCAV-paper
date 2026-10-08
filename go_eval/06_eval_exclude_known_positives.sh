#!/bin/bash
# Fig 2 GO temporal evaluation with known-positive exclusion (2026-09-29).
#
# A post-2021 annotation date does not make a (protein, term) pair new: about a
# third of the original validation pairs were proteins that proFAB already lists
# as positives for that term (Swiss-Prot <= 2020-01, true-path propagated), most
# of them in the very span the CAV was trained on. --exclude-known-positives drops
# pairs whose protein is in the term's proFAB random train or test positive span,
# then re-applies --min-n.
#
# The term set is held fixed to the shipped Fig 2 terms (static/fig2_go_terms_*.txt,
# MF 226 / BP 395 / CC 258): 175+ further GO CAVs were trained in June 2026, after
# the original eval, and would otherwise silently enlarge the benchmark with
# proteins that have neither embeddings nor DeepGo-SE predictions.
#
# Uses the CAVs currently on disk (retrained in the shared scaler frame 2026-09-16).
# Writes to fresh *_noknown out-dirs so no cached July scores are reused, and
# stages figure data in figures/figure_data_staging/noknown/ for preview; nothing
# in figures/figure_data/ is touched. Embeddings are reused from the original
# eval_embeddings/ (superset of the proteins needed).
set -euo pipefail
source "$(dirname "$0")/../config/paths.sh"
cd "$(dirname "$0")"

STAGE="$BIOCAV_PAPER_ROOT/figures/figure_data_staging/noknown"
SCALER="$BIOCAV_PAPER_ROOT/shared_protein/reference_population/scaler_v1.pkl"
mkdir -p "$STAGE"

for ont in mf bp cc; do
    OUT=results/temporal_eval_${ont}_noknown
    python $BIOCAV_REPO/specific_scripts/eval_go_temporal.py \
        --gold-standard  static/go_annotations/${ont}.gpad \
        --go-base-dirs   "$PROFAB_GO_SETS/go_dataset_part*" \
        --out-dir        $OUT/ \
        --eval-embed-dir results/temporal_eval_${ont}/eval_embeddings/ \
        --scaler-pkl     "$SCALER" \
        --min-n          5 \
        --date-cutoff    20211231 \
        --go-terms-file  static/fig2_go_terms_${ont}.txt \
        --exclude-known-positives

    python $BIOCAV_REPO/specific_scripts/summarize_temporal_eval.py \
        --results  $OUT/eval_temporal_results.tsv \
        --go-obo   static/go_annotations/go.obo \
        --out-dir  $OUT/ \
        --figure-data-dir "$STAGE" \
        --label $ont

    python $BIOCAV_REPO/specific_scripts/compare_tool_temporal.py \
        --tool-predictions  outside_tools/deepgose_preds_${ont}.tsv \
        --results           $OUT/eval_temporal_results.tsv \
        --per-term-summary  $OUT/eval_temporal_per_term_summary.tsv \
        --out-dir           $OUT/ \
        --output            $OUT/tool_comparison_deepgose.tsv \
        --tool-format long_tsv \
        --val-pkl results/temporal_eval_${ont}/eval_embeddings/val_proteins.span.fasta.pkl \
        --go-base-dirs $PROFAB_GO_SETS/go_dataset_part* \
        --scaler-pkl "$SCALER" \
        --figure-data-dir "$STAGE" \
        --label $ont \
        --llr-threshold 2.3
done

FIG_DATA_DIR="$STAGE" python rank_all_terms/compute_ranks_all_terms.py
