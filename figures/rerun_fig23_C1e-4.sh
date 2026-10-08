#!/bin/bash
# -----------------------------------------------------------------------------
# Rebuild every Fig 2 / Fig 3 data file with the C=1e-4 converged GO/EC CAVs
# (2026-10-01). Each step is the panel's normal pipeline, unchanged, pointed at a
# symlink overlay of proFAB in which only concept_v1.npy differs
# (distribution/retrain_go_ec_C1e-4/04_build_overlay.py). Writes ONLY to $STAGE and to
# fresh *_C1e-4 results dirs; promotion into figure_data/ is a separate step.
#
#   LIB=new (default)  overlay with the C=1e-4 CAVs, all 6,678 GO + 953 EC (incl. the two
#                      new-only terms, 2026-10-01) -> figure_data_staging/C1e-4
#   LIB=old            shipped proFAB CAVs (regression) -> $STAGE given by caller
#
# Fig 2 B, D, E   go_eval temporal eval (06_eval_exclude_known_positives.sh steps)
# Fig 2 C         go_eval/rank_C1e-4/rerank_go.py (BLAST extended by go_eval/rank_C1e-4/extend_blast.sh)
# Fig 3 A         ec_eval/rank_C1e-4/rerank_ec.py
# Fig 3 B         hierarchy_decay/ec_hierarchy_decay.py (5 proteins; C0HLV0 -> P0DUE1 2026-10-06)
# Fig 3 D, E      go_ec_crossval/go_ec_cosine_validation.py (+ 20k background sample, rng 0)
# Fig 3 C         motif_clustering/combined_cav_umap_3d.py --dims 2 --coords-csv-out
# -----------------------------------------------------------------------------
set -o pipefail
PAPER=/groups/clairemcwhite/cav_workspace/bioCAV-paper
source "$PAPER/config/paths.sh"
set -eu
PY=$CONDA_ENV_DIR/bin/python
LIB=${LIB:-new}
if [ "$LIB" = new ]; then
    BASE=/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay
    STAGE=${STAGE:-$PAPER/figures/figure_data_staging/C1e-4}
    TAG=C1e-4
else
    BASE=/xdisk/clairemcwhite/proFAB
    STAGE=${STAGE:?set STAGE for LIB=old}
    TAG=regress_old
fi
GO_SETS=$BASE/GO_sets; EC_SETS=$BASE/EC_sets
SCALER=$PAPER/shared_protein/reference_population/scaler_v1.pkl
mkdir -p "$STAGE"
STEPS=${STEPS:-"go_eval rank_go rank_ec decay cosine umap"}
echo "LIB=$LIB BASE=$BASE STAGE=$STAGE STEPS=$STEPS"

for step in $STEPS; do case $step in
go_eval)  # Fig 2 B, D, E (and the eval-set pairs file)
    cd "$PAPER/go_eval"
    for ont in mf bp cc; do
        OUT=results/temporal_eval_${ont}_noknown_$TAG
        $PY $BIOCAV_REPO/specific_scripts/eval_go_temporal.py \
            --gold-standard static/go_annotations/${ont}.gpad \
            --go-base-dirs $GO_SETS/go_dataset_part* \
            --out-dir $OUT/ \
            --eval-embed-dir results/temporal_eval_${ont}/eval_embeddings/ \
            --scaler-pkl "$SCALER" --min-n 5 --date-cutoff 20211231 \
            --go-terms-file static/fig2_go_terms_${ont}.txt \
            --exclude-known-positives
        $PY $BIOCAV_REPO/specific_scripts/summarize_temporal_eval.py \
            --results $OUT/eval_temporal_results.tsv --go-obo static/go_annotations/go.obo \
            --out-dir $OUT/ --figure-data-dir "$STAGE" --label $ont
        $PY $BIOCAV_REPO/specific_scripts/compare_tool_temporal.py \
            --tool-predictions outside_tools/deepgose_preds_${ont}.tsv \
            --results $OUT/eval_temporal_results.tsv \
            --per-term-summary $OUT/eval_temporal_per_term_summary.tsv \
            --out-dir $OUT/ --output $OUT/tool_comparison_deepgose.tsv --tool-format long_tsv \
            --val-pkl results/temporal_eval_${ont}/eval_embeddings/val_proteins.span.fasta.pkl \
            --go-base-dirs $GO_SETS/go_dataset_part* --scaler-pkl "$SCALER" \
            --figure-data-dir "$STAGE" --label $ont --llr-threshold 2.3
    done ;;
rank_go)  # Fig 2 C; keeps the 06 pairs file as go_specificity_ranks_{ont}.eval.csv
    for ont in mf bp cc; do
        [ -f "$STAGE/go_specificity_ranks_${ont}.csv" ] && \
            cp "$STAGE/go_specificity_ranks_${ont}.csv" "$STAGE/go_specificity_ranks_${ont}.eval.csv"
    done
    $PY "$PAPER/go_eval/rank_C1e-4/rerank_go.py" --lib $LIB --out "$STAGE" ;;
rank_ec)  # Fig 3 A
    $PY "$PAPER/ec_eval/rank_C1e-4/rerank_ec.py" --lib $LIB --out "$STAGE" ;;
decay)    # Fig 3 B
    cd "$PAPER"
    $PY hierarchy_decay/ec_hierarchy_decay.py \
        --val-pkl ec_eval/static/ec_annotations/EC-Bench/test_ec.fasta.pkl \
        --gold ec_eval/static/ec_annotations/EC-Bench/test_ec_long.tsv \
        --ec-base-dirs $EC_SETS/ec_dataset_part* --scaler-pkl "$SCALER" \
        --proteins C0HLQ6:1.2.1.12 P0DQP1:2.1.1.72 P0DUB9:3.1.1.3 P0DUE1:3.2.2.6 A0A345BJN1:4.3.1.24 \
        --out-dir hierarchy_decay/results_$TAG
    cp hierarchy_decay/results_$TAG/ec_hierarchy_decay.csv "$STAGE/" ;;
cosine)   # Fig 3 D, E
    cd "$PAPER/go_ec_crossval"
    $PY go_ec_cosine_validation.py --ec2go static/ec2go.txt \
        --go-base-dirs $GO_SETS/go_dataset_part* --ec-base-dirs $EC_SETS/ec_dataset_part* \
        --out-dir results_$TAG
    $PY - results_$TAG "$STAGE" <<'PY'
import sys, numpy as np, pandas as pd
res, stage = sys.argv[1:]
pd.read_csv(f"{res}/go_ec_cosine_validation.csv").to_csv(f"{stage}/go_ec_cosine_pairs.csv", index=False)
bg = np.load(f"{res}/go_ec_background_sims.npy")
s = np.random.default_rng(0).choice(bg, size=min(20000, len(bg)), replace=False)  # = plot_go_ec_cosine_validation.py
pd.DataFrame({"sim": s}).to_csv(f"{stage}/go_ec_cosine_background_sample.csv", index=False)
PY
    ;;
umap)     # Fig 3 C
    cd "$PAPER/go_eval"
    if [ "$LIB" = new ]; then  # namespace include-lists for the overlay's GO CAVs (+ GO:0032436)
        LISTS=results/go_namespace_lists_$TAG
        $PY scripts/build_go_namespace_lists.py --profab-go "$GO_SETS" \
            --go-obo static/go_annotations/go.obo --out-dir $LISTS
    else
        LISTS=results/go_namespace_lists
    fi
    mkdir -p results/figures_$TAG
    $PY ../motif_clustering/combined_cav_umap_3d.py \
        --ec-pattern "$EC_SETS/ec_dataset_part*/ecNo_*/random_positive_train_max1000_cav/concept_v1.npy" \
        --go-pattern "$GO_SETS/go_dataset_part*/GO_*/random_positive_train_max1000_cav/concept_v1.npy" \
        --mf-list $LISTS/mf_ids.txt --bp-list $LISTS/bp_ids.txt --cc-list $LISTS/cc_ids.txt \
        --ec-enzyme-dat ../ec_eval/static/ec_annotations/enzyme.dat \
        --ec-enzclass ../ec_eval/static/ec_annotations/enzclass.txt \
        --go-obo static/go_annotations/go.obo \
        --dims 2 --out results/figures_$TAG/combined_umap_2d.html \
        --coords-csv-out "$STAGE/go_ec_combined_umap.csv" ;;
*) echo "unknown step $step" >&2; exit 1 ;;
esac; done
echo "done -> $STAGE"
