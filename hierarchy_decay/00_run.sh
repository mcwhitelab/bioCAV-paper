source "$(dirname "$0")/../config/paths.sh"
cd "$(dirname "$0")/.."

$CONDA_ENV_DIR/bin/python hierarchy_decay/ec_hierarchy_decay.py \
    --val-pkl      ec_eval/static/ec_annotations/EC-Bench/test_ec.fasta.pkl \
    --gold         ec_eval/static/ec_annotations/EC-Bench/test_ec_long.tsv \
    --ec-base-dirs "$PROFAB_EC_SETS/ec_dataset_part*" \
    --scaler-pkl   shared_protein/reference_population/scaler_v1.pkl \
# 2026-09-29: A0A452E9Y6 (1.11.1.7) and P0DPS1 (2.7.7.7) replaced because they are
# training positives of CAVs in their own lineage (1.11/1.11.1/1.11.1.7; 2.7.7.7).
# All five below are leak-free at every level: in neither the train nor the test span of
# any lineage CAV. Rule: full 1-4 CAV lineage, leak-free, one per level-2 EC class
# where possible, first accession alphabetically. Class 1 has only C0HLQ6 leak-free.
    --proteins     C0HLQ6:1.2.1.12 P0DQP1:2.1.1.72 P0DUB9:3.1.1.3 \
                   P0DUE1:3.2.2.6 A0A345BJN1:4.3.1.24 \
    --out-dir      hierarchy_decay/results

$CONDA_ENV_DIR/bin/python hierarchy_decay/plot_hierarchy_decay.py \
    --ec-csv hierarchy_decay/results/ec_hierarchy_decay.csv \
    --out    hierarchy_decay/results/hierarchy_decay.png

# GO analysis (parallel to EC, using is_a ancestor chains instead of dot
# notation) is available but not part of the current figure -- rerun with:
#   $CONDA_ENV_DIR/bin/python hierarchy_decay/go_hierarchy_decay.py \
#       --val-pkl      go_eval/results/temporal_eval_mf/eval_embeddings/val_proteins.span.fasta.pkl \
#       --go-obo       go_eval/static/go_annotations/go.obo \
#       --go-base-dirs "$PROFAB_GO_SETS/go_dataset_part*" \
#       --scaler-pkl   shared_protein/reference_population/scaler_v1.pkl \
#       --proteins     O69729:GO:0000155 A3KPE8:GO:0000175 G8BDQ3:GO:0000400 \
#                      F8VPJ6:GO:0000979 O94913:GO:0000993 \
#       --out-dir      hierarchy_decay/results

# Candidates were selected via:
#   find gold-standard proteins with level-4 EC labels whose full 1-2-3-4
#   lineage all have trained CAVs (see ec_hierarchy_decay.py's inline logic)
#   python hierarchy_decay/find_go_chains.py \
#       --go-obo go_eval/static/go_annotations/go.obo \
#       --pairs  figures/figure_data/go_specificity_ranks_mf.csv \
#       --go-base-dirs "$PROFAB_GO_SETS/go_dataset_part*"
