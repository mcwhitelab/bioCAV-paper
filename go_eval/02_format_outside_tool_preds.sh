source "$(dirname "$0")/../config/paths.sh"

# Format API jsonl
dgw_dir=$AHMAD_GO_RESULTS_DIR
for go in bp cc mf;
do
python scripts/format_deepgoweb_api.py $dgw_dir/val_proteins_${go}.deepgo.jsonl outside_tools/val_proteins_${go}.deepgo.csv
done

# Format DeepGOSE
#
# Run on the model first (external, produces the negative-test preds below):
#   sbatch /xdisk/clairemcwhite/test_deepgose_val.sbatch    # val proteins
#   sbatch /xdisk/clairemcwhite/test_deepgose.sbatch        # negative test proteins

for f in $PROFAB_GO_SETS/go_dataset_part*/*/*gz; do gunzip "$f"; done

# All the deepgose preds on negative test proteins
cat $PROFAB_GO_SETS/go_dataset_part*/*/random_negative*.span_preds_bp.tsv | sort -u > outside_tools/neg_test_deepgose_preds_bp.tsv
cat $PROFAB_GO_SETS/go_dataset_part*/*/random_negative*.span_preds_mf.tsv | sort -u > outside_tools/neg_test_deepgose_preds_mf.tsv
cat $PROFAB_GO_SETS/go_dataset_part*/*/random_negative*.span_preds_cc.tsv | sort -u > outside_tools/neg_test_deepgose_preds_cc.tsv

# Combine positives and negatives
cat outside_tools/neg_test_deepgose_preds_bp.tsv outside_tools/val_proteins.span_preds_bp.tsv > outside_tools/deepgose_preds_bp.tsv
cat outside_tools/neg_test_deepgose_preds_mf.tsv outside_tools/val_proteins.span_preds_mf.tsv > outside_tools/deepgose_preds_mf.tsv
cat outside_tools/neg_test_deepgose_preds_cc.tsv outside_tools/val_proteins.span_preds_cc.tsv > outside_tools/deepgose_preds_cc.tsv
