source "$(dirname "$0")/../config/paths.sh"

# Re-download the current GO consortium EC<->GO mapping (comment out to reuse
# the checked-in snapshot in static/ec2go.txt):
# curl -s -o "$(dirname "$0")/static/ec2go.txt" https://current.geneontology.org/ontology/external2go/ec2go

$CONDA_ENV_DIR/bin/python "$(dirname "$0")/go_ec_cosine_validation.py" \
    --ec2go        "$(dirname "$0")/static/ec2go.txt" \
    --go-base-dirs "$PROFAB_GO_SETS/go_dataset_part*" \
    --ec-base-dirs "$PROFAB_EC_SETS/ec_dataset_part*" \
    --out-dir      "$(dirname "$0")/results"

$CONDA_ENV_DIR/bin/python "$(dirname "$0")/plot_go_ec_cosine_validation.py" \
    --results-dir "$(dirname "$0")/results"

# Supplemental table: annotated EC partner's training-set size by rank
# bucket -- shows the rare rank>10 failures are the most severely
# under-trained EC CAVs, not a biological mismatch.
$CONDA_ENV_DIR/bin/python "$(dirname "$0")/make_rank_vs_training_size_table.py" \
    --pairs-csv    "$(dirname "$0")/results/go_ec_cosine_validation.csv" \
    --ec-base-dirs "$PROFAB_EC_SETS/ec_dataset_part*" \
    --out-dir      "$(dirname "$0")/results"

$CONDA_ENV_DIR/bin/python "$(dirname "$0")/render_supp_table.py" \
    --csv "$(dirname "$0")/results/supp_table_rank_vs_training_size.csv"
