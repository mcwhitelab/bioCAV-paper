source "$(dirname "$0")/../config/paths.sh"

# Validation proteins with temporal split are downloaded manually from
# QuickGO and placed at:
#   results/temporal_eval_{mf,bp,cc}/eval_embeddings/val_proteins.span.fasta

for ont in bp cc mf;
do

python $BIOCAV_REPO/specific_scripts/summarize_temporal_eval.py \
    --results  results/temporal_eval_${ont}/eval_temporal_results.tsv \
    --go-obo static/go_annotations/go.obo \
    --out-dir  results/temporal_eval_${ont}/ \
    --figure-data-dir "$BIOCAV_PAPER_ROOT/figures/figure_data" \
    --label $ont

done
