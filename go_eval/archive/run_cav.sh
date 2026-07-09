go=$1
source "$(dirname "$0")/../../config/paths.sh"
data_dir=$BIOCAV_PAPER_ROOT/protein_analysis/static/GO_sets/$go
python $BIOCAV_REPO/scripts/train_cav_from_embeddings.py \
    --pos $data_dir/temporal_positive_train.fasta.pkl \
    --neg $data_dir/temporal_negative_train.fasta.pkl \
    --out ./cavs/$go

## 2. Test — classification accuracy
python $BIOCAV_REPO/scripts/evaluate_cav.py \
    --pos $data_dir/temporal_positive_test.fasta.pkl \
    --neg $data_dir/temporal_negative_test.fasta.pkl \
    --cav-dir ./cavs/$go \
    --out ./cavs/$go/test_eval.json

## 3. Validate — localization (requires aa_embeddings in pkl)
python $BIOCAV_REPO/scripts/run_attribution.py \
    --pkl $data_dir/temporal_positive_validation.fasta.pkl \
    --cav-dir ./cavs/$go/ \
    --out ./cavs/$go/validation_localization.tsv


#python $BIOCAV_REPO/scripts/scan_sequence.py \
#    --pkl $data_dir/temporal_positive_validation.fasta.pkl \
#    --cav-dir ./cavs/$go/ \
#    --window-size 20 \
#    --out ./cavs/GO_0005887/validation_scan.tsv
