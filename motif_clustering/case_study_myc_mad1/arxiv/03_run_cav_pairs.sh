source "$(dirname "$0")/../../../config/paths.sh"
data_dir=$BIOCAV_PAPER_ROOT/protein_analysis/fastas/pairs/myc_mad/


python $BIOCAV_REPO/scripts/train_cav_from_embeddings.py \
    --pos $data_dir/MYC_HUMAN.fasta.pkl \
    --neg  reference_population/neg_10000.pkl \
    --out $data_dir/MYC_concept \
    --cv-folds 0 \
    --pca-pkl reference_population/global_pca_v1.pkl


## 3. Validate — localization (requires aa_embeddings in pkl)
python $BIOCAV_REPO/scripts/run_attribution.py \
    --pkl $data_dir/MAD1_HUMAN.fasta.pkl \
    --cav-dir $data_dir/MYC_concept \
    --out $data_dir/validation_MYC_to_MAD1.tsv

python $BIOCAV_REPO/scripts/train_cav_from_embeddings.py \
    --pos $data_dir/MAD1_HUMAN.fasta.pkl \
    --neg  reference_population/neg_10000.pkl \
    --out $data_dir/MAD1_concept \
    --cv-folds 0 \
    --pca-pkl reference_population/global_pca_v1.pkl




## 3. Validate — localization (requires aa_embeddings in pkl)
python $BIOCAV_REPO/scripts/run_attribution.py \
    --pkl $data_dir/MYC_HUMAN.fasta.pkl \
    --cav-dir $data_dir/MAD1_concept \
    --out $data_dir/validation_MAD1_to_MYC.tsv
