# Archived: exploratory whole-PF00010-domain CAV, superseded by the direct
# MYC/MAD1 concept training in interface_search/. Not used to produce the
# motif_clustering CAVs (those come from a separate, grad-student-owned
# pipeline) and not part of the interface_search case study.
source "$(dirname "$0")/../../../config/paths.sh"
data_dir=$BIOCAV_PAPER_ROOT/protein_analysis/fastas/PF00010/


#python $BIOCAV_REPO/scripts/train_cav_from_embeddings.py \
#    --pos $data_dir/MYC_HUMAN.fasta.pkl \
#    --neg  reference_population/neg_10000.pkl \
#    --out $data_dir/MYC_concept \
#    --cv-folds 0


## 3. Validate — localization (requires aa_embeddings in pkl)
#python $BIOCAV_REPO/scripts/run_attribution.py \
#    --pkl $data_dir/MAD1_HUMAN.fasta.pkl \
#    --cav-dir $data_dir/MYC_concept \
#    --out $data_dir/validation_MYC_to_MAD1.tsv

python $BIOCAV_REPO/scripts/train_cav_from_embeddings.py \
    --pos $data_dir/pos_1000.fasta.pkl \
    --neg  reference_population/neg_10000.pkl \
    --out $data_dir/PF00010_concept \
    --cv-folds 5 \
    --pca-pkl reference_population/global_pca_v1.pkl



#python $BIOCAV_REPO/scripts/scan_sequence.py \
#    --pkl $data_dir/temporal_positive_validation.fasta.pkl \
#    --cav-dir ./cavs/$go/ \
#    --window-size 20 \
#    --out ./cavs/GO_0005887/validation_scan.tsv
