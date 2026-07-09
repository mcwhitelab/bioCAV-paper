source "$(dirname "$0")/../../config/paths.sh"
data_dir=$BIOCAV_PAPER_ROOT/protein_analysis/fastas/pairs/myc_mad/


python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav  $data_dir/MYC_concept \
    --pkl  fastas/uniprot_human_all.fasta.pkl \
    --out  $data_dir/top_proteins_MYC.tsv \
    --k -1

python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav  $data_dir/MAD1_concept \
    --pkl  fastas/uniprot_human_all.fasta.pkl \
    --out  $data_dir/top_proteins_MAD1.tsv \
    --k -1

