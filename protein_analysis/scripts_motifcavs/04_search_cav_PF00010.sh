source "$(dirname "$0")/../../config/paths.sh"
data_dir=$BIOCAV_PAPER_ROOT/protein_analysis/fastas/PF00010/


python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav  $data_dir/PF00010_concept \
    --pkl  fastas/uniprot_human_all.fasta.pkl \
    --out  $data_dir/top_proteins_PF00010.tsv \
    --k -1


