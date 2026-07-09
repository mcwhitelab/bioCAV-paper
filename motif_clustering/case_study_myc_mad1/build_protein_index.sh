source "$(dirname "$0")/../../config/paths.sh"
pkl=fastas/uniprot_human_all.fasta.pkl
python $BIOCAV_REPO/specific_scripts/build_protein_index.py \
    --pkl  $pkl \
    --out  ${pkl}.index

