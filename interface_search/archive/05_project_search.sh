source "$(dirname "$0")/../../config/paths.sh"
data_dir=$BIOCAV_PAPER_ROOT/protein_analysis/fastas/pairs/myc_mad/

# Step 1: project HLH out of MYC
python $BIOCAV_REPO/specific_scripts/project_away_cav.py \
    --target $data_dir/MYC_concept/ \
    --away   $BIOCAV_PAPER_ROOT/protein_analysis/fastas/PF00010/PF00010_concept \
    --out    $data_dir/MYC_minus_HLH/

# Step 2: query proteins with the residual
python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav  $data_dir/MYC_minus_HLH/ \
    --pkl  fastas/uniprot_human_all.fasta.pkl \
    --out  $data_dir/top_proteins_MYC_minus_HLH.tsv \
    --k -1

# Step 1: project HLH out of MYC
python $BIOCAV_REPO/specific_scripts/project_away_cav.py \
    --target $data_dir/MAD1_concept/ \
    --away   $BIOCAV_PAPER_ROOT/protein_analysis/fastas/PF00010/PF00010_concept \
    --out    $data_dir/MAD1_minus_HLH/

# Step 2: query proteins with the residual
python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav  $data_dir/MAD1_minus_HLH/ \
    --pkl  fastas/uniprot_human_all.fasta.pkl \
    --out  $data_dir/top_proteins_MAD1_minus_HLH.tsv \
    --k -1

# emxphase with HLH
python $BIOCAV_REPO/specific_scripts/combine_cav_scores.py \
    --specific  $data_dir/top_proteins_MAD1.tsv \
    --filter    fastas/PF00010/top_proteins_PF00010.tsv \
    --out       $data_dir/top_proteins_MAD1_x_bHLH.tsv \
    --threshold 0.0 \
    --beta      500.0
python $BIOCAV_REPO/specific_scripts/combine_cav_scores.py \
    --specific  $data_dir/top_proteins_MYC.tsv \
    --filter    fastas/PF00010/top_proteins_PF00010.tsv \
    --out       $data_dir/top_proteins_MYC_x_bHLH.tsv \
    --threshold 0.0 \
    --beta      500.0

