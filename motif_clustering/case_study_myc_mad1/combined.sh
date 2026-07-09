#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=embed_fasta
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL

# -------------
### Code
# -------------
source ~/.bashrc
# $0 is unreliable under sbatch (script runs from a spool copy); SLURM sets
# SLURM_SUBMIT_DIR to the directory sbatch was invoked from instead.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../../config/paths.sh"

conda activate /groups/clairemcwhite/envs/$CONDA_ENV

for f in fastas/pairs/myc_mad/*fasta;
do

#  select layer 26, as in arxiv paper
python $HF_EMBED_SCRIPT -f $f -o ${f}.pkl --get_aa_embeddings --get_sequence_embedding --strat mean -l  -11    -m $ESM_MODEL_DIR -b 1 --max_length 2048
done


#Prepare pooled embedding vectors from a .pkl file for CAV training.

#Reads a .pkl + .pkl.info embedding file and a spans file, pools each sample's
#embeddings according to its span, and saves the result as a .npy file.

#Run this separately for positive and negative sets, then pass the resulting
#.npy files to train_cav_from_embeddings.py.

#Spans file format (tab-separated, no header, '#' lines are comments):
#    accession                       # whole-sequence
#    accession  TAB  pos             # single position (e.g. PTM)
#    accession  TAB  start  TAB  end # window [start, end)

#Examples
#--------
# Motif positives (window spans):
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
workdir=$BIOCAV_PAPER_ROOT/protein_analysis
for f in $workdir/fastas/pairs/myc_mad/*fasta
do
echo $f
python $BIOCAV_REPO/scripts/prepare_embeddings.py \
    --pkl ${f}.pkl \
    --info ${f}.pkl.seqnames \
    --spans ${f}.spans \
    --out ${f}.spans.npy
done


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

pkl=fastas/uniprot_human_all.fasta.pkl
python $BIOCAV_REPO/specific_scripts/build_protein_index.py \
    --pkl  $pkl \
    --out  ${pkl}.index

