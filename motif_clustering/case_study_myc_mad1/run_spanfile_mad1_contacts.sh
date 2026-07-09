#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=span_cav
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL

# -------------
### Inputs — replace these placeholders before running
# -------------
# SPAN_FILE   : tab-separated spans file defining positive examples
#               (accession [TAB pos] [TAB end] or accession TAB p1,p2,...)
# INPUT_FASTA : fasta containing the sequences referenced in SPAN_FILE
# SEARCH_FASTA: fasta of sequences to rank against the trained CAV

SPAN_FILE=fastas/pairs/myc_mad/MAD1_HUMAN.fasta.spans
INPUT_FASTA=fastas/pairs/myc_mad/MAD1_HUMAN.fasta
SEARCH_FASTA=/xdisk/clairemcwhite/clairemcwhite/uniprot_human_all.fasta

# -------------
### Paths
# -------------
source ~/.bashrc
# $0 is unreliable under sbatch (script runs from a spool copy); SLURM sets
# SLURM_SUBMIT_DIR to the directory sbatch was invoked from instead.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../../config/paths.sh"
conda activate /groups/clairemcwhite/envs/$CONDA_ENV
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

ref_neg=reference_population/neg_10000.pkl
ref_pca=reference_population/global_pca_v1.pkl

# Derive output names from inputs
concept_dir=${SPAN_FILE%.spans}_concept
search_pkl=${SEARCH_FASTA}.pkl
out_tsv=${SPAN_FILE%.spans}_query_results.tsv

# -------------
### Step 1: Embed input fasta
# -------------
if [ -f "${INPUT_FASTA}.pkl" ]; then
    echo "Skipping embedding: ${INPUT_FASTA}.pkl  already exists"
else
python $HF_EMBED_SCRIPT \
    -f $INPUT_FASTA \
    -o ${INPUT_FASTA}.pkl \
    --get_aa_embeddings \
    --get_sequence_embedding \
    --strat mean \
    -l -11 \
    -m $ESM_MODEL_DIR \
    -b 1 \
    --max_length 2048
fi
# -------------
### Step 2: Pool embeddings by span → positive embedding matrix
# -------------
python $BIOCAV_REPO/scripts/prepare_embeddings.py \
    --pkl  ${INPUT_FASTA}.pkl \
    --info ${INPUT_FASTA}.pkl.seqnames \
    --spans $SPAN_FILE \
    --out  ${SPAN_FILE%.spans}_pos.npy

# -------------
### Step 3: Train CAV
# -------------
python $BIOCAV_REPO/scripts/train_cav_from_embeddings.py \
    --pos    ${SPAN_FILE%.spans}_pos.npy \
    --neg    $ref_neg \
    --out    $concept_dir \
    --cv-folds 0 \
    --pca-pkl  $ref_pca

# -------------
### Step 4: Embed search fasta
# -------------

# -------------
### Step 4: Embed search fasta (skipped if pkl already exists)
# -------------
if [ -f "$search_pkl" ]; then
    echo "Skipping embedding: $search_pkl already exists"
else
    python $HF_EMBED_SCRIPT \
        -f $SEARCH_FASTA \
        -o $search_pkl \
        --get_sequence_embedding \
        --strat mean \
        -l -11 \
        -m $ESM_MODEL_DIR \
        -b 1 \
        --max_length 2048
fi

# -------------
### Step 5: Query search fasta against CAV
# -------------
#echo $old
#python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
#    --cav  $concept_dir \
#    --pkl  $search_pkl \
#    --out  $out_tsv \
#    --k -1


python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav   $concept_dir \
    --pkl   $search_pkl \
    --out   $out_tsv \
    --k     -1 \
    --sliding-window \
    --spans $SPAN_FILE \
    --fasta $SEARCH_FASTA


#echo $new
#python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
#    --cav  $concept_dir \
#    --pkl  $search_pkl \
#    --out  $out_tsv \
#    --k -1 \
#    --length-correct \
#    --span-length 54 \
#    --fasta $SEARCH_FASTA \
#    --ref-pkl $ref_neg
