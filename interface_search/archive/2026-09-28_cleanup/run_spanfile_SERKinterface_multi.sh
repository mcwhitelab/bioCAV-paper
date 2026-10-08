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
### Inputs
# -------------
# SPAN_FILE   : multi-example spans file — one line per LRR-RK/SERK interface
#               (RLK5, PSKR1, FLS2), each a window over its BAK1/SERK-facing
#               ectodomain region. Training on all three positives at once
#               (rather than one span + project_away_cav.py) lets the
#               receptor-specific idiosyncrasies average out, leaving a
#               direction closer to the shared SERK-binding signal itself.
# INPUT_FASTA : fasta containing the sequences referenced in SPAN_FILE
#               (concatenation of RLK5_ARATH.fasta + PSKR1_ARATH.fasta +
#               FLS2_ARATH.fasta)
# SEARCH_FASTA: fasta of sequences to rank against the trained CAV

SPAN_FILE=fastas/pairs/lrr/SERK_interface_multi.fasta.span
INPUT_FASTA=fastas/pairs/lrr/SERK_interface_multi.fasta
SEARCH_FASTA=fastas/uniprotkb_IPR001611_AND_model_organism_2026_07_10.fasta

# -------------
### Paths
# -------------
source ~/.bashrc

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"
conda activate $CONDA_ENV_DIR
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

ref_neg="$BIOCAV_PAPER_ROOT/shared_protein/reference_population/neg_10000.pkl"
ref_scaler="$BIOCAV_PAPER_ROOT/shared_protein/reference_population/scaler_v1.pkl"

# Derive output names from inputs
concept_dir=${SPAN_FILE%.span}_concept
search_pkl=${SEARCH_FASTA}.pkl
out_tsv=${SPAN_FILE%.span}_query_results.tsv

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
### Step 2: Pool embeddings by span → positive embedding matrix (3 rows: RLK5, PSKR1, FLS2)
# -------------
python $BIOCAV_REPO/scripts/prepare_embeddings.py \
    --pkl  ${INPUT_FASTA}.pkl \
    --info ${INPUT_FASTA}.pkl.seqnames \
    --spans $SPAN_FILE \
    --out  ${SPAN_FILE%.span}_pos.npy

# -------------
### Step 3: Train CAV
# -------------
python $BIOCAV_REPO/scripts/train_cav_from_embeddings.py \
    --pos    ${SPAN_FILE%.span}_pos.npy \
    --neg    $ref_neg \
    --out    $concept_dir \
    --cv-folds 0 \
    --scaler-pkl $ref_scaler

# -------------
### Step 4: Embed search fasta (skipped if pkl already exists)
# -------------
if [ -f "$search_pkl" ]; then
    echo "Skipping embedding: $search_pkl already exists"
else
    python $HF_EMBED_SCRIPT \
        -f $SEARCH_FASTA \
        -o $search_pkl \
        --get_aa_embeddings \
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
python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav   $concept_dir \
    --pkl   $search_pkl \
    --out   $out_tsv \
    --k     -1 \
    --sliding-window \
    --spans $SPAN_FILE \
    --fasta $SEARCH_FASTA
