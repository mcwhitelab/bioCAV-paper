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
# AWAY_FASTA  : fasta for a second interface CAV whose common signal with
#               RLK5-SERK should be projected away (CLV1-CLE, another
#               LRR-receptor-kinase interface) — leaves the RLK5-SERK-specific
#               residual.

SPAN_FILE=fastas/pairs/lrr/RLK5_ARATH.fasta.span
INPUT_FASTA=fastas/pairs/lrr/RLK5_ARATH.fasta
SEARCH_FASTA=fastas/uniprotkb_IPR001611_AND_model_organism_2026_07_10.fasta
AWAY_FASTA=fastas/pairs/lrr/CLV1_ARATH.fasta
AWAY_SPAN_FILE=${AWAY_FASTA}.span

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
away_concept_dir=${AWAY_SPAN_FILE%.span}_concept
residual_dir=${concept_dir}_minus_CLV1
residual_out_tsv=${SPAN_FILE%.span}_minus_CLV1_query_results.tsv

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

# -------------
### Step 6: Train away CAV (CLV1-CLE interface)
# -------------
# Mirrors Steps 1-3 above, but for the interface direction being projected away.

if [ -f "${AWAY_FASTA}.pkl" ]; then
    echo "Skipping embedding: ${AWAY_FASTA}.pkl already exists"
else
python $HF_EMBED_SCRIPT \
    -f $AWAY_FASTA \
    -o ${AWAY_FASTA}.pkl \
    --get_aa_embeddings \
    --get_sequence_embedding \
    --strat mean \
    -l -11 \
    -m $ESM_MODEL_DIR \
    -b 1 \
    --max_length 2048
fi

python $BIOCAV_REPO/scripts/prepare_embeddings.py \
    --pkl  ${AWAY_FASTA}.pkl \
    --info ${AWAY_FASTA}.pkl.seqnames \
    --spans $AWAY_SPAN_FILE \
    --out  ${AWAY_SPAN_FILE%.span}_pos.npy

python $BIOCAV_REPO/scripts/train_cav_from_embeddings.py \
    --pos    ${AWAY_SPAN_FILE%.span}_pos.npy \
    --neg    $ref_neg \
    --out    $away_concept_dir \
    --cv-folds 0 \
    --scaler-pkl $ref_scaler

# -------------
### Step 7: Project away CLV1-CLE signal from the RLK5-SERK CAV
# -------------
# Leaves the RLK5-SERK-specific residual (common LRR-receptor-kinase signal removed).
# See bioCAV/specific_scripts/project_away_cav.py for details.

python $BIOCAV_REPO/specific_scripts/project_away_cav.py \
    --target $concept_dir \
    --away   $away_concept_dir \
    --out    $residual_dir

python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav   $residual_dir \
    --pkl   $search_pkl \
    --out   $residual_out_tsv \
    --k     -1 \
    --sliding-window \
    --spans $SPAN_FILE \
    --fasta $SEARCH_FASTA
