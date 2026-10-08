#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=span_cav
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL

# -------------
### Inputs
# -------------
# Colicin E9 DNase / immunity-protein specificity system (Im9 cognate,
# Im2/Im7/Im8 homologous non-cognate). Unlike the LRR receptor ectodomains,
# the ~86-88 aa Im proteins ARE almost entirely the toxin-binding domain, so
# SPAN_FILE uses whole-sequence mode (accession only, no residue window) —
# a sub-window like helix II would hard-code the expected answer rather
# than letting the CAV's own residue-level signal reveal it (see
# run_pdb_bfactor_Im9.sh for the sliding-window structural validation onto
# 1EMV that checks whether signal concentrates on helix II, as reported by
# Keeble/Wallis-style homologue-scanning mutagenesis of Im9 for E9 DNase).
#
# SPAN_FILE   : whole-sequence span for Im9 (positive training example)
# INPUT_FASTA : fasta containing Im9
# SEARCH_FASTA: small demonstration search set — Im9 (self, sanity check),
#               Im2, Im7, Im8 — the homologous immunity-protein family that
#               all bind related E-colicin DNases but with very different
#               affinities for E9 specifically.

SPAN_FILE=fastas/pairs/colicin/IMM9_ECOLX.fasta.span
INPUT_FASTA=fastas/pairs/colicin/IMM9_ECOLX.fasta
SEARCH_FASTA=fastas/pairs/colicin/IMM_family_search.fasta

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
### Step 2: Pool embeddings by span (whole-sequence) → positive embedding matrix
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
### Step 5: Query search fasta against CAV (whole-sequence scoring — no
### sliding window, since these proteins are shorter than any sensible
### window and are already ~entirely the domain of interest)
# -------------
python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav   $concept_dir \
    --pkl   $search_pkl \
    --out   $out_tsv \
    --k     -1
