#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=im9_mutants
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL

# -------------
### Purpose
# -------------
# Causal ablation test: does either Im9 CAV (generic or contrastive-vs-
# homologs) respond to the actual specificity-determining side chains, or
# just to whole-sequence identity? Score WT Im9 against three known
# interface point mutants (V34A, D51A, Y55A) and a chimera with Im2's
# residues swapped in at the helix-II specificity triad (L33D/V34N/T38R).
# If the CAVs track causal biophysics, scores should drop for the mutants/
# chimera relative to WT; if they only track whole-sequence identity, the
# score should barely move (these are single/triple substitutions in an
# 86 aa protein).

SEARCH_FASTA=fastas/pairs/colicin/IMM9_mutants.fasta

# -------------
### Paths
# -------------
source ~/.bashrc

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"
conda activate $CONDA_ENV_DIR
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

search_pkl=${SEARCH_FASTA}.pkl

# -------------
### Step 1: Embed mutant fasta
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
### Step 2: Score against both CAVs
# -------------
python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav   fastas/pairs/colicin/IMM9_ECOLX.fasta_concept \
            fastas/pairs/colicin/IMM9_vs_homologs_concept \
    --pkl   $search_pkl \
    --out   fastas/pairs/colicin/IMM9_mutants_query_results.tsv \
    --k     -1
