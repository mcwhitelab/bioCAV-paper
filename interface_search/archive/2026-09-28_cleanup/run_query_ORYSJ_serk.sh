#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=query_orysj
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
# Query the FLS2, PSKR1, and multi-example SERK-interface CAVs against
# fastas/ORYSJ_LRR.fasta — the same rice-LRR search fasta already used for
# RLK5 (see fastas/pairs/lrr/RLK5_ARATH.fasta_query_ORYSJ_results.long.tsv),
# so results are comparable across all four CAVs. Assumes the CAVs are
# already trained (fastas/pairs/lrr/*_concept/ directories exist).

SEARCH_FASTA=fastas/ORYSJ_LRR.fasta

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
### Step 1: Embed search fasta (skipped if pkl already exists)
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
### Step 2: Query each CAV against the search fasta
# -------------
query() {
    local span_file=$1
    local concept_dir=$2
    local out_tsv=${span_file%.span}_query_ORYSJ_results.tsv

    python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
        --cav   $concept_dir \
        --pkl   $search_pkl \
        --out   $out_tsv \
        --k     -1 \
        --sliding-window \
        --spans $span_file \
        --fasta $SEARCH_FASTA
}

query fastas/pairs/lrr/FLS2_ARATH.fasta.span              fastas/pairs/lrr/FLS2_ARATH.fasta_concept
query fastas/pairs/lrr/PSKR1_ARATH.fasta.span              fastas/pairs/lrr/PSKR1_ARATH.fasta_concept
query fastas/pairs/lrr/SERK_interface_multi.fasta.span     fastas/pairs/lrr/SERK_interface_multi.fasta_concept
