#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=query_tax399
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
# Query the RLK5, FLS2, PSKR1, and multi-example SERK-interface CAVs against
# fastas/uniprotkb_IPR001611_AND_taxonomy_id_399_2026_07_10.fasta (in place
# of the ORYSJ_LRR.fasta search used previously), so results are comparable
# across all four CAVs. Assumes the CAVs are already trained
# (fastas/pairs/lrr/*_concept/ directories exist). The search fasta's .pkl
# already exists, so no re-embedding is needed.

SEARCH_FASTA=fastas/uniprotkb_IPR001611_AND_taxonomy_id_399_2026_07_10.fasta

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
    local out_tsv=${span_file%.span}_query_taxonomy399_results.tsv

    python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
        --cav   $concept_dir \
        --pkl   $search_pkl \
        --out   $out_tsv \
        --k     -1 \
        --sliding-window \
        --spans $span_file \
        --fasta $SEARCH_FASTA
}

query fastas/pairs/lrr/RLK5_ARATH.fasta.span                fastas/pairs/lrr/RLK5_ARATH.fasta_concept
query fastas/pairs/lrr/FLS2_ARATH.fasta.span                fastas/pairs/lrr/FLS2_ARATH.fasta_concept
query fastas/pairs/lrr/PSKR1_ARATH.fasta.span                fastas/pairs/lrr/PSKR1_ARATH.fasta_concept
query fastas/pairs/lrr/SERK_interface_multi.fasta.span       fastas/pairs/lrr/SERK_interface_multi.fasta_concept
