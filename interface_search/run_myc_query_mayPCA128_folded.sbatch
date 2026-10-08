#!/bin/bash
# Diagnostic: score the May PCA-128 MYC CAV (folded to 1152-d, ranking-identical)
# on the July proteome embedding, to separate the PCA effect from the re-embed.
#SBATCH --job-name=myc_query_pcafold
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=6:00:00
#SBATCH --output=slurm-myc_query_pcafold_pcafold_%j.out

source ~/.bashrc
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"
conda activate $CONDA_ENV_DIR
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
cd "$SCRIPT_DIR"

SPAN_FILE=fastas/pairs/myc_mad/MYC_HUMAN.fasta.spans
SEARCH_FASTA=/xdisk/clairemcwhite/clairemcwhite/uniprot_human_all.fasta
python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav   ${SPAN_FILE%.spans}_concept_mayPCA128_folded \
    --pkl   ${SEARCH_FASTA}.pkl \
    --out   ${SPAN_FILE%.spans}_query_results_mayPCA128_folded.tsv \
    --k     -1 \
    --sliding-window \
    --spans $SPAN_FILE \
    --fasta $SEARCH_FASTA
