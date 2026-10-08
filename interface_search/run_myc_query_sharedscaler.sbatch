#!/bin/bash
# Re-run only step 5 (proteome search) of run_spanfile_myc_contacts.sh with the
# July-10 shared-scaler MYC CAV (residues 366-420).  The original job 22850716
# trained this CAV, then OOM'd loading the 94 GB proteome pkl at 64G, so the
# existing MYC_HUMAN.fasta_query_results.tsv is still the May-6 PCA-128 CAV's
# output.  Writes to a new file so both can be compared.
#SBATCH --job-name=myc_query
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=6:00:00
#SBATCH --output=slurm-myc_query_%j.out

source ~/.bashrc
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"
conda activate $CONDA_ENV_DIR
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
cd "$SCRIPT_DIR"

SPAN_FILE=fastas/pairs/myc_mad/MYC_HUMAN.fasta.spans
SEARCH_FASTA=/xdisk/clairemcwhite/clairemcwhite/uniprot_human_all.fasta
python $BIOCAV_REPO/specific_scripts/query_proteins_by_cav.py \
    --cav   ${SPAN_FILE%.spans}_concept \
    --pkl   ${SEARCH_FASTA}.pkl \
    --out   ${SPAN_FILE%.spans}_query_results_sharedscaler.tsv \
    --k     -1 \
    --sliding-window \
    --spans $SPAN_FILE \
    --fasta $SEARCH_FASTA
