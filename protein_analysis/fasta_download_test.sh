#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=embed_fasta
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:02:00
#SBATCH --mail-type=ALL

# -------------
### Code
# -------------
source ~/.bashrc
# $0 is unreliable under sbatch (script runs from a spool copy); SLURM sets
# SLURM_SUBMIT_DIR to the directory sbatch was invoked from instead.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"

conda activate /groups/clairemcwhite/envs/$CONDA_ENV



f=static/GO_sets/GO_0005887/temporal_negative_train_ids.list

s=1000
python scripts/download_fasta.py $f -o ${f}.${s}.fasta --subsample $s --random --seed 42
