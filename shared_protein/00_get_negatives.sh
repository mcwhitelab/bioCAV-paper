#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=embed_fasta
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --gres=gpu:1
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

f=$REFERENCE_NEG_FASTA

#  select layer 26, as in arxiv paper
python $HF_EMBED_SCRIPT -f $f -o reference_population/neg_10000.pkl -ss mean -s -l  -11    -m $ESM_MODEL_DIR -b 1 --max_length 2048

#python $HF_EMBED_SCRIPT -f $f -o reference_population_maxpooling/neg_10000.pkl -ss maxpool -s -l  -11    -m $ESM_MODEL_DIR -b 1 --max_length 2048

