#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=embed_fasta
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
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

for f in fastas/PF00010/*fasta;
do

#  select layer 26, as in arxiv paper
python $HF_EMBED_SCRIPT -f $f -o ${f}.pkl --get_aa_embeddings --get_sequence_embedding --strat mean -l  -11    -m $ESM_MODEL_DIR -b 1 --max_length 2048
done


