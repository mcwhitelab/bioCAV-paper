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

f=fastas/uniprot_human_all.fasta
o=/xdisk/clairemcwhite/clairemcwhite/fastas_pca/uniprot_human_all.fasta.pkl
pca=reference_population/global_pca_v1.converted.pkl
#  select layer 26, as in arxiv paper
#python $HF_EMBED_SCRIPT -f $f -o $o -ss mean -a -l  -11    -m $ESM_MODEL_DIR -b 1 --max_length 1024
echo $pca
python -c "import pickle; pickle.load(open('$pca', 'rb')); print('ok')"


python $HF_EMBED_SCRIPT -f $f -o $o -ss mean -a -l  -11    -m $ESM_MODEL_DIR -b 1 --max_length 1024 --aa_trim max_length  --aa_pcamatrix_pkl $pca

#=fastas_maxpooling/uniprot_human_all.fasta

#select layer 26, as in arxiv paper
#python $HF_EMBED_SCRIPT -f $f -o ${f}.pkl -ss maxpool -s -l  -11    -m $ESM_MODEL_DIR -b 1 --max_length 2048
 
