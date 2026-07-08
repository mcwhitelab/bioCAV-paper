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

conda activate /groups/clairemcwhite/envs/core_pkgs4

for f in fastas/PF00010/*fasta;
do

#  select layer 26, as in arxiv paper
python /groups/clairemcwhite/claire_workspace/github/mcwlab_utils/hf_embed_new.py -f $f -o ${f}.pkl --get_aa_embeddings --get_sequence_embedding --strat mean -l  -11    -m /groups/clairemcwhite/models/ESMplusplus_large -b 1 --max_length 2048 
done


