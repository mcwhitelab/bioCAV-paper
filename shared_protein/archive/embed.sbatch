#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=embed_fasta
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL

# -------------
### Code
# -------------
source ~/.bashrc


conda activate /groups/clairemcwhite/envs/core_pkgs


#GO=GO_0005887
GO=GO_0005886
for f in static/GO_sets/$GO/*t[re][as]*fasta; do
    python ../github/mcwlab_utils/hf_embed_new.py -f $f -o ${f}.pkl -ss mean -s -l  -9 -10 -11 -12   -m /groups/clairemcwhite/models/ESMplusplus_large -b 1 --max_length 2048 --subsample 200
done

for f in static/GO_sets/$GO/*valid*fasta; do
    python ../github/mcwlab_utils/hf_embed_new.py -f $f -o ${f}.pkl -ss mean -a -l  -9 -10 -11 -12   -m /groups/clairemcwhite/models/ESMplusplus_large -b 1 --max_length 2048 --subsample 200
done

