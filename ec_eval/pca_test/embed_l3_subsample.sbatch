#!/bin/bash
# PCA test (2026-09-29): embed up to 100 sampled training positives per EC level-3
# CAV, with the exact flags of bioCAV train_cav_from_span.sh step 2.
#SBATCH --job-name=ec_l3_embed
#SBATCH --output=slurm_logs/ec_l3_embed_%j.out
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --gres=gpu:1
set -o pipefail
source ~/.bashrc
conda activate /groups/clairemcwhite/envs/core_pkgs4
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
python /groups/clairemcwhite/claire_workspace/github/mcwlab_utils/hf_embed_new.py \
    -f l3_subsample.fasta -o l3_subsample.fasta.pkl \
    --get_aa_embeddings --get_sequence_embedding --strat mean -l -11 \
    -m /groups/clairemcwhite/models/ESMplusplus_large -b 1 --max_length 2048
