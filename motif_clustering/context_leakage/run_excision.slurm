#!/bin/bash
#SBATCH --account=clairemcwhite
#SBATCH --job-name=cav_excision
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=gpu_standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --gpus=1
#SBATCH --time=04:00:00

cd /groups/clairemcwhite/cav_workspace/bioCAV-paper/motif_clustering/context_leakage

/groups/clairemcwhite/envs/core_pkgs4/bin/python run_excision.py \
    --n-pos 300 --n-neg 1000 --batch-size 16 --out results
