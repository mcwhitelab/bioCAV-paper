#!/bin/bash
# Pfam CAVs, MYC PCA recipe k=128 with the SHARED scaler (control arm); see train_pfam_pca128.py.
# Output: /xdisk/clairemcwhite/clairemcwhite/pfam_cavs_pca128_sharedscaler/cavs/<PF>/ . Idempotent.
# -----------------------------------------------------------------------------
#SBATCH --job-name=pfam_cav_pca128sh
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --array=0-39
#SBATCH --output=slurm_logs/pca128sh_%A_%a.out
#SBATCH --error=slurm_logs/pca128sh_%A_%a.err

# NB: no `set -u` — /etc/bashrc references an unbound BASHRCSOURCED.
set -o pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
/groups/clairemcwhite/envs/core_pkgs2/bin/python train_pfam_pca128.py \
    --shard "$SLURM_ARRAY_TASK_ID" --n-shards 40 --scaler shared
