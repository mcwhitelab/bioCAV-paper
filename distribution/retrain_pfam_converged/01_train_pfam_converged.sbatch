#!/bin/bash
# -----------------------------------------------------------------------------
# Retrain all 27,472 Pfam CAVs to lbfgs convergence (tol=1e-10) with --scaler $SCALER
# (shared|own) scaler frame.  CPU-only; embeddings read in place.
# Output: /xdisk/clairemcwhite/clairemcwhite/pfam_cavs_converged_{shared,own}scaler/cavs/<PF>/
# Submit: sbatch -J pfam_conv_<arm> --export=ALL,SCALER=<shared|own>[,CVAL=1e-3] 01_train_pfam_converged.sbatch
# The original library and earlier retrains are not touched.  Idempotent.
# -----------------------------------------------------------------------------
#SBATCH --job-name=pfam_cav_conv
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-39
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --error=slurm_logs/%x_%A_%a.err

# NB: no `set -u` — /etc/bashrc references an unbound BASHRCSOURCED.
set -o pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
/groups/clairemcwhite/envs/core_pkgs2/bin/python train_pfam_converged.py \
    --shard "$SLURM_ARRAY_TASK_ID" --n-shards 40 --scaler "${SCALER:-shared}" --C "${CVAL:-1}"
