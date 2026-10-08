#!/bin/bash
# -----------------------------------------------------------------------------
# Retrain all 27,472 Pfam CAVs in the shared reference scaler frame.
# CPU-only: positives/negatives are already embedded and are read in place from
# /xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/embeddings.
# Output: /xdisk/clairemcwhite/clairemcwhite/pfam_cavs_shared_scaler/cavs/<PF>/
# The original library is not touched.  Idempotent: resubmit to pick up failures.
# Array is 0-based because the script shards with motifs[shard::n_shards].
# -----------------------------------------------------------------------------
#SBATCH --job-name=pfam_cav_shared
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-39
#SBATCH --output=slurm_logs/retrain_%A_%a.out
#SBATCH --error=slurm_logs/retrain_%A_%a.err

# NB: no `set -u` — /etc/bashrc references an unbound BASHRCSOURCED.
set -o pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
/groups/clairemcwhite/envs/core_pkgs2/bin/python retrain_pfam_cavs.py \
    --shard "$SLURM_ARRAY_TASK_ID" --n-shards 40
