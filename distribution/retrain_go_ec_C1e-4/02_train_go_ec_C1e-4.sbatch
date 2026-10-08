#!/bin/bash
# GO/EC retrain at C=1e-4, step 2 (CPU): train all 7,631 GO/EC CAVs from the pooled
# unique-sequence embeddings of step 1.  Submitted with --dependency=afterok:<step 1>.
# Output: /xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_sharedscaler/  Idempotent.
#SBATCH --job-name=goec_train_C1e-4
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=0-39
#SBATCH --output=slurm_logs/train_%A_%a.out
#SBATCH --error=slurm_logs/train_%A_%a.err
set -o pipefail
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
/groups/clairemcwhite/envs/core_pkgs4/bin/python 02_train_go_ec_C1e-4.py --shard "$SLURM_ARRAY_TASK_ID" --n-shards 40
