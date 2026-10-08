#!/bin/bash
#SBATCH --job-name=mixedlm_de
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --error=slurm_logs/%x_%A_%a.err
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --array=1-29

mkdir -p slurm_logs

CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" COMMANDS.sh)
echo "Running: ${CMD}"
eval "${CMD}"
