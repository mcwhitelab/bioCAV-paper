#!/bin/bash
#SBATCH --job-name=ec_cavs
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --error=slurm_logs/%x_%A_%a.err
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --array=1-48%20
#SBATCH --mail-type=ALL

mkdir -p slurm_logs

CHUNK_SIZE=20
START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))

for i in $(seq $START $END); do
    CMD=$(sed -n "${i}p" COMMANDS_files/cav_EC_COMMANDS.sh)
    [ -z "$CMD" ] && continue
    echo "Running line ${i}: ${CMD}"
    eval "${CMD}"
done

##TCH --partition=gpu_standard
##TCH --gres=gpu:1

