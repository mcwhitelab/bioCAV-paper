#!/bin/bash
#SBATCH --job-name=go_cavs_eval
#SBATCH --output=slurm_logs/eval_cc.out
#SBATCH --error=slurm_logs/eval_cc.err
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL


source ~/.bashrc
# $0 is unreliable under sbatch (script runs from a spool copy); SLURM sets
# SLURM_SUBMIT_DIR to the directory sbatch was invoked from instead.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../../config/paths.sh"
conda activate /groups/clairemcwhite/envs/$CONDA_ENV/

for ont in bp cc mf;
do
#python $BIOCAV_REPO/specific_scripts/eval_go_temporal.py \
#    --gold-standard  static/go_annotations/${ont}.gpad \
#    --go-base-dirs   "$PROFAB_GO_SETS/go_dataset_part*" \
#    --out-dir        results/temporal_eval_${ont}/ \
#    --eval-embed-dir results/temporal_eval_${ont}/eval_embeddings/ \
#    --scaler-pkl     reference_population/scaler_v1.pkl \
#    --min-n          5 \
#    --date-cutoff    20211231

python $BIOCAV_REPO/specific_scripts/summarize_temporal_eval.py \
    --results  results/temporal_eval_${ont}/eval_temporal_results.tsv \
    --go-obo static/go_annotations/go.obo \
    --out-dir  results/temporal_eval_${ont}/


done





#
