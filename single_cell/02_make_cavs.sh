#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=embed_fasta
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mem=128G

# -------------
### Code
# -------------
source ~/.bashrc
# $0 is unreliable under sbatch (script runs from a spool copy); SLURM sets
# SLURM_SUBMIT_DIR to the directory sbatch was invoked from instead.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"

conda activate /groups/clairemcwhite/envs/$CONDA_ENV



#ds=9d8e5dca-03a3-457d-b7fb-844c75735c83

ds=b617ee1b-f8c8-4de9-b82b-e803ab93550d
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
python $BIOCAV_REPO/specific_scripts/run_cav_pipeline.py \
    --library    cav_library/$ds \
    --dataset-id $ds \
    --model      $GENEFORMER_MODEL_DIR/Geneformer-V2-104M \
    --token-dict $GENEFORMER_MODEL_DIR/geneformer/token_dictionary_gc104M.pkl \
    --reference-pkl reference_population/embeddings/cells.pkl \
    --pca-pkl       reference_population/global_pca_v1.pkl \
    --reference-n 10000 \
    --pos-n 500 \
    --from-step embed
