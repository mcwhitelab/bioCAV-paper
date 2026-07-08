
source "$(dirname "$0")/../config/paths.sh"

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
python $BIOCAV_REPO/specific_scripts/build_reference_population.py \
    --n-cells   10000 \
    --out       reference_population/ \
    --model      $GENEFORMER_MODEL_DIR/Geneformer-V2-104M \
    --token-dict $GENEFORMER_MODEL_DIR/geneformer/token_dictionary_gc104M.pkl \

