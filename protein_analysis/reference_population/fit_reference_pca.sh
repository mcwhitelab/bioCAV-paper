source "$(dirname "$0")/../../config/paths.sh"

python $BIOCAV_REPO/scripts/fit_global_pca.py \
    --pkl  reference_population/neg_10000.pkl \
    --out  reference_population/global_pca_v1.pkl
