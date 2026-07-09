export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
source "$(dirname "$0")/../config/paths.sh"

LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d

python $BIOCAV_REPO/specific_scripts/de_pseudobulk.py \
    --h5ad          $LIB/data/cells.h5ad \
    --lib-dir       $LIB/ \
    --group-col     cell_type \
    --context-col   tissue \
    --condition-col disease \
    --method        deseq2 \
    --out-dir       $LIB/results/de_vs_cav/ \
    --cav-dir  $LIB/results/gene_correlation/  \
    --donor-col donor_id  \
    --compare-only \
       --lfc-threshold 1.5 \
    --label-top-n 15 \
    --highlight-genes PDCD1 LAG3 HAVCR2 TIGIT TOX CTLA4 \
                      FAP ACTA2 POSTN \
                      VEGFA CD274 MKI67 TGFB1
