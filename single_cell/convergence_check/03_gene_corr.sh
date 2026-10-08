#!/bin/bash
source ~/.bashrc >/dev/null 2>&1; conda activate /groups/clairemcwhite/envs/core_pkgs4
cd /groups/clairemcwhite/cav_workspace/bioCAV-paper/single_cell
source ../config/paths.sh; export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
LIB=cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d
for arm in ${ARMS:-stored orig conv1 conv4}; do
  python scripts/paired_donor_gene_corr.py --h5ad $LIB/data/cells.h5ad \
    --coords convergence_check/coords_${arm}.tsv --lib-dir $LIB \
    --pairs-csv static/paired_donor_pairs.csv \
    --out-dir convergence_check/gene_corr_${arm} \
    --only neutrophil__breast__normal_vs_breast_cancer epithelial_cell__lung10x__normal_vs_lung_cancer fibroblast__colorectum__normal_vs_colorectal_cancer
done
echo DONE
