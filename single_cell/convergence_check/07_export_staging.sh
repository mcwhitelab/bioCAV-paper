#!/bin/bash
# Export Fig 6 sc_*.csv for the gate (stored CAVs) and the consistent estimator (np_conv4) into staging.
source ~/.bashrc >/dev/null 2>&1; conda activate /groups/clairemcwhite/envs/core_pkgs4
cd /groups/clairemcwhite/cav_workspace/bioCAV-paper/single_cell
source ../config/paths.sh; export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
ST=../figures/figure_data_staging/singlecell_C1e-4
for arm in stored np_conv4; do
  python scripts/export_fig5_case_studies.py $ST/$arm \
    --coords convergence_check/coords_${arm}.tsv --gene-corr-dir convergence_check/gene_corr_${arm}
done
echo DONE
