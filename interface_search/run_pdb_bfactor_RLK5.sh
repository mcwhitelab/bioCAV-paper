#! /bin/bash

# -------------
### Inputs — download AlphaFold models before running, save under fastas/pairs/lrr/
# -------------
#   RLK5_ARATH (P47735): https://alphafold.ebi.ac.uk/entry/P47735
#   XA21_ORYSJ  (Q2R2D5): https://alphafold.ebi.ac.uk/entry/Q2R2D5
#   BRI1_ARATH  (O22476): https://alphafold.ebi.ac.uk/entry/O22476

SPAN_FILE=fastas/pairs/lrr/RLK5_ARATH.fasta.span
CONCEPT_DIR=fastas/pairs/lrr/RLK5_ARATH.fasta_concept

# -------------
### Paths
# -------------
source ~/.bashrc
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"
conda activate $CONDA_ENV_DIR

paint() {
    local pdb_file=$1
    local out_pdb=$2

    if [ ! -f "$pdb_file" ]; then
        echo "PDB not found: $pdb_file — skipping"
        return
    fi

    python $BIOCAV_REPO/specific_scripts/pdb_bfactor_from_cav.py \
        --pdb          $pdb_file \
        --cav-dir      $CONCEPT_DIR \
        --spans        $SPAN_FILE \
        --embed-script $HF_EMBED_SCRIPT \
        --model        $ESM_MODEL_DIR \
        --out          $out_pdb
}

paint fastas/pairs/lrr/AF-P47735-F1-model_v6.pdb fastas/pairs/lrr/RLK5_ARATH_cav_bfactor.pdb
paint fastas/pairs/lrr/AF-Q2R2D5-F1-model_v6.pdb fastas/pairs/lrr/XA21_ORYSJ_cav_bfactor.pdb
paint fastas/pairs/lrr/AF-O22476-F1-model_v6.pdb fastas/pairs/lrr/BRI1_ARATH_cav_bfactor.pdb
paint fastas/pairs/lrr/AF-O48849-F1-model_v6.pdb fastas/pairs/lrr/RLP23_ARATH_cav_bfactor.pdb
paint fastas/pairs/lrr/AF-Q9FII5-F1-model_v6.pdb fastas/pairs/lrr/TDR_ARATH_cav_bfactor.pdb

