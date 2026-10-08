#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=pdb_bfactor
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL

# -------------
### Inputs
# -------------
# Same AlphaFold models as run_pdb_bfactor_RLK5.sh, repainted with the
# multi-example SERK-interface CAV (RLK5 + PSKR1 + FLS2 positives) instead
# of the RLK5-only CAV.
#   RLK5_ARATH (P47735): https://alphafold.ebi.ac.uk/entry/P47735
#   XA21_ORYSJ  (Q2R2D5): https://alphafold.ebi.ac.uk/entry/Q2R2D5
#   BRI1_ARATH  (O22476): https://alphafold.ebi.ac.uk/entry/O22476

SPAN_FILE=fastas/pairs/lrr/SERK_interface_multi.fasta.span
CONCEPT_DIR=fastas/pairs/lrr/SERK_interface_multi.fasta_concept

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

paint fastas/pairs/lrr/AF-P47735-F1-model_v6.pdb fastas/pairs/lrr/RLK5_ARATH_serkmulti_cav_bfactor.pdb
paint fastas/pairs/lrr/AF-Q2R2D5-F1-model_v6.pdb fastas/pairs/lrr/XA21_ORYSJ_serkmulti_cav_bfactor.pdb
paint fastas/pairs/lrr/AF-O22476-F1-model_v6.pdb fastas/pairs/lrr/BRI1_ARATH_serkmulti_cav_bfactor.pdb
paint fastas/pairs/lrr/AF-O48849-F1-model_v6.pdb fastas/pairs/lrr/RLP23_ARATH_serkmulti_cav_bfactor.pdb
paint fastas/pairs/lrr/AF-Q9FII5-F1-model_v6.pdb fastas/pairs/lrr/TDR_ARATH_serkmulti_cav_bfactor.pdb
