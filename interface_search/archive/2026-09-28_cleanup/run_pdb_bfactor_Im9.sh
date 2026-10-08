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
# Paint the Im9-trained (E9-binding) CAV onto both the cognate complex
# (1EMV, E9 DNase-Im9) and the non-cognate complex (2WPT, E9 DNase-Im2).
# Chain A is the Im protein in both structures (86 aa); chain B is E9
# DNase — default --chain (first chain) picks the Im protein automatically.
# --window-size 1 gives unsmoothed per-residue scores, since we want to see
# whether signal concentrates specifically on helix II (~res 30-44,
# specificity residues ~33/34/38/41) rather than spread across the fold.

CONCEPT_DIR=fastas/pairs/colicin/IMM9_ECOLX.fasta_concept

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
        --window-size  1 \
        --embed-script $HF_EMBED_SCRIPT \
        --model        $ESM_MODEL_DIR \
        --out          $out_pdb
}

# Cognate: E9 DNase-Im9 (1.7 A)
paint fastas/pairs/colicin/1EMV.pdb fastas/pairs/colicin/Im9_cognate_1EMV_cav_bfactor.pdb

# Non-cognate: E9 DNase-Im2 (1.77 A) — same overall binding arrangement,
# ~7 orders of magnitude weaker affinity for E9
paint fastas/pairs/colicin/2WPT.pdb fastas/pairs/colicin/Im2_noncognate_2WPT_cav_bfactor.pdb
