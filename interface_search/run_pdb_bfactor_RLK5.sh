#! /bin/bash

# -------------
### Inputs — replace PDB_FILE before running
# -------------
# PDB_FILE: AlphaFold (or other) structure for RLK5_ARATH (UniProt P47735).
#           Download e.g. from https://alphafold.ebi.ac.uk/entry/P47735
#           and save it under fastas/pairs/lrr/.

PDB_FILE=fastas/pairs/lrr/AF-P47735-F1-model_v6.pdb
INPUT_FASTA=fastas/pairs/lrr/RLK5_ARATH.fasta
SPAN_FILE=fastas/pairs/lrr/RLK5_ARATH.fasta.span
CONCEPT_DIR=fastas/pairs/lrr/RLK5_ARATH.fasta_concept
OUT_PDB=fastas/pairs/lrr/RLK5_ARATH_cav_bfactor.pdb

# -------------
### Paths
# -------------
source ~/.bashrc
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"
conda activate $CONDA_ENV_DIR

if [ ! -f "$PDB_FILE" ]; then
    echo "PDB_FILE not found: $PDB_FILE"
    echo "Download the AlphaFold model for P47735 (RLK5_ARATH) and place it there."
    exit 1
fi

python $BIOCAV_REPO/specific_scripts/pdb_bfactor_from_cav.py \
    --pdb       $PDB_FILE \
    --pkl       ${INPUT_FASTA}.pkl \
    --cav-dir   $CONCEPT_DIR \
    --spans     $SPAN_FILE \
    --out       $OUT_PDB


# Now paint on XA21

PDB_FILE=fastas/pairs/lrr/AF-Q2R2D5-F1-model_v6.pdb
OUT_PDB=fastas/pairs/lrr/XA21_ORYSJ_cav_bfactor.pdb

python $BIOCAV_REPO/specific_scripts/pdb_bfactor_from_cav.py \
    --pdb       $PDB_FILE \
    --pkl       ${INPUT_FASTA}.pkl \
    --cav-dir   $CONCEPT_DIR \
    --spans     $SPAN_FILE \
    --out       $OUT_PDB


PDB_FILE=fastas/pairs/lrr/AF-O22476-F1-model_v6.pdb
OUT_PDB=fastas/pairs/lrr/BRI1_ARATH_cav_bfactor.pdb

python $BIOCAV_REPO/specific_scripts/pdb_bfactor_from_cav.py \
    --pdb       $PDB_FILE \
    --pkl       ${INPUT_FASTA}.pkl \
    --cav-dir   $CONCEPT_DIR \
    --spans     $SPAN_FILE \
    --out       $OUT_PDB



