#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=chlamy_window_cav
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL

# -------------
### Inputs
# -------------
# peps12 (367 aa) scored higher against the RpoZ CAVs than expected once
# scored per-window instead of whole-sequence mean pool (see peps12
# sliding-window results). This job repeats that windowing across the whole
# Chlamydomonas proteome (18,579 seqs) so peps12's windowed rank can be
# compared against everyone else's windowed best score, not their
# whole-sequence score — an apples-to-apples background.
#
# WINDOW = 91 : mean length of the RpoZ family (rpoz_midsize.fasta /
#               rpoz_family_reviewed.fasta both average ~91 aa)

RUNDIR=fastas/pairs/rpoz_chlamy
FASTA=$RUNDIR/chlamy_proteome.fasta
WINDOW=91
OUT=$RUNDIR/chlamy_window_scores.pkl

# -------------
### Paths
# -------------
source ~/.bashrc

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"
conda activate $CONDA_ENV_DIR
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

python "$RUNDIR/window_score_proteome.py" \
    --fasta "$FASTA" \
    --model "$ESM_MODEL_DIR" \
    --layer -11 \
    --window "$WINDOW" \
    --chunk-size 50 \
    --cav full="$RUNDIR/rpoz_cav.pkl" \
    --cav midsize="$RUNDIR/rpoz_midsize_cav.pkl" \
    --cav single="$RUNDIR/rpoz_single_cav.pkl" \
    --out "$OUT"
