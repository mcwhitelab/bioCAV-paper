#! /bin/bash

# -------------
### Directives
# -------------
#SBATCH --job-name=chlamy_mixed_cav
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
# Last-ditch narrower hypothesis: instead of scoring peps12 against CAVs
# trained on whole RpoZ sequences, restrict the CAV to just the RpoZ
# binding-interface coil (Q31MH9 / RPOZ_SYNE7 residues 8-46, 1-indexed,
# 39 aa) mapped onto each homolog via pairwise alignment
# (make_interface_spans.py). peps12 is hypothesized to bind the same
# pocket as RpoZ despite diverged overall structure, so a narrower,
# interface-only concept may pick up convergent signal that whole-sequence
# or whole-domain-length windowing missed.
#
# WINDOW = 39 : length of the Q31MH9 interface span itself (residues 8-46)

RUNDIR=fastas/pairs/rpoz_chlamy
FASTA=$RUNDIR/chlamy_proteome.fasta
WINDOW=39
OUT=$RUNDIR/chlamy_window_scores_mixed.pkl

# -------------
### Paths
# -------------
source ~/.bashrc

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../config/paths.sh"
conda activate $CONDA_ENV_DIR
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

python "$RUNDIR/window_score_proteome_scaled.py" \
    --fasta "$FASTA" \
    --model "$ESM_MODEL_DIR" \
    --layer -11 \
    --window "$WINDOW" \
    --chunk-size 50 \
    --cav mixed="$RUNDIR/mixed_interface_concept" \
    --out "$OUT"
