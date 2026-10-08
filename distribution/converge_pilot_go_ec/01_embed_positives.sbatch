#!/bin/bash
# -----------------------------------------------------------------------------
# GO/EC convergence pilot: re-embed each sampled term's proFAB random train
# positives and pool over the spans, exactly as train_cav_from_span.sh and
# distribution/retrain_go/01_retrain_go_cavs.sbatch do (same hf_embed flags,
# same prepare_embeddings.py).  Unlike those, the pooled positives are KEPT:
#   /xdisk/clairemcwhite/clairemcwhite/converge_pilot_go_ec/pos/<LIB>_<term>.pos.npy
# No CAV is trained here and no shipped CAV directory is touched.
# Idempotent: a term whose pos.npy exists is skipped.
# -----------------------------------------------------------------------------
#SBATCH --job-name=conv_pilot_embed
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=3:00:00
#SBATCH --gres=gpu:1
#SBATCH --array=1-10
#SBATCH --output=slurm_logs/embed_%A_%a.out
#SBATCH --error=slurm_logs/embed_%A_%a.err

# NB: no `set -u` around the conda bootstrap (/etc/bashrc BASHRCSOURCED).
set -o pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../../config/paths.sh"
source ~/.bashrc
conda activate "$CONDA_ENV_DIR"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

python - <<'PY' || { echo "FATAL: no usable GPU on $(hostname)" >&2; exit 1; }
import sys, torch
sys.exit(0 if torch.cuda.is_available() and torch.cuda.device_count() > 0 else 1)
PY
echo "GPU OK on $(hostname): $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"

OUT=/xdisk/clairemcwhite/clairemcwhite/converge_pilot_go_ec/pos
NTASKS=${SLURM_ARRAY_TASK_COUNT:-1}; TASKID=${SLURM_ARRAY_TASK_ID:-1}
WORK="${TMPDIR:-/tmp}/conv_pilot_${SLURM_ARRAY_JOB_ID:-local}_${TASKID}"
mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT
mapfile -t ROWS < "$SCRIPT_DIR/embed_list.tsv"
ok=0; skipped=0; failed=0
for i in "${!ROWS[@]}"; do
    (( i % NTASKS != TASKID - 1 )) && continue
    IFS=$'\t' read -r LIB TERM SPAN <<< "${ROWS[$i]}"
    NAME="${LIB}_${TERM//[:.]/_}"
    POS="$OUT/${NAME}.pos.npy"
    [ -s "$POS" ] && { skipped=$((skipped+1)); continue; }
    FASTA="${SPAN}.fasta"
    [ -f "$FASTA" ] || { echo "FAIL $NAME: no $FASTA" >&2; failed=$((failed+1)); continue; }
    PKL="$WORK/${NAME}.pkl"
    python "$HF_EMBED_SCRIPT" -f "$FASTA" -o "$PKL" --get_aa_embeddings --get_sequence_embedding \
        --strat mean -l -11 -m "$ESM_MODEL_DIR" -b 1 --max_length 2048 \
        || { echo "FAIL $NAME: embed" >&2; failed=$((failed+1)); rm -f "$PKL"*; continue; }
    python "$BIOCAV_REPO/scripts/prepare_embeddings.py" --pkl "$PKL" --info "${PKL}.seqnames" \
        --spans "$SPAN" --out "$POS" \
        || { echo "FAIL $NAME: pool" >&2; failed=$((failed+1)); rm -f "$PKL"* "$POS"; continue; }
    rm -f "$PKL"*
    ok=$((ok+1)); echo "OK $NAME"
done
echo "Task $TASKID done: ok=$ok skipped=$skipped failed=$failed"
[ "$failed" -eq 0 ]
