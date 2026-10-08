#!/bin/bash
# -----------------------------------------------------------------------------
# Retrain the GO CAVs that were fit in their own scaler frame so that they use
# the shared reference scaler instead.
#
# Why: 3,930 of 6,677 GO CAVs were built by 02_make_cav_commands.sh, which never
# passed --scaler-pkl, so train_cav_from_embeddings.py fit a fresh StandardScaler
# on (positives + 10k negatives).  But eval_go_temporal.py and
# 03_compare_to_outside_tools.sh standardise with the shared scaler_v1.pkl, so
# those CAVs are scored in a frame they were not trained in.  The remaining 2,747
# already use the shared scaler and are left untouched.
#
# The pooled positive embeddings were deleted by train_cav_from_span.sh step 5,
# so each concept must be re-embedded.  Embedding flags here are byte-identical
# to the original run — only --scaler-pkl is added.
#
# Run 00_backup_originals.sh first; this overwrites the CAV directories in place.
# Idempotent: a concept already in the shared frame is skipped, so the array can
# be resubmitted to pick up failures.
# -----------------------------------------------------------------------------
#SBATCH --job-name=go_cav_rescale
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --gres=gpu:1
#SBATCH --array=1-40%10
#SBATCH --output=slurm_logs/retrain_%A_%a.out
#SBATCH --error=slurm_logs/retrain_%A_%a.err

# NB: no `set -u` around the conda bootstrap — /etc/bashrc references an unbound
# variable (BASHRCSOURCED) and aborts the task in one second under `set -u`.
set -o pipefail

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../../config/paths.sh"

source ~/.bashrc
conda activate "$CONDA_ENV_DIR"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

REF_DIR="$BIOCAV_PAPER_ROOT/shared_protein/reference_population"
REF_NEG="$REF_DIR/neg_10000.pkl"
REF_SCALER="$REF_DIR/scaler_v1.pkl"
SPAN_LIST="$SCRIPT_DIR/retrain_spans_npos.tsv"   # span<TAB>n_positive from the original report

for f in "$REF_NEG" "$REF_SCALER" "$SPAN_LIST"; do
    [ -f "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }
done

# Fail loudly on a CPU fallback.  hf_embed_new.py prints "Moving model to cpu"
# and carries on when CUDA is unavailable; at ~20x slower that silently blows
# the walltime instead of erroring.
python - <<'PY' || { echo "FATAL: no usable GPU on $(hostname)" >&2; exit 1; }
import sys, torch
sys.exit(0 if torch.cuda.is_available() and torch.cuda.device_count() > 0 else 1)
PY
echo "GPU OK on $(hostname): $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"

NTASKS=${SLURM_ARRAY_TASK_COUNT:-1}
TASKID=${SLURM_ARRAY_TASK_ID:-1}

# Per-concept scratch on node-local disk: aa_embeddings for one concept can reach
# ~1 GB, and we delete them immediately, so they never touch /xdisk.
WORK="${TMPDIR:-/tmp}/go_cav_rescale_${SLURM_ARRAY_JOB_ID:-local}_${TASKID}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

ok=0; skipped=0; failed=0
mapfile -t SPANS < "$SPAN_LIST"
echo "Task $TASKID/$NTASKS over ${#SPANS[@]} concepts; scratch=$WORK"

for i in "${!SPANS[@]}"; do
    (( i % NTASKS != TASKID - 1 )) && continue
    SPAN="${SPANS[$i]%%$'\t'*}"
    EXPECT_NPOS="${SPANS[$i]##*$'\t'}"
    CAV_DIR="${SPAN%.span}_cav"
    FASTA="${SPAN}.fasta"
    NAME="$(basename "$(dirname "$SPAN")")"

    if [ ! -f "$FASTA" ]; then
        echo "FAIL  $NAME: no $FASTA" >&2; failed=$((failed+1)); continue
    fi

    # Idempotence: n_samples_seen_ == 10000 means this CAV is already in the
    # shared frame (own-frame scalers saw 10000 + n_positive).
    if [ -f "$CAV_DIR/scaler_v1.pkl" ] && python - "$CAV_DIR/scaler_v1.pkl" <<'PY'
import sys, joblib
sys.exit(0 if float(joblib.load(sys.argv[1]).n_samples_seen_) == 10000 else 1)
PY
    then
        echo "SKIP  $NAME: already shared frame"; skipped=$((skipped+1)); continue
    fi

    PKL="$WORK/${NAME}.pkl"
    POS="$WORK/${NAME}.pos.npy"

    # Step 1: embed — flags identical to the original train_cav_from_span.sh run.
    python "$HF_EMBED_SCRIPT" \
        -f "$FASTA" \
        -o "$PKL" \
        --get_aa_embeddings \
        --get_sequence_embedding \
        --strat mean \
        -l -11 \
        -m "$ESM_MODEL_DIR" \
        -b 1 \
        --max_length 2048 || { echo "FAIL  $NAME: embed" >&2; failed=$((failed+1)); rm -f "$PKL"*; continue; }

    # Step 2: pool per-residue embeddings over the concept spans.
    python "$BIOCAV_REPO/scripts/prepare_embeddings.py" \
        --pkl "$PKL" --info "${PKL}.seqnames" --spans "$SPAN" --out "$POS" \
        || { echo "FAIL  $NAME: pool" >&2; failed=$((failed+1)); rm -f "$PKL"* "$POS"; continue; }

    # Step 3: retrain in the SHARED frame.  --cv-folds 5 matches the original
    # config (regularization_C 1.0, random_seed 42 are the script defaults).
    python "$BIOCAV_REPO/scripts/train_cav_from_embeddings.py" \
        --pos "$POS" \
        --neg "$REF_NEG" \
        --out "$CAV_DIR" \
        --cv-folds 5 \
        --scaler-pkl "$REF_SCALER" \
        || { echo "FAIL  $NAME: train" >&2; failed=$((failed+1)); rm -f "$PKL"* "$POS"; continue; }

    # Guard: the rescale must not change which positives survived.  A drop here
    # would mean the retrained CAV is weaker than the one it replaced, not merely
    # re-standardised.
    GOT_NPOS=$(python -c "import json,sys;print(json.load(open(sys.argv[1]))['concept_metrics']['n_positive'])" "$CAV_DIR/report_v1.json" 2>/dev/null)
    if [ "$GOT_NPOS" != "$EXPECT_NPOS" ]; then
        echo "MISMATCH $NAME: n_positive $EXPECT_NPOS -> $GOT_NPOS" >&2
        failed=$((failed+1))
    else
        ok=$((ok+1))
        echo "OK    $NAME (n_positive=$GOT_NPOS)"
    fi
    rm -f "$PKL"* "$POS" "$CAV_DIR/pos.npy" "$CAV_DIR/neg.npy"
done

echo "Task $TASKID done: ok=$ok skipped=$skipped failed=$failed"
[ "$failed" -eq 0 ]
