#!/bin/bash
# -----------------------------------------------------------------------------
# GO/EC retrain at C=1e-4 (2026-09-29), step 1: embed every unique positive
# sequence once (538 chunks of 500 from 00_build_unique_chunks.py) and pool it,
# with the flags of train_cav_from_span.sh / retrain_go/01_retrain_go_cavs.sbatch
# unchanged.  Positives are whole-sequence, so a sequence's pooled vector does not
# depend on which term it came from.  hf_embed drops sequences > 2048 aa exactly as
# the per-term runs did.  Pooled vectors are KEPT on xdisk:
#   /xdisk/clairemcwhite/clairemcwhite/go_ec_positive_embeddings/pooled/chunk_NNNN.{pos.npy,ids}
# (.ids = the U<md5> IDs of the rows, in order).  Per-residue pkls live on node-local
# scratch and are deleted per chunk.  Idempotent: resubmit to pick up failures.
# -----------------------------------------------------------------------------
#SBATCH --job-name=goec_embed
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --gres=gpu:1
#SBATCH --array=1-20
#SBATCH --output=slurm_logs/embed_%A_%a.out
#SBATCH --error=slurm_logs/embed_%A_%a.err

# NB: no `set -u` around the conda bootstrap (/etc/bashrc BASHRCSOURCED).
set -o pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../../config/paths.sh"
source ~/.bashrc
conda activate "$CONDA_ENV_DIR"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

# Fail loudly on a CPU fallback (hf_embed_new.py silently moves to CPU, ~50x slower).
python - <<'PY' || { echo "FATAL: no usable GPU on $(hostname)" >&2; exit 1; }
import sys, torch
sys.exit(0 if torch.cuda.is_available() and torch.cuda.device_count() > 0 else 1)
PY
echo "GPU OK on $(hostname): $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"

BASE=/xdisk/clairemcwhite/clairemcwhite/go_ec_positive_embeddings
OUT=$BASE/pooled; mkdir -p "$OUT"
NTASKS=${SLURM_ARRAY_TASK_COUNT:-1}; TASKID=${SLURM_ARRAY_TASK_ID:-1}
WORK="${TMPDIR:-/tmp}/goec_embed_${SLURM_ARRAY_JOB_ID:-local}_${TASKID}"
mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT
mapfile -t CHUNKS < <(ls "$BASE"/chunks/chunk_*.fasta | sort)
ok=0; skipped=0; failed=0
for i in "${!CHUNKS[@]}"; do
    (( i % NTASKS != TASKID - 1 )) && continue
    FASTA="${CHUNKS[$i]}"; N=$(basename "$FASTA" .fasta); SPAN="${FASTA%.fasta}.span"
    POS="$OUT/$N.pos.npy"; IDS="$OUT/$N.ids"
    [ -s "$POS" ] && [ -s "$IDS" ] && { skipped=$((skipped+1)); continue; }
    PKL="$WORK/$N.pkl"
    python "$HF_EMBED_SCRIPT" -f "$FASTA" -o "$PKL" --get_aa_embeddings --get_sequence_embedding \
        --strat mean -l -11 -m "$ESM_MODEL_DIR" -b 1 --max_length 2048 > "$WORK/$N.embed.log" 2>&1 \
        || { echo "FAIL $N: embed" >&2; tail -3 "$WORK/$N.embed.log" >&2; failed=$((failed+1)); rm -f "$PKL"*; continue; }
    python "$BIOCAV_REPO/scripts/prepare_embeddings.py" --pkl "$PKL" --info "${PKL}.seqnames" \
        --spans "$SPAN" --out "$POS.tmp.npy" > "$WORK/$N.pool.log" 2>&1 \
        || { echo "FAIL $N: pool" >&2; failed=$((failed+1)); rm -f "$PKL"* "$POS.tmp.npy"; continue; }
    # row IDs: span order, restricted to sequences the embedder kept (it drops > 2048 aa)
    python - "$SPAN" "${PKL}.seqnames" "$POS.tmp.npy" "$IDS.tmp" <<'PY' || { echo "FAIL $N: ids" >&2; failed=$((failed+1)); rm -f "$PKL"* "$POS.tmp.npy" "$IDS.tmp"; continue; }
import sys, numpy as np
span, names, pos, ids = sys.argv[1:]
kept = {l.split()[0].lstrip(">") for l in open(names) if l.strip()}
rows = [l.split()[0] for l in open(span) if l.strip() and l.split()[0] in kept]
n = np.load(pos).shape[0]
assert n == len(rows), (n, len(rows))
open(ids, "w").write("\n".join(rows) + "\n")
PY
    mv "$POS.tmp.npy" "$POS"; mv "$IDS.tmp" "$IDS"
    rm -f "$PKL"*
    ok=$((ok+1)); echo "OK $N ($(wc -l < "$IDS") rows)"
done
echo "Task $TASKID done: ok=$ok skipped=$skipped failed=$failed"
[ "$failed" -eq 0 ]
