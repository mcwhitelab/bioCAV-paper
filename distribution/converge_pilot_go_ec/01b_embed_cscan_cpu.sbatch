#!/bin/bash
# -----------------------------------------------------------------------------
# GO/EC C-scan subset (2026-09-29): the GPU pilot (job 24046079) was cancelled at the
# user's request; this embeds the smaller subset in terms_cscan.tsv / embed_list_cscan.tsv
# (30 GO + 15 EC terms, 6,687 positives) on CPU nodes, one term per array task.
# Embedding flags are identical to 01_embed_positives.sbatch / train_cav_from_span.sh;
# hf_embed_new.py runs on CPU when no GPU is present.  Pooled positives are KEPT:
#   /xdisk/clairemcwhite/clairemcwhite/converge_pilot_go_ec/pos/<LIB>_<term>.pos.npy
# No CAV is trained here and no shipped CAV directory is touched.  Idempotent.
# -----------------------------------------------------------------------------
#SBATCH --job-name=cscan_embed_cpu
#SBATCH --account=mcwhite
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --array=1-45
#SBATCH --output=slurm_logs/cscan_embed_%A_%a.out
#SBATCH --error=slurm_logs/cscan_embed_%A_%a.err
# NB: no `set -u` around the conda bootstrap (/etc/bashrc BASHRCSOURCED).
set -o pipefail
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../../config/paths.sh"
source ~/.bashrc
conda activate "$CONDA_ENV_DIR"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "CPU embedding on $(hostname), $SLURM_CPUS_PER_TASK cores"

OUT=/xdisk/clairemcwhite/clairemcwhite/converge_pilot_go_ec/pos
NTASKS=${SLURM_ARRAY_TASK_COUNT:-1}; TASKID=${SLURM_ARRAY_TASK_ID:-1}
WORK="${TMPDIR:-/tmp}/conv_pilot_${SLURM_ARRAY_JOB_ID:-local}_${TASKID}"
mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT
mapfile -t ROWS < "$SCRIPT_DIR/embed_list_cscan.tsv"
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
