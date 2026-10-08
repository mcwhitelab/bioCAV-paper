#!/bin/bash
# -----------------------------------------------------------------------------
# DeepGo-SE predictions for the proFAB NEGATIVE test proteins (Figure 2
# comparator). Run alongside deepgose_val_proteins.sbatch: the AUC/AUPR in
# Figure 2 D/E need the tool scored on the same negatives as the CAVs, and
# scoring positives only is what invalidated the dropped DeepGOWeb run.
#
# PROVENANCE: the shipped predictions were produced on 2026-06-20 from
# /xdisk/clairemcwhite/test_deepgose.sbatch. This file is the same job with
# paths taken from config/paths.sh instead of hard-coded; it has not been
# re-run. See deepgose_val_proteins.sbatch for the full note.
#
# Only --device and --batch-size are set; the prediction threshold is the
# predict.py default of 0.1.
#
# One pass over every per-concept negative FASTA in the proFAB GO sets
# (~6,600 concepts), writing <fasta>_preds_{bp,cc,mf}.tsv.gz beside each input.
# -----------------------------------------------------------------------------
#SBATCH --job-name=deepgose_neg
#SBATCH --output=slurm_logs/deepgose_neg.out
#SBATCH --error=slurm_logs/deepgose_neg.err
#SBATCH --account=mcwhite
#SBATCH --partition=gpu_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL

set -o pipefail
source ~/.bashrc

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source "$SCRIPT_DIR/../../config/paths.sh"

DEEPGOSE_SIF="${DEEPGOSE_SIF:-/xdisk/clairemcwhite/deepgose_latest.sif}"
DEEPGOSE_DATA="${DEEPGOSE_DATA:-/xdisk/clairemcwhite/data}"
TORCH_CACHE="${TORCH_CACHE:-/xdisk/clairemcwhite/torch_cache}"

for r in "$DEEPGOSE_SIF" "$DEEPGOSE_DATA"; do
    [ -e "$r" ] || { echo "FATAL: missing $r" >&2; exit 1; }
done

FASTA_FILES=("$PROFAB_GO_SETS"/go_dataset_part*/*/random_negative_test_max1000.span.fasta)
TOTAL=${#FASTA_FILES[@]}
PROCESSED=0
SKIPPED=0

echo "Found $TOTAL FASTA file(s) to process"

for FASTA in "${FASTA_FILES[@]}"; do
    [ -f "$FASTA" ] || continue

    STEM="${FASTA%.*}"
    if [[ -f "${STEM}_preds_bp.tsv.gz" && -f "${STEM}_preds_cc.tsv.gz" && -f "${STEM}_preds_mf.tsv.gz" ]]; then
        ((SKIPPED++))
        continue
    fi

    echo "[PROCESSING] $FASTA"
    APPTAINERENV_TORCH_HOME="$TORCH_CACHE" \
    apptainer exec --nv \
        -B "$DEEPGOSE_DATA":/workspace/deepgo2/data \
        -B /xdisk/clairemcwhite:/xdisk/clairemcwhite \
        --pwd /workspace/deepgo2 \
        "$DEEPGOSE_SIF" \
        python predict.py \
            -if "$FASTA" \
            --device cuda \
            --batch-size 1 || { echo "[FAIL] $FASTA" >&2; exit 1; }

    ((PROCESSED++))
done

echo "Total: $TOTAL | processed: $PROCESSED | skipped: $SKIPPED"
