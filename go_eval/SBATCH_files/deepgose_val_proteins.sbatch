#!/bin/bash
# -----------------------------------------------------------------------------
# DeepGo-SE predictions for the GO temporal-holdout VALIDATION proteins
# (Figure 2 comparator).
#
# PROVENANCE: the run that produced the predictions shipped in
# go_eval/outside_tools/ was executed on 2026-06-20 from
# /xdisk/clairemcwhite/test_deepgose_val.sbatch, which hard-coded the paths of
# a previous repository layout (claire_workspace/cav/protein_analysis/). That
# tree still exists and its outputs are md5-identical to the copies in this
# repository. This file is the same job repointed at the current layout via
# config/paths.sh; it has not been re-run, because doing so would only
# reproduce byte-identical output.
#
# Only --device and --batch-size are set; every other DeepGo-SE parameter is
# the predict.py default, INCLUDING the prediction threshold of 0.1
# (--threshold/-t, default=0.1). Do not confuse this with the DeepGOWeb API
# run, which used threshold 0 and was dropped from the paper.
#
# Outputs land beside each input FASTA as <fasta>_preds_{bp,cc,mf}.tsv.gz.
# go_eval/02_format_outside_tool_preds.sh gunzips and concatenates them.
# -----------------------------------------------------------------------------
#SBATCH --job-name=deepgose_val
#SBATCH --output=slurm_logs/deepgose_val.out
#SBATCH --error=slurm_logs/deepgose_val.err
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

# Container and model data live outside the repository; override if they move.
DEEPGOSE_SIF="${DEEPGOSE_SIF:-/xdisk/clairemcwhite/deepgose_latest.sif}"
DEEPGOSE_DATA="${DEEPGOSE_DATA:-/xdisk/clairemcwhite/data}"
TORCH_CACHE="${TORCH_CACHE:-/xdisk/clairemcwhite/torch_cache}"

for r in "$DEEPGOSE_SIF" "$DEEPGOSE_DATA"; do
    [ -e "$r" ] || { echo "FATAL: missing $r" >&2; exit 1; }
done

FASTA_FILES=("$BIOCAV_PAPER_ROOT"/go_eval/results/temporal_eval_*/eval_embeddings/val_proteins.span.fasta)
TOTAL=${#FASTA_FILES[@]}
PROCESSED=0
SKIPPED=0

echo "Found $TOTAL FASTA file(s) to process"

for FASTA in "${FASTA_FILES[@]}"; do
    [ -f "$FASTA" ] || { echo "[WARN] not found: $FASTA" >&2; continue; }

    # predict.py writes <fasta-stem>_preds_<ont>.tsv.gz next to the input, so
    # the skip test must look for exactly those names. (The original script
    # checked a different, never-written filename, so it never skipped.)
    STEM="${FASTA%.*}"
    if [[ -f "${STEM}_preds_bp.tsv.gz" && -f "${STEM}_preds_cc.tsv.gz" && -f "${STEM}_preds_mf.tsv.gz" ]]; then
        ((SKIPPED++))
        echo "[SKIP] already processed: $FASTA"
        continue
    fi

    echo "[PROCESSING] $FASTA"
    APPTAINERENV_TORCH_HOME="$TORCH_CACHE" \
    apptainer exec --nv \
        -B "$DEEPGOSE_DATA":/workspace/deepgo2/data \
        -B /xdisk/clairemcwhite:/xdisk/clairemcwhite \
        -B "$BIOCAV_PAPER_ROOT":"$BIOCAV_PAPER_ROOT" \
        --pwd /workspace/deepgo2 \
        "$DEEPGOSE_SIF" \
        python predict.py \
            -if "$FASTA" \
            --device cuda \
            --batch-size 1 || { echo "[FAIL] $FASTA" >&2; exit 1; }

    ((PROCESSED++))
    echo "[COMPLETED] $FASTA"
done

echo "Total: $TOTAL | processed: $PROCESSED | skipped: $SKIPPED"
