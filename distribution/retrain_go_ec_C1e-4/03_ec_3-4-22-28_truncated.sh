#!/bin/bash
# -----------------------------------------------------------------------------
# GO/EC retrain at C=1e-4, completeness fix (2026-10-01): EC 3.4.22.28 (picornain 3C).
# All 76 training positives are picornavirus polyproteins of 2,052-2,432 aa, so the
# main run (which drops sequences > 2048 aa, as the shipped library did) had zero
# positives and no CAV; the shipped library has none either.  Here every positive is
# truncated to its first 2046 aa (2046 + BOS/EOS = the 2048-token model limit, so the
# tokenizer cuts nothing further) with hf_embed's own -t option; all other embed,
# pooling and training settings are those of steps 1-2.  The 3C protease domain lies
# around residues 1,500-1,750 of these polyproteins, so it is retained.
# Output: the usual layout under go_ec_cavs_C1e-4_sharedscaler/, report_v1.json flags
# "truncated_to": 2046.  Run directly on a GPU node (76 sequences).
# -----------------------------------------------------------------------------
set -o pipefail   # no -e until after the conda bootstrap (~/.bashrc unalias fails under -e)
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/../../config/paths.sh"
source ~/.bashrc
conda activate "$CONDA_ENV_DIR"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
set -e
python -c "import sys, torch; sys.exit(0 if torch.cuda.is_available() else 1)" \
    || { echo "FATAL: no usable GPU" >&2; exit 1; }

TERM_DIR=/xdisk/clairemcwhite/proFAB/EC_sets/ec_dataset_part2/ecNo_3-4-22-28
SPAN=$TERM_DIR/random_positive_train_max1000.span
WORK=/xdisk/clairemcwhite/clairemcwhite/go_ec_positive_embeddings/ec_3-4-22-28_trunc2046
mkdir -p "$WORK"

# 1. truncated FASTA (headers = span accessions, span order) + one-column span file
python - "$SPAN" "$WORK" "$HERE" <<'PY'
import importlib.util, sys
span, work, here = sys.argv[1:]
spec = importlib.util.spec_from_file_location("t", f"{here}/02_train_go_ec_C1e-4.py")
t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)
recs = t.fasta_seqs(span + ".fasta")
ids = [l.split()[0] for l in open(span) if l.strip()]
ids = [a for a in ids if a in recs]
with open(f"{work}/pos.fasta", "w") as f:
    for a in ids:
        f.write(f">{a}\n{recs[a][:2046]}\n")
open(f"{work}/pos.span", "w").write("\n".join(ids) + "\n")
print(f"{len(ids)} positives; lengths {min(len(recs[a]) for a in ids)}-{max(len(recs[a]) for a in ids)} -> 2046")
PY

# 2. embed (step-1 flags, but -t 2046 in place of --max_length 2048) and pool
python "$HF_EMBED_SCRIPT" -f "$WORK/pos.fasta" -o "$WORK/pos.pkl" --get_aa_embeddings \
    --get_sequence_embedding --strat mean -l -11 -m "$ESM_MODEL_DIR" -b 1 -t 2046 \
    > "$WORK/embed.log" 2>&1
python "$BIOCAV_REPO/scripts/prepare_embeddings.py" --pkl "$WORK/pos.pkl" \
    --info "$WORK/pos.pkl.seqnames" --spans "$WORK/pos.span" --out "$WORK/pos.npy" \
    > "$WORK/pool.log" 2>&1
rm -f "$WORK/pos.pkl"   # per-residue embeddings; the pooled pos.npy is kept

# 3. train with the step-2 settings
python "$HERE/03_train_truncated.py" "$WORK" "$SPAN"
