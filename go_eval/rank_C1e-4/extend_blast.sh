#!/bin/bash
# -----------------------------------------------------------------------------
# Extend the DIAMOND nearest-neighbour baselines (Fig 2C GO, Fig 3A EC) to the two
# new-only C=1e-4 CAVs (GO:0032436, EC 3.4.22.28; 2026-10-01), so BLAST is ranked over
# the same candidate terms as CAV and DeepGo-SE. The original searches are untouched.
#
# The original searches ran with --max-target-seqs 0, i.e. every hit at e <= 1e-3 against
# the union of all CAV training positives. So for a new term's positives that are ALREADY
# in that database, the original hits are already complete, and adding the term only needs
# its (term, accession) member rows. That is exactly what having the term in the original
# search would give. This holds for all positives of both terms (GO:0032436 280/280, EC
# 3.4.22.28 76/76; they are positives of other terms too). Only positives MISSING from
# the original database would be searched (none here), with the original settings
# (DIAMOND 2.0.9, blastp --more-sensitive --max-target-seqs 0 --evalue 1e-3,
# outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore). They use
# --dbsize = the original database's letters, so e-values (and the 1e-3 cutoff) stay on
# the original scale. Missing accessions belong to no other term, so their hits cannot
# leak into other terms.
# (Revised 2026-10-01. A first version searched all new positives against a small
# database with tagged subject IDs, which only approximated this.)
# The EC 3.4.22.28 positives are searched at full length (truncation applies only to
# embedding). Outputs: <blast_baseline>/C1e-4_extra/{hits, members}.
# -----------------------------------------------------------------------------
set -euo pipefail
PAPER=/groups/clairemcwhite/cav_workspace/bioCAV-paper
PY=/groups/clairemcwhite/envs/core_pkgs4/bin/python
module load diamond/2.0.9 2>/dev/null || { source /etc/profile.d/lmod.sh 2>/dev/null; module load diamond/2.0.9; }
DMD_ARGS=(--more-sensitive --max-target-seqs 0 --evalue 1e-3 --threads 8
          --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore)

run() {  # term_label span queries orig_hits orig_members outdir hits_name members_name orig_db orig_fasta
    local label=$1 span=$2 q=$3 hits0=$4 mem0=$5 out=$6 hits_name=$7 mem_name=$8 db0=$9 fa0=${10}
    local letters; letters=$(diamond dbinfo -d "$db0" | awk '/Letters/ {print $2}')
    rm -rf "$out"; mkdir -p "$out"
    $PY - "$span" "$label" "$out" "$fa0" <<'PY'
import sys
span, label, out, fa0 = sys.argv[1:]
def acc(h):
    p = h[1:].split()[0].split("|")
    return p[1] if len(p) >= 2 and p[0] in ("sp", "tr") else p[0]
recs, h = {}, None
for l in open(span + ".fasta"):
    l = l.strip()
    if l.startswith(">"): h = acc(l); recs[h] = []
    elif h: recs[h].append(l)
ids = [l.split()[0] for l in open(span) if l.strip()]
ids = [a for a in ids if a in recs]
in_db = {l[1:].split()[0] for l in open(fa0) if l.startswith(">")}
missing = [a for a in ids if a not in in_db]
with open(f"{out}/extra_members.tsv", "w") as f:
    for a in ids: f.write(f"{label}\t{a}\n")
with open(f"{out}/missing_positives.fa", "w") as f:
    for a in missing: f.write(f">{a}\n{''.join(recs[a])}\n")
print(f"{label}: {len(ids)} training positives, {len(ids) - len(missing)} already in the original "
      f"database, {len(missing)} to search")
PY
    : > "$out/extra_hits.tsv"
    if [ -s "$out/missing_positives.fa" ]; then
        diamond makedb --in "$out/missing_positives.fa" -d "$out/missing_positives" --quiet
        diamond blastp -q "$q" -d "$out/missing_positives.dmnd" -o "$out/extra_hits.raw.tsv" \
            "${DMD_ARGS[@]}" --dbsize "$letters" --quiet
        awk -F'\t' '$1 != $2' "$out/extra_hits.raw.tsv" > "$out/extra_hits.tsv"
    fi
    cat "$hits0" "$out/extra_hits.tsv" > "$out/$hits_name"
    cat "$mem0"  "$out/extra_members.tsv" > "$out/$mem_name"
    echo "$label: $(wc -l < "$out/extra_hits.tsv") new hits (dbsize $letters) -> $out/$hits_name"
}

GB=$PAPER/go_eval/blast_baseline
run GO:0032436 /xdisk/clairemcwhite/proFAB/GO_sets/go_dataset_part5/GO_0032436/random_positive_train_max1000.span \
    $GB/val_queries.fa $GB/val_vs_train.tsv $GB/term_members.tsv $GB/C1e-4_extra val_vs_train.tsv term_members.tsv \
    $GB/train_positives_all.dmnd $GB/train_positives_all.fa
EB=$PAPER/ec_eval/blast_baseline
run 3.4.22.28 /xdisk/clairemcwhite/proFAB/EC_sets/ec_dataset_part2/ecNo_3-4-22-28/random_positive_train_max1000.span \
    $PAPER/ec_eval/static/ec_annotations/EC-Bench/test_ec.fasta $EB/test_vs_ectrain.tsv $EB/ec_term_members.tsv \
    $EB/C1e-4_extra test_vs_ectrain.tsv ec_term_members.tsv $EB/ec_train_positives_all.dmnd $EB/ec_train_positives_all.fa
