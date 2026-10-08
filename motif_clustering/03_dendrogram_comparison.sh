#!/usr/bin/env bash
# Head-to-head of four dendrogram constructions on the 27k CAV library.
# All four are cut with fcluster at matched cluster counts, so the ARI /
# per-clan-F1 curves are directly comparable across methods.
set -u

source "$(dirname "$0")/../config/paths.sh"
CAVS="/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs/*/L25_concept_v1.npy"
OUT=results/dendrogram_comparison
mkdir -p "$OUT"

run () {                       # run <label> <extra cav_cluster args…>
    label=$1; shift
    echo "########## $label ##########"
    python "$(dirname "$0")/cav_cluster.py" \
        --cav-pattern      "$CAVS" \
        --pfam-annotations results/figures/pfamA.txt \
        --clan-annotations results/figures/Pfam-A.clans.tsv \
        --cut-mode maxclust --no-shuffle-ids \
        --out "$OUT/$label.csv" "$@" || return 1

    python "$(dirname "$0")/score_dendrogram.py" \
        --clusters     "$OUT/$label.csv" \
        --cav-pattern  "$CAVS" \
        --label        "$label" \
        --out          "$OUT/${label}.scores.csv" \
        --per-clan-out "$OUT/${label}.per_clan.csv"
}

run cosine_average    --distance cosine    --linkage average   # historical baseline
run cosine_ward       --distance cosine    --linkage ward      # linkage-only fix
run diffusion_ward    --distance diffusion --linkage ward      # walktrap reweighting
run snn_average       --distance snn       --linkage average   # cheap sanity check

echo "########## summary ##########"
python - "$OUT" <<'PY'
import sys, glob, pandas as pd
files = sorted(glob.glob(f"{sys.argv[1]}/*.scores.csv"))
if not files:
    sys.exit("no score files")
df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
print("\n-- cut quality vs clans --")
print(df[["run","column","n_clusters_all","ARI","AMI","median_per_clan_F1"]]
      .to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print("\n-- best per run --")
best = df.loc[df.groupby("run")["ARI"].idxmax(),
              ["run","column","ARI","median_per_clan_F1"]]
print(best.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print("\n-- leaf-order contiguity & cophenetic (one row per run) --")
cols = ["run","median_norm_span","perm_median_norm_span","contiguity_ratio",
        "frac_clans_norm_span_lt2","cophenetic_vs_rawcosine_spearman"]
cols = [c for c in cols if c in df.columns]
print(df.groupby("run").first().reset_index()[cols]
      .to_string(index=False, float_format=lambda x: f"{x:.4f}"))
df.to_csv(f"{sys.argv[1]}/all_scores.csv", index=False)
print(f"\nwrote {sys.argv[1]}/all_scores.csv")
PY
