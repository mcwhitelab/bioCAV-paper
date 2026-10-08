#!/usr/bin/env python
"""
Fig 3A: rank the true EC against the trained EC CAVs of the same EC level.

The original ranks (compare_ec_tools.py::compute_ec_ranks) pool all 952 EC CAVs
across levels 1-4, so a level-3/4 true EC competes with its own, less specific
ancestors (e.g. 1.1 and 1 for 1.1.1). This mirrors the GO fix
(go_eval/rank_all_terms/compute_ranks_all_terms.py, *_rank_ont_filt):

  rank_level_filt — candidates are the EC CAVs at the true EC's own level
                    (ancestors are therefore excluded by construction), minus
                    every other EC the protein carries in test_ec_long.tsv,
                    truncated to that level ("filtered rank").
  rank_all        — all 952 EC CAVs, recomputed (sanity check vs. shipped cav_rank)

CAV score = shared-scaler-standardised embedding · unit concept vector.
Pairs, llr and clean_supcon_correct are taken unchanged from the shipped
figure_data/ec_specificity_ranks.csv.
"""
import glob
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

PAPER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PAPER.parent / "bioCAV" / "core"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402

EC_SETS = Path("/xdisk/clairemcwhite/proFAB/EC_sets")
EC_EVAL = PAPER / "ec_eval"
FIG_DATA = PAPER / "figures" / "figure_data"
SCALER = PAPER / "shared_protein" / "reference_population" / "scaler_v1.pkl"
GOLD = EC_EVAL / "static" / "ec_annotations" / "EC-Bench" / "test_ec_long.tsv"
VAL_PKL = EC_EVAL / "static" / "ec_annotations" / "EC-Bench" / "test_ec.fasta.pkl"
CAV_SUBDIR = "random_positive_train_max1000_cav"
SKIP = {"sp", "tr", "sw", "ref"}


def normalise_pid(pid):
    parts = [p for p in str(pid).split()[0].split("|") if p and p not in SKIP]
    return parts[0] if parts else str(pid)


def strip_wildcards(ec):
    return ".".join(p for p in str(ec).split(".") if p not in ("-", ""))


def main():
    ids, vecs = [], []
    for d in sorted(glob.glob(str(EC_SETS / "ec_dataset_part*" / "ecNo_*"))):
        cav = Path(d) / CAV_SUBDIR / "concept_v1.npy"
        if cav.exists():
            ids.append(Path(d).name[len("ecNo_"):].replace("-", "."))
            vecs.append(np.load(cav))
    cav_matrix = np.stack(vecs)
    level = np.array([len(e.split(".")) for e in ids])
    col = {e: i for i, e in enumerate(ids)}
    print(f"Loaded {len(ids)} EC CAVs; per level: "
          f"{dict(zip(*np.unique(level, return_counts=True)))}")

    gold = pd.read_csv(GOLD, sep="\t")
    true_ecs = gold.assign(ec=gold["ec_number"].map(strip_wildcards)) \
                   .groupby("protein_id")["ec"].apply(set).to_dict()

    embs, sids = load_sequence_embeddings(VAL_PKL)
    idx = {}
    for i, s in enumerate(sids):
        idx.setdefault(normalise_pid(s), i)
    S = joblib.load(SCALER).transform(embs) @ cav_matrix.T

    old = pd.read_csv(FIG_DATA / "ec_specificity_ranks.csv", dtype={"ec_number": str})
    rows = []
    for r in old.itertuples(index=False):
        i, j = idx[r.protein_id], col[r.ec_number]
        s = S[i]
        lvl = len(r.ec_number.split("."))
        others = {".".join(e.split(".")[:lvl]) for e in true_ecs.get(r.protein_id, ())
                  if len(e.split(".")) >= lvl} - {r.ec_number}
        keep = level == lvl
        keep[[col[e] for e in others if e in col]] = False
        rows.append({
            "protein_id": r.protein_id, "ec_number": r.ec_number, "level": r.level,
            "llr": r.llr, "clean_supcon_correct": r.clean_supcon_correct,
            "cav_rank_level_filt": int((s[keep] > s[j]).sum()) + 1,
            "cav_rank_all": int((s > s[j]).sum()) + 1,
            "cav_rank_all_shipped": r.cav_rank,
            "n_cavs_level": int((level == lvl).sum()),
            "n_dropped_filt": len([e for e in others if e in col]),
        })
    out = pd.DataFrame(rows)
    dest = FIG_DATA / "ec_specificity_ranks_samelevel.csv"
    out.to_csv(dest, index=False)

    print(f"{len(out)} pairs -> {dest.name}")
    print(f"cav_rank_all == shipped: {(out.cav_rank_all == out.cav_rank_all_shipped).mean():.3f}")
    for lvl, g in out.groupby("level"):
        for k in ["all_shipped", "level_filt"]:
            c = g[f"cav_rank_{k}"]
            print(f"  level {lvl} (n={len(g)}, {g.n_cavs_level.iloc[0]} same-level CAVs) {k:12s} "
                  f"rank1 {100*(c==1).mean():5.1f}%  top3 {100*(c<=3).mean():5.1f}%  "
                  f"top10 {100*(c<=10).mean():5.1f}%  median {c.median():.0f}")
    print("pairs with other same-level true ECs dropped:", int((out.n_dropped_filt > 0).sum()))


if __name__ == "__main__":
    main()
