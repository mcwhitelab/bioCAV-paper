#!/usr/bin/env python
"""
Fig 3A with the GO-side fixes applied (exploratory, 2026-09-29).

Pairs: figure_data/ec_specificity_ranks.csv (the shipped Fig 3A pairs, EC level
3 and 4), as used by rank_same_level/compute_ec_ranks_same_level.py.

1. Known-positive exclusion. Drop a pair when the protein is already a proFAB
   positive for that EC: random train or test span of the EC's own CAV. The old
   ec_train_overlap_pairs.tsv list only covered exact-EC pairs of the gold
   standard, so it missed most level-3 pairs.
2. Strict candidate filter (_pf). Candidates are the same-level EC CAVs, minus
   the protein's other true ECs in test_ec_long.tsv (truncated to the level),
   minus every same-level EC for which proFAB lists the protein as a positive.
3. Scores and ranks, all over identical candidates:
     CAV   = shared-scaler-standardised embedding . unit concept vector (CAVs on disk)
     BLAST = best DIAMOND bitscore to the EC CAV's training positives
             (test_vs_ectrain.tsv, self-hits removed, 0 = no hit);
             tie-averaged rank
4. Single-prediction comparators from EC-Bench results_test_100_go_3.csv:
   CLEAN-SupCon (shipped clean_supcon_correct) and EC-Bench's own blastp
   (correct if any predicted EC truncated to the level equals the true EC).

The llr gate for CAV "no prediction" is carried from the shipped file.
Writes figures/figure_data_staging/noknown/ec_ranks_noknown.csv.
"""
import glob
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
sys.path.insert(0, str(PAPER.parent / "bioCAV" / "core"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402

EC_SETS = Path("/xdisk/clairemcwhite/proFAB/EC_sets")
BENCH = PAPER / "ec_eval" / "static" / "ec_annotations" / "EC-Bench"
FIG_DATA = PAPER / "figures" / "figure_data"
STAGE = PAPER / "figures" / "figure_data_staging" / "noknown"
SCALER = PAPER / "shared_protein" / "reference_population" / "scaler_v1.pkl"
CAV_SUBDIR = "random_positive_train_max1000_cav"
SKIP = {"sp", "tr", "sw", "ref"}
HITS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]


def norm(pid):
    parts = [p for p in str(pid).split()[0].split("|") if p and p not in SKIP]
    return parts[0] if parts else str(pid)


def strip_wildcards(ec):
    return ".".join(p for p in str(ec).split(".") if p not in ("-", ""))


def read_span(path):
    return {norm(x) for x in path.read_text().split()} if path.exists() else set()


def main():
    ids, vecs, known = [], [], {}
    for d in sorted(glob.glob(str(EC_SETS / "ec_dataset_part*" / "ecNo_*"))):
        d = Path(d)
        cav = d / CAV_SUBDIR / "concept_v1.npy"
        if not cav.exists():
            continue
        ec = d.name[len("ecNo_"):].replace("-", ".")
        ids.append(ec)
        vecs.append(np.load(cav))
        known[ec] = (read_span(d / "random_positive_train_max1000.span") |
                     read_span(d / "random_positive_test_max1000.span"))
    level = np.array([len(e.split(".")) for e in ids])
    col = {e: i for i, e in enumerate(ids)}
    pos_of = {}
    for ec, ps in known.items():
        for p in ps:
            pos_of.setdefault(p, set()).add(ec)

    gold = pd.read_csv(BENCH / "test_ec_long.tsv", sep="\t")
    true_ecs = gold.assign(ec=gold["ec_number"].map(strip_wildcards)) \
                   .groupby("protein_id")["ec"].apply(set).to_dict()

    embs, sids = load_sequence_embeddings(BENCH / "test_ec.fasta.pkl")
    idx = {}
    for i, s in enumerate(sids):
        idx.setdefault(norm(s), i)
    S = joblib.load(SCALER).transform(embs) @ np.stack(vecs).T

    members = pd.read_csv(HERE / "ec_term_members.tsv", sep="\t", names=["ec", "sseqid"],
                          dtype=str)
    hits = pd.read_csv(HERE / "test_vs_ectrain.tsv", sep="\t", names=HITS)
    hits = hits[hits.qseqid != hits.sseqid]
    best = hits.groupby(["qseqid", "sseqid"], as_index=False)[["bitscore", "pident"]].max()
    qt = best.merge(members, on="sseqid").groupby(["qseqid", "ec"], as_index=False) \
             .agg(bitscore=("bitscore", "max"), pident=("pident", "max"))
    B = np.zeros_like(S)
    P = np.zeros_like(S)
    for q, e, b, pi in qt.itertuples(index=False):
        if q in idx and e in col:
            B[idx[q], col[e]] = b
            P[idx[q], col[e]] = pi

    tools = pd.read_csv(BENCH / "results_test_100_go_3.csv")
    tools.columns = [c.strip("﻿") for c in tools.columns]
    blastp = dict(zip(tools["id"], tools["blastp"].astype(str)))

    old = pd.read_csv(FIG_DATA / "ec_specificity_ranks.csv", dtype={"ec_number": str})
    rows = []
    for r in old.itertuples(index=False):
        p, ec = r.protein_id, r.ec_number
        lvl = len(ec.split("."))
        i, j = idx[p], col[ec]
        others = {".".join(e.split(".")[:lvl]) for e in true_ecs.get(p, ())
                  if len(e.split(".")) >= lvl}
        others |= {e for e in pos_of.get(p, ()) if len(e.split(".")) == lvl}
        others -= {ec}
        keep = level == lvl
        keep[[col[e] for e in others if e in col]] = False
        cs, bs = S[i][keep], B[i][keep]
        b_gt = int((bs > B[i, j]).sum())
        b_ties = int((bs == B[i, j]).sum()) - 1
        preds = [x.strip() for x in blastp.get(p, "").replace(",", ";").split(";")
                 if x.strip() and x != "nan"]
        rows.append({
            "protein_id": p, "ec_number": ec, "level": lvl, "llr": r.llr,
            "known_positive": p in known[ec],
            "cav_rank_level_filt_pf": int((cs > S[i, j]).sum()) + 1,
            "blast_rank_mid": b_gt + 1 + b_ties / 2,
            "blast_predicted": B[i, j] > 0,
            "blast_true_pident": P[i, j],
            "clean_supcon_correct": r.clean_supcon_correct,
            "ecbench_blastp_correct": any(".".join(x.split(".")[:lvl]) == ec for x in preds),
            "n_cand": int(keep.sum()),
        })
    out = pd.DataFrame(rows)
    STAGE.mkdir(parents=True, exist_ok=True)
    out.to_csv(STAGE / "ec_ranks_noknown.csv", index=False)

    print(f"{len(out)} pairs; known positives dropped: "
          f"{out.groupby('level').known_positive.sum().to_dict()}")
    k = out[~out.known_positive]
    for lvl, g in k.groupby("level"):
        cav1 = ((g.cav_rank_level_filt_pf == 1) & (g.llr > 0)).mean()
        cav3 = ((g.cav_rank_level_filt_pf <= 3) & (g.llr > 0)).mean()
        b1 = ((g.blast_rank_mid < 1.5) & g.blast_predicted).mean()
        b3 = ((g.blast_rank_mid <= 3) & g.blast_predicted).mean()
        print(f"  level {lvl} n={len(g)} | CAV rank1 {100*cav1:.1f}% top3 {100*cav3:.1f}% | "
              f"BLAST(ours) rank1 {100*b1:.1f}% top3 {100*b3:.1f}% | "
              f"CLEAN {100*g.clean_supcon_correct.mean():.1f}% | "
              f"EC-Bench blastp {100*g.ecbench_blastp_correct.mean():.1f}% | "
              f"no BLAST hit to true EC {100*(~g.blast_predicted).mean():.0f}%")


if __name__ == "__main__":
    main()
