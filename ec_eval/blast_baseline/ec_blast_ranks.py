#!/usr/bin/env python
"""
Homology baselines for Fig 3A (2026-09-29, exploratory), on the same 209
(protein, EC) pairs as figure_data/ec_specificity_ranks_samelevel.csv.

v1 — EC-Bench's published blastp predictions (results_test_100_go_3.csv,
     column "blastp"), scored exactly like clean_supcon_correct
     (compare_ec_tools.matches_exact: match at the true EC's own level).
v2 — DIAMOND nearest-neighbour ranking, mirroring go_eval/blast_baseline/
     blast_ranks.py: score of an EC = best bitscore of the protein to that EC's
     CAV training positives (random_positive_train_max1000.span.fasta; self-hits
     removed; --more-sensitive, e <= 1e-3). The true EC is ranked among the
     same candidates as cav_rank_level_filt (same EC level, minus the protein's
     other true ECs at that level), further excluding EC numbers deleted or
     transferred in enzyme.dat. blast_rank_opt = 1 + #{score > true};
     blast_rank_mid adds half the ties. "No prediction" = no hit to the true EC.

Also recomputes cav_rank_level_live = cav_rank_level_filt without the deleted/
transferred ECs, for the Fig 3A preview. Writes ec_blast_ranks.csv here.
"""
import glob
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
sys.path.insert(0, str(PAPER / "ec_eval" / "rank_same_level"))
sys.path.insert(0, str(PAPER.parent / "bioCAV" / "specific_scripts"))
import compute_ec_ranks_same_level as E  # noqa: E402
from compare_ec_tools import matches_exact  # noqa: E402

ENZYME_DAT = PAPER / "ec_eval" / "static" / "ec_annotations" / "enzyme.dat"
TOOLS = PAPER / "ec_eval" / "static" / "ec_annotations" / "EC-Bench" / "results_test_100_go_3.csv"
HITS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]


def dead_ecs():
    dead, cur = set(), None
    for line in open(ENZYME_DAT):
        if line.startswith("ID   "):
            cur = line[5:].strip()
        elif line.startswith("DE   ") and cur and ("Transferred entry" in line or "Deleted entry" in line):
            dead.add(cur)
    return dead


def main():
    ids, vecs = [], []
    for d in sorted(glob.glob(str(E.EC_SETS / "ec_dataset_part*" / "ecNo_*"))):
        cav = Path(d) / E.CAV_SUBDIR / "concept_v1.npy"
        if cav.exists():
            ids.append(Path(d).name[len("ecNo_"):].replace("-", "."))
            vecs.append(np.load(cav))
    col = {e: i for i, e in enumerate(ids)}
    level = np.array([len(e.split(".")) for e in ids])
    dead = dead_ecs()
    live = np.array([e not in dead for e in ids])
    print("deleted/transferred EC CAVs:", [e for e in ids if e in dead])

    gold = pd.read_csv(E.GOLD, sep="\t")
    true_ecs = gold.assign(ec=gold["ec_number"].map(E.strip_wildcards)) \
                   .groupby("protein_id")["ec"].apply(set).to_dict()
    embs, sids = E.load_sequence_embeddings(E.VAL_PKL)
    idx = {}
    for i, s in enumerate(sids):
        idx.setdefault(E.normalise_pid(s), i)
    S = joblib.load(E.SCALER).transform(embs) @ np.stack(vecs).T

    members = pd.read_csv(HERE / "ec_term_members.tsv", sep="\t", names=["ec", "sseqid"])
    hits = pd.read_csv(HERE / "test_vs_ec_train.tsv", sep="\t", names=HITS)
    hits = hits[hits.qseqid != hits.sseqid]
    qt = hits.groupby(["qseqid", "sseqid"], as_index=False)[["bitscore", "pident"]].max() \
             .merge(members, on="sseqid") \
             .groupby(["qseqid", "ec"], as_index=False).agg(bits=("bitscore", "max"), pident=("pident", "max"))
    by_q = {q: d for q, d in qt.groupby("qseqid")}

    tools = pd.read_csv(TOOLS, dtype=str).fillna("-").rename(columns={"id": "protein_id"})
    blastp = tools.drop_duplicates("protein_id").set_index("protein_id")["blastp"]

    pairs = pd.read_csv(E.FIG_DATA / "ec_specificity_ranks_samelevel.csv", dtype={"ec_number": str})
    rows = []
    for r in pairs.itertuples(index=False):
        p, ec = r.protein_id, r.ec_number
        lvl = len(ec.split("."))
        others = {".".join(e.split(".")[:lvl]) for e in true_ecs.get(p, ())
                  if len(e.split(".")) >= lvl} - {ec}
        keep = (level == lvl) & live
        keep[[col[e] for e in others if e in col]] = False
        j = col[ec]
        s = S[idx[p]]
        b = np.zeros(len(ids)); pid_ = np.zeros(len(ids))
        if p in by_q:
            d = by_q[p][by_q[p].ec.isin(col)]
            b[d.ec.map(col).values] = d.bits.values
            pid_[d.ec.map(col).values] = d.pident.values
        cand = b[keep]
        gt, ties = int((cand > b[j]).sum()), int((cand == b[j]).sum()) - 1
        rows.append({
            "protein_id": p, "ec_number": ec, "level": r.level,
            "cav_rank_level_live": int((s[keep] > s[j]).sum()) + 1,
            "blast_true_bits": b[j], "blast_true_pident": pid_[j], "blast_predicted": b[j] > 0,
            "blast_rank_opt": gt + 1, "blast_rank_mid": gt + 1 + ties / 2, "blast_n_ties": ties,
            "ecbench_blastp_pred": blastp.get(p, "-"),
            "ecbench_blastp_correct": matches_exact(blastp.get(p, "-"), ec),
            "n_cand": int(keep.sum()),
        })
    out = pairs[["protein_id", "ec_number", "llr", "clean_supcon_correct", "cav_rank_level_filt"]] \
        .merge(pd.DataFrame(rows), on=["protein_id", "ec_number"])
    out.to_csv(HERE / "ec_blast_ranks.csv", index=False)
    print(f"{len(out)} pairs; CAV rank changed by live mask: "
          f"{(out.cav_rank_level_live != out.cav_rank_level_filt).sum()}")
    for lvl, g in out.groupby("level"):
        pr = g.blast_predicted
        f = lambda m: f"{100 * m.mean():.1f}"
        print(f"  level {lvl} n={len(g)} | CAV r1 {f((g.llr > 0) & (g.cav_rank_level_live == 1))} "
              f"top10 {f((g.llr > 0) & (g.cav_rank_level_live <= 10))} | DIAMOND hit {f(pr)} "
              f"r1 mid {f(pr & (g.blast_rank_mid < 1.5))} opt {f(pr & (g.blast_rank_opt == 1))} "
              f"top10 mid {f(pr & (g.blast_rank_mid <= 10))} | EC-Bench blastp correct "
              f"{f(g.ecbench_blastp_correct)} | CLEAN-SupCon correct {f(g.clean_supcon_correct)}")


if __name__ == "__main__":
    main()
