#!/usr/bin/env python
"""
Does a PCA step before CAV training help EC level-3 ranking? (exploratory, 2026-09-29)

Arms, all over the 153 EC level-3 CAVs:
  prod        production CAVs on disk (all <=1000 positives, no PCA)
  sub_nopca   retrained on the same <=100 sampled positives per CAV
              (l3_subsample.fasta.pkl), balanced LR C=1 lbfgs seed 42, shared scaler
  sub_pca{k}  as sub_nopca, but LR in the top-k PCs of the 10k standardised shared
              negatives (exact PCA), direction folded back to 1152-d
              (the MYC recipe: build_pca_sweep.py)

Evaluation, rank of the true EC among the level-3 CAVs; strict filter (other ECs the
protein is a proFAB positive for, or carries in EC-Bench, are removed):
  profab_test  held-out proFAB random_positive_test positives (disjoint from the
               training split), capped at 50 per EC, for the 33 ECs whose test
               embeddings are cached (all 153 remain candidates)
  ecbench      the 107 leak-free EC-Bench level-3 pairs
Plus per-EC AUROC of the CAV's test positives vs its test negatives.
Ranks are ungated (no LLR "no prediction" rule), so they are comparable across arms.
Writes pca_test_summary.tsv and pca_test_per_ec.tsv.
"""
import glob
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
sys.path.insert(0, str(PAPER.parent / "bioCAV" / "core"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402

EC_SETS = Path("/xdisk/clairemcwhite/proFAB/EC_sets")
REF = PAPER / "shared_protein" / "reference_population"
BENCH = PAPER / "ec_eval" / "static" / "ec_annotations" / "EC-Bench"
STAGE = PAPER / "figures" / "figure_data_staging" / "noknown"
KS = [32, 64, 128, 256, 512]
CAP = 50
SKIP = {"sp", "tr", "sw", "ref"}


def norm(pid):
    parts = [p for p in str(pid).split()[0].split("|") if p and p not in SKIP]
    return parts[0] if parts else str(pid)


def load_pkl(path):
    E, ids = load_sequence_embeddings(str(path))
    return np.asarray(E, dtype=np.float64), [norm(s) for s in ids]


def read_span(path):
    return {norm(x) for x in path.read_text().split()} if path.exists() else set()


def fit_dir(Xpos, Xneg):
    X = np.vstack([Xpos, Xneg])
    y = np.r_[np.ones(len(Xpos)), np.zeros(len(Xneg))]
    lr = LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs",
                            max_iter=1000, random_state=42).fit(X, y)
    return lr.coef_.ravel()


def main():
    scaler = joblib.load(REF / "scaler_v1.pkl")
    neg, _ = load_pkl(REF / "neg_10000.pkl")
    Zneg = scaler.transform(neg)

    dirs, ecs = {}, []
    for d in sorted(glob.glob(str(EC_SETS / "ec_dataset_part*" / "ecNo_*"))):
        d = Path(d)
        ec = d.name[len("ecNo_"):]
        if ec.count("-") == 2 and (d / "random_positive_train_max1000_cav" / "concept_v1.npy").exists():
            ecs.append(ec.replace("-", "."))
            dirs[ecs[-1]] = d
    col = {e: i for i, e in enumerate(ecs)}
    print(f"{len(ecs)} level-3 EC CAVs", flush=True)

    # proFAB positives (train+test) for every level-3 EC, for the strict filter
    known = {e: read_span(dirs[e] / "random_positive_train_max1000.span") |
                read_span(dirs[e] / "random_positive_test_max1000.span") for e in ecs}
    pos_of = {}
    for e, ps in known.items():
        for p in ps:
            pos_of.setdefault(p, set()).add(e)

    # training subsample embeddings
    sub, sub_ids = load_pkl(HERE / "l3_subsample.fasta.pkl")
    sidx = {p: i for i, p in enumerate(sub_ids)}
    members = pd.read_csv(HERE / "l3_subsample_members.tsv", sep="\t", names=["ec", "p"])
    Zsub = scaler.transform(sub)

    # CAV matrices, one column per EC
    W = {"prod": np.stack([np.load(dirs[e] / "random_positive_train_max1000_cav" /
                                   "concept_v1.npy").ravel() for e in ecs], axis=1)}
    pcas = {k: PCA(n_components=k, svd_solver="full").fit(Zneg) for k in KS}
    for arm in ["sub_nopca"] + [f"sub_pca{k}" for k in KS]:
        W[arm] = np.zeros((Zneg.shape[1], len(ecs)))
    for j, e in enumerate(ecs):
        rows = [sidx[p] for p in members.loc[members.ec == e, "p"] if p in sidx]
        Xp = Zsub[rows]
        W["sub_nopca"][:, j] = fit_dir(Xp, Zneg)
        for k, pca in pcas.items():
            w = pca.components_.T @ fit_dir(pca.transform(Xp), pca.transform(Zneg))
            W[f"sub_pca{k}"][:, j] = w
        if (j + 1) % 25 == 0:
            print(f"  trained {j+1}/{len(ecs)}", flush=True)
    for a in W:
        W[a] = W[a] / np.linalg.norm(W[a], axis=0, keepdims=True)

    # evaluation sets
    rng = np.random.default_rng(42)
    ev_rows, per_ec = [], []
    test_neg = {}
    for e in ecs:
        # held-out test embeddings are cached only for the ECs the original eval
        # scored (33 of 153); the others still compete as candidates
        if not (dirs[e] / "random_positive_test_max1000.span.fasta.pkl").exists() or \
                not (dirs[e] / "random_negative_test_max1000.span.fasta.pkl").exists():
            continue
        Ep, ip = load_pkl(dirs[e] / "random_positive_test_max1000.span.fasta.pkl")
        En, _ = load_pkl(dirs[e] / "random_negative_test_max1000.span.fasta.pkl")
        test_neg[e] = scaler.transform(En)
        keep = rng.permutation(len(ip))[:CAP]
        Zp = scaler.transform(Ep)
        per_ec.append((e, Zp, test_neg[e]))
        for i in keep:
            ev_rows.append(("profab_test", ip[i], e, Zp[i]))
    ecb = pd.read_csv(STAGE / "ec_ranks_noknown.csv", dtype={"ec_number": str})
    ecb = ecb[(~ecb.known_positive) & (ecb.level == 3)]
    Eb, ib = load_pkl(BENCH / "test_ec.fasta.pkl")
    bidx = {p: i for i, p in enumerate(ib)}
    Zb = scaler.transform(Eb)
    gold = pd.read_csv(BENCH / "test_ec_long.tsv", sep="\t")
    bench_true = gold.assign(e3=gold.ec_number.astype(str).str.split(".").str[:3].str.join(".")) \
                     .groupby("protein_id").e3.apply(set).to_dict()
    for r in ecb.itertuples(index=False):
        ev_rows.append(("ecbench", r.protein_id, r.ec_number, Zb[bidx[r.protein_id]]))

    out = []
    for arm, M in W.items():
        for src, p, e, z in ev_rows:
            s = z @ M
            drop = (pos_of.get(p, set()) | bench_true.get(p, set())) - {e}
            keep = np.ones(len(ecs), bool)
            keep[[col[x] for x in drop if x in col]] = False
            rank = int((s[keep] > s[col[e]]).sum()) + 1
            out.append((arm, src, p, e, rank))
    res = pd.DataFrame(out, columns=["arm", "source", "protein_id", "ec", "rank"])

    aurocs = []
    for arm, M in W.items():
        for e, Zp, Zn in per_ec:
            j = col[e]
            y = np.r_[np.ones(len(Zp)), np.zeros(len(Zn))]
            aurocs.append((arm, e, roc_auc_score(y, np.r_[Zp @ M[:, j], Zn @ M[:, j]])))
    au = pd.DataFrame(aurocs, columns=["arm", "ec", "auroc"])

    summ = res.groupby(["source", "arm"]).agg(
        n=("rank", "size"), rank1=("rank", lambda r: 100 * (r == 1).mean()),
        top3=("rank", lambda r: 100 * (r <= 3).mean()),
        top10=("rank", lambda r: 100 * (r <= 10).mean()),
        median_rank=("rank", "median")).reset_index()
    macro = res[res.source == "profab_test"].assign(r1=lambda d: d["rank"] == 1) \
        .groupby(["arm", "ec"]).r1.mean().groupby("arm").mean().mul(100).rename("profab_macro_rank1")
    summ = summ.merge(macro.reset_index(), on="arm", how="left") \
               .merge(au.groupby("arm").auroc.median().rename("median_auroc").reset_index(), on="arm")
    cos = {a: float(np.median(np.sum(W[a] * W["prod"], axis=0))) for a in W}
    summ["median_cos_to_prod"] = summ.arm.map(cos)
    order = ["prod", "sub_nopca"] + [f"sub_pca{k}" for k in KS]
    summ["arm"] = pd.Categorical(summ.arm, order)
    summ = summ.sort_values(["source", "arm"])
    summ.to_csv(HERE / "pca_test_summary.tsv", sep="\t", index=False, float_format="%.3f")
    res.to_csv(HERE / "pca_test_ranks.tsv", sep="\t", index=False)
    au.to_csv(HERE / "pca_test_per_ec_auroc.tsv", sep="\t", index=False, float_format="%.4f")
    pd.set_option("display.width", 200)
    print(summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
