#!/usr/bin/env python
"""
Does a PCA step before CAV training help GO CC ranking? (exploratory, 2026-09-29)
GO counterpart of ec_eval/pca_test/pca_test.py.

Arms, all over the 707 CC GO CAVs:
  prod        production CAVs on disk (all <=1000 positives, no PCA)
  sub_nopca   retrained on <=30 sampled positives per CAV (cc_subsample.fasta.pkl),
              balanced LR C=1 lbfgs seed 42, shared scaler
  sub_pca{k}  same positives, LR in the top-k PCs of the 10k standardised shared
              negatives (exact PCA), direction folded back to 1152-d

Evaluation: rank of the true term among CC CAVs, with the strict (_pf) filter of
rank_all_terms/compute_ranks_all_terms.py: minus the true term's ancestors, the
protein's GPAD terms and its proFAB-positive terms (both ancestor-closed).
  temporal     the 4,055 known-positive-excluded CC validation pairs (Fig 2)
  profab_test  held-out proFAB test positives, capped at 10 per term, for the 258
               terms whose test embeddings are cached (all 707 remain candidates)
Plus per-term AUROC of test positives vs test negatives. Ranks are ungated.
Writes pca_test_cc_summary.tsv, pca_test_cc_ranks.tsv, pca_test_cc_per_term_auroc.tsv.
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
sys.path.insert(0, str(HERE.parent / "rank_all_terms"))
from compute_ranks_all_terms import (  # noqa: E402
    CAV_SUBDIR, GO_EVAL, GO_SETS, OBO, ONTS, ancestors, load_profab_positives,
    normalise_pid, parse_obo,
)
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402  (path set by import above)

REF = PAPER / "shared_protein" / "reference_population"
STAGE = PAPER / "figures" / "figure_data_staging" / "noknown"
KS = [32, 64, 128, 256, 512]
CAP = 10


def load_pkl(path):
    E, ids = load_sequence_embeddings(str(path))
    return np.asarray(E, dtype=np.float64), [normalise_pid(s) for s in ids]


def fit_dir(Xpos, Xneg):
    X = np.vstack([Xpos, Xneg])
    y = np.r_[np.ones(len(Xpos)), np.zeros(len(Xneg))]
    lr = LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs",
                            max_iter=1000, random_state=42).fit(X, y)
    return lr.coef_.ravel()


def main():
    ns, parents = parse_obo(OBO)
    scaler = joblib.load(REF / "scaler_v1.pkl")
    Zneg = scaler.transform(load_pkl(REF / "neg_10000.pkl")[0])

    dirs, terms = {}, []
    for d in sorted(glob.glob(str(GO_SETS / "go_dataset_part*" / "GO_*"))):
        d = Path(d)
        g = d.name.replace("GO_", "GO:")
        if ns.get(g) == "cc" and (d / CAV_SUBDIR / "concept_v1.npy").exists():
            terms.append(g)
            dirs[g] = d
    col = {g: i for i, g in enumerate(terms)}
    print(f"{len(terms)} CC CAVs", flush=True)

    sub, sub_ids = load_pkl(HERE / "cc_subsample.fasta.pkl")
    sidx = {p: i for i, p in enumerate(sub_ids)}
    Zsub = scaler.transform(sub)
    members = pd.read_csv(HERE / "cc_subsample_members.tsv", sep="\t", names=["g", "p"])
    by_term = members.groupby("g").p.apply(list).to_dict()

    W = {"prod": np.stack([np.load(dirs[g] / CAV_SUBDIR / "concept_v1.npy").ravel()
                           for g in terms], axis=1)}
    pcas = {k: PCA(n_components=k, svd_solver="full").fit(Zneg) for k in KS}
    Zneg_k = {k: p.transform(Zneg) for k, p in pcas.items()}
    for arm in ["sub_nopca"] + [f"sub_pca{k}" for k in KS]:
        W[arm] = np.zeros((Zneg.shape[1], len(terms)))
    for j, g in enumerate(terms):
        Xp = Zsub[[sidx[p] for p in by_term.get(g, []) if p in sidx]]
        W["sub_nopca"][:, j] = fit_dir(Xp, Zneg)
        for k, pca in pcas.items():
            W[f"sub_pca{k}"][:, j] = pca.components_.T @ fit_dir(pca.transform(Xp), Zneg_k[k])
        if (j + 1) % 100 == 0:
            print(f"  trained {j+1}/{len(terms)}", flush=True)
    for a in W:
        W[a] = W[a] / np.linalg.norm(W[a], axis=0, keepdims=True)

    # known-true terms per protein: GPADs (any date) and proFAB positives, ancestor-closed
    gpad = pd.concat([
        pd.read_csv(GO_EVAL / "static" / "go_annotations" / f"{o}.gpad", sep="\t",
                    comment="!", header=None, usecols=[0, 1, 3], names=["db", "pid", "go"])
        for o in ONTS])
    gpad = gpad[gpad["db"] == "UniProtKB"]
    gpad["pid"] = gpad["pid"].astype(str).str.split("-").str[0]
    gpad_terms = gpad.groupby("pid")["go"].apply(set).to_dict()

    rng = np.random.default_rng(42)
    ev, per_term = [], []
    for g in terms:
        # held-out test embeddings are cached only for the terms Fig 2 evaluated
        # (258 of 707); the others still compete as candidates
        if not (dirs[g] / "random_negative_test_max1000.span.fasta.pkl").exists() or \
                not (dirs[g] / "random_positive_test_max1000.span.fasta.pkl").exists():
            continue
        Ep, ip = load_pkl(dirs[g] / "random_positive_test_max1000.span.fasta.pkl")
        En, _ = load_pkl(dirs[g] / "random_negative_test_max1000.span.fasta.pkl")
        Zp, Zn = scaler.transform(Ep), scaler.transform(En)
        per_term.append((g, Zp, Zn))
        for i in rng.permutation(len(ip))[:CAP]:
            ev.append(("profab_test", ip[i], g, Zp[i]))
    pairs = pd.read_csv(STAGE / "go_specificity_ranks_allterms_cc.csv")
    Ev, iv = load_pkl(GO_EVAL / "results" / "temporal_eval_cc" / "eval_embeddings" /
                      "val_proteins.span.fasta.pkl")
    vidx = {}
    for i, p in enumerate(iv):
        vidx.setdefault(p, i)
    Zv = scaler.transform(Ev)
    for r in pairs.itertuples(index=False):
        ev.append(("temporal", r.protein_id, r.go_term, Zv[vidx[r.protein_id]]))
    profab_pos = load_profab_positives({p for _, p, _, _ in ev})
    print(f"{len(ev)} evaluation pairs", flush=True)

    drop_cache = {}
    out = []
    for src, p, g, z in ev:
        key = (p, g)
        if key not in drop_cache:
            drop = ancestors(g, parents)
            for t in gpad_terms.get(p, set()) | profab_pos.get(p, set()):
                drop |= {t} | ancestors(t, parents)
            drop -= {g}
            keep = np.ones(len(terms), bool)
            keep[[col[t] for t in drop if t in col]] = False
            drop_cache[key] = keep
        keep = drop_cache[key]
        for arm, M in W.items():
            s = z @ M
            out.append((arm, src, p, g, int((s[keep] > s[col[g]]).sum()) + 1))
    res = pd.DataFrame(out, columns=["arm", "source", "protein_id", "go_term", "rank"])

    au = []
    for arm, M in W.items():
        for g, Zp, Zn in per_term:
            j = col[g]
            y = np.r_[np.ones(len(Zp)), np.zeros(len(Zn))]
            au.append((arm, g, roc_auc_score(y, np.r_[Zp @ M[:, j], Zn @ M[:, j]])))
    au = pd.DataFrame(au, columns=["arm", "go_term", "auroc"])

    summ = res.groupby(["source", "arm"]).agg(
        n=("rank", "size"), rank1=("rank", lambda r: 100 * (r == 1).mean()),
        top3=("rank", lambda r: 100 * (r <= 3).mean()),
        top10=("rank", lambda r: 100 * (r <= 10).mean()),
        median_rank=("rank", "median")).reset_index()
    summ = summ.merge(au.groupby("arm").auroc.median().rename("median_auroc").reset_index(), on="arm")
    cos = {a: float(np.median(np.sum(W[a] * W["prod"], axis=0))) for a in W}
    summ["median_cos_to_prod"] = summ.arm.map(cos)
    order = ["prod", "sub_nopca"] + [f"sub_pca{k}" for k in KS]
    summ["arm"] = pd.Categorical(summ.arm, order)
    summ = summ.sort_values(["source", "arm"])
    summ.to_csv(HERE / "pca_test_cc_summary.tsv", sep="\t", index=False, float_format="%.3f")
    res.to_csv(HERE / "pca_test_cc_ranks.tsv", sep="\t", index=False)
    au.to_csv(HERE / "pca_test_cc_per_term_auroc.tsv", sep="\t", index=False, float_format="%.4f")
    pd.set_option("display.width", 200)
    print(summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
