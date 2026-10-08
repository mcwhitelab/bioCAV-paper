#!/usr/bin/env python
"""
Pfam CAVs, MYC PCA recipe at k=128, with independent (per-motif) scalers.  Exploratory
(2026-09-29): does a PCA step restore motif clustering after the scaler change?

Per motif (same embeddings and truncation filter as ../retrain_pfam/retrain_pfam_cavs.py):
  1. --scaler own (default): StandardScaler fit on that motif's positives + 10k negatives
     (the original library's per-motif scaler, i.e. Ahmad's create_scaler_and_pca without PCA).
     --scaler shared: the shared reference scaler (scaler_v1.pkl, 10k negatives only), as in
     the shared-scaler retrain; the control arm.
  2. Exact PCA (svd_solver=full), k=128, fit on the 10k standardized negatives only
     (the MYC recipe, interface_search/.../pca_sweep_mad1/build_pca_sweep.py).
  3. Balanced LR, C=1, lbfgs, seed 42, max_iter 5000, in PC space.  No CV.
  4. Direction folded back to 1152-d (components_.T @ coef), unit-normalised.
With --scaler own the CAV lives in the motif's own frame, so its scaler is saved beside it.
Writes cavs/<PF>/L25_concept_v1.npy (cav_cluster.py layout), L25_scaler_v1.pkl and
L25_pca128_report.json (train AUROC, variance explained, cosine to old + retrained CAV).

Embeddings are read in place; nothing is copied.  Idempotent (skips motifs with a report).
Usage: train_pfam_pca128.py --shard I --n-shards N [--scaler own|shared] [--out DIR] [--motifs ...]
"""
import argparse, json, sys, time, traceback
from pathlib import Path

import joblib
import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

SRC = Path("/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all")
EMB = SRC / "embeddings" / "ESMplusplus_large"
SHARED = Path("/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_shared_scaler/cavs")
DEFAULT_OUT = {"own": Path("/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_pca128_ownscaler"),
               "shared": Path("/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_pca128_sharedscaler")}
SHARED_SCALER = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper/"
                     "shared_protein/reference_population/scaler_v1.pkl")
K, LAYER = 128, 25
TOKEN_CAP_RESIDUES = 1998  # max_length=2000 tokens minus BOS/EOS


def load(motif):
    d = EMB / motif
    pos, neg = np.load(d / f"L{LAYER}_pos.npy"), np.load(d / f"L{LAYER}_neg.npy")
    assert np.isfinite(neg).all()
    meta = json.load(open(d / "pos" / f"L{LAYER}_meta.json"))
    assert len(meta) == len(pos)
    keep = np.isfinite(pos).all(axis=1) & np.array([m["window_span"][0] < TOKEN_CAP_RESIDUES for m in meta])
    return pos[keep], neg, int((~keep).sum())


def cos_to(path, w):
    if not path.exists():
        return None
    r = np.load(path).ravel()
    return float(w @ r / np.linalg.norm(r))


def train_one(motif, out, shared_scaler=None):
    d = out / "cavs" / motif
    rep = d / f"L{LAYER}_pca128_report.json"
    if rep.exists():
        return "skip"
    pos, neg, n_drop = load(motif)
    X = np.vstack([pos, neg]).astype(np.float64)
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    scaler = shared_scaler if shared_scaler is not None else StandardScaler().fit(X)
    Xs = scaler.transform(X)
    pca = PCA(n_components=K, svd_solver="full").fit(Xs[len(pos):])
    Z = pca.transform(Xs)
    lr = LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs",
                            max_iter=5000, random_state=42).fit(Z, y)
    w = pca.components_.T @ lr.coef_.ravel()
    w /= np.linalg.norm(w)
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / f"L{LAYER}_concept_v1.npy", w.astype(np.float32))
    if shared_scaler is None:
        joblib.dump(scaler, d / f"L{LAYER}_scaler_v1.pkl")
    json.dump({
        "motif": motif, "k": K, "scaler": "own" if shared_scaler is None else str(SHARED_SCALER), "n_pos": len(pos), "n_neg": len(neg),
        "n_dropped_truncated_pos": n_drop,
        "pca_var_explained": float(pca.explained_variance_ratio_.sum()),
        "train_auroc": float(roc_auc_score(y, lr.decision_function(Z))),
        "lr_n_iter": int(lr.n_iter_[0]),
        "cos_to_old": cos_to(SRC / "cavs" / motif / f"L{LAYER}_concept_v1.npy", w),
        "cos_to_shared_retrain": cos_to(SHARED / motif / f"L{LAYER}_concept_v1.npy", w),
        "date": time.strftime("%Y-%m-%d"),
    }, open(rep, "w"), indent=1)
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--scaler", choices=["own", "shared"], default="own")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--motifs", nargs="*")
    a = ap.parse_args()
    a.out = a.out or DEFAULT_OUT[a.scaler]
    shared = None
    if a.scaler == "shared":
        shared = joblib.load(SHARED_SCALER)
        assert shared.n_samples_seen_ == 10000, "expected the negatives-only shared scaler"
    motifs = a.motifs or sorted(p.name for p in EMB.iterdir() if p.name.startswith("PF"))
    motifs = motifs[a.shard::a.n_shards]
    print(f"shard {a.shard}/{a.n_shards}: {len(motifs)} motifs -> {a.out}", flush=True)
    fails = a.out / "failed" / f"shard_{a.shard}.tsv"
    fails.parent.mkdir(parents=True, exist_ok=True)
    counts, t0 = {"ok": 0, "skip": 0, "fail": 0}, time.time()
    for i, m in enumerate(motifs, 1):
        try:
            counts[train_one(m, a.out, shared)] += 1
        except Exception as e:
            counts["fail"] += 1
            with open(fails, "a") as fh:
                fh.write(f"{m}\t{type(e).__name__}: {str(e).splitlines()[0] if str(e) else ''}\n")
            traceback.print_exc()
        if i % 50 == 0:
            print(f"{i}/{len(motifs)} {counts} {(time.time() - t0) / i:.1f}s/motif", flush=True)
    print(f"done {counts}", flush=True)


if __name__ == "__main__":
    main()
