#!/usr/bin/env python
"""
Pfam CAVs trained to convergence (2026-09-29).

The original library and both retrains used sklearn's default lbfgs tol=1e-4, which
stops after 16-50 iterations on these nearly separable problems, far from the optimum
(converged vs default cos ~0.7).  This retrains every motif with tol=1e-10 and a large
max_iter so the direction is the unique L2-regularised optimum and does not depend on
the optimiser path.

Per motif (same embeddings and truncation filter as ../retrain_pfam_pca128/train_pfam_pca128.py):
  1. --scaler shared (default): the shared reference scaler (scaler_v1.pkl, 10k negatives
     only).  --scaler own: a StandardScaler fit on that motif's positives + 10k negatives,
     as in the original library; it is saved beside the CAV.  At convergence the two agree
     closely (own vs shared cos 0.995-0.999 on PF00001/PF00010/PF06943); both are run so
     the clustering can be compared.
  2. Balanced LR, C=1, lbfgs, tol=1e-10, max_iter=100000, seed 42, full 1152-d, no PCA.
  3. coef_ unit-normalised.
Writes cavs/<PF>/L25_concept_v1.npy (cav_cluster.py layout) and L25_converged_report.json
(n_iter, hit_max_iter, train AUROC, cosine to the old and shared-scaler-retrain CAVs).

Embeddings are read in place; nothing is copied.  Idempotent (skips motifs with a report).
Converged means it: lbfgs at tol=1e-10 matches an exact newton-cholesky solve to cos
0.9999999997 (PF00001/PF00010/PF00069).
C sweep (2026-09-29): converged C=1 clusters clans far worse than the early-stopped
library, so --C tests whether stronger L2 recovers it as a reproducible optimum.
Usage: train_pfam_converged.py --shard I --n-shards N [--scaler shared|own] [--C 1e-3] [--out DIR] [--motifs ...]
"""
import argparse, json, sys, time, traceback
from pathlib import Path

import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "retrain_pfam_pca128"))
from train_pfam_pca128 import load, cos_to, SRC, EMB, SHARED, SHARED_SCALER, LAYER  # noqa: E402

DEFAULT_OUT = {"shared": Path("/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_converged_sharedscaler"),
               "own": Path("/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_converged_ownscaler")}
TOL, MAX_ITER = 1e-10, 100000


def train_one(motif, out, shared_scaler=None, C=1.0):
    d = out / "cavs" / motif
    rep = d / f"L{LAYER}_converged_report.json"
    if rep.exists():
        return "skip"
    pos, neg, n_drop = load(motif)
    X = np.vstack([pos, neg]).astype(np.float64)
    scaler = shared_scaler if shared_scaler is not None else StandardScaler().fit(X)
    X = scaler.transform(X)
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    t0 = time.time()
    lr = LogisticRegression(C=C, class_weight="balanced", solver="lbfgs", tol=TOL,
                            max_iter=MAX_ITER, random_state=42).fit(X, y)
    w = lr.coef_.ravel() / np.linalg.norm(lr.coef_)
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / f"L{LAYER}_concept_v1.npy", w.astype(np.float32))
    if shared_scaler is None:
        joblib.dump(scaler, d / f"L{LAYER}_scaler_v1.pkl")
    n_iter = int(lr.n_iter_[0])
    json.dump({
        "motif": motif, "scaler": "own" if shared_scaler is None else str(SHARED_SCALER), "C": C, "tol": TOL, "max_iter": MAX_ITER,
        "n_iter": n_iter, "hit_max_iter": n_iter >= MAX_ITER, "fit_seconds": round(time.time() - t0, 2),
        "n_pos": len(pos), "n_neg": len(neg), "n_dropped_truncated_pos": n_drop,
        "train_auroc": float(roc_auc_score(y, lr.decision_function(X))),
        "cos_to_old": cos_to(SRC / "cavs" / motif / f"L{LAYER}_concept_v1.npy", w),
        "cos_to_shared_retrain": cos_to(SHARED / motif / f"L{LAYER}_concept_v1.npy", w),
        "date": time.strftime("%Y-%m-%d"),
    }, open(rep, "w"), indent=1)
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--scaler", choices=["shared", "own"], default="shared")
    ap.add_argument("--C", type=float, default=1.0,
                    help="L2 inverse regularisation; C != 1 writes to ..._C<value>_<scaler>scaler by default")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--motifs", nargs="*")
    a = ap.parse_args()
    if a.out is None:
        a.out = DEFAULT_OUT[a.scaler]
        if a.C != 1.0:
            a.out = a.out.with_name(f"pfam_cavs_converged_C{a.C:g}_{a.scaler}scaler")
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
            counts[train_one(m, a.out, shared, a.C)] += 1
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
