#!/usr/bin/env python3
"""Refit every single-cell CAV in the b617ee1b library three ways, in the
library's own PCA-128 frame (per-CAV scaler_v1/pca_v1, which are the shared
global ones):
  orig   : C=1, lbfgs, tol=1e-4 (default), max_iter=1000, balanced  (stored recipe)
  conv1  : C=1, tol=1e-10, max_iter=1e5
  conv4  : C=1e-4, tol=1e-10, max_iter=1e5   (library decision, C*n ~ 1)
Writes unit concept vectors to cavs_<arm>/<name>/concept_v1.npy and a summary TSV.
"""
import sys, json, warnings
from pathlib import Path
import numpy as np, pandas as pd, joblib
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning

LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
ARMS = {"orig": dict(C=1.0, tol=1e-4, max_iter=1000),
        "conv1": dict(C=1.0, tol=1e-10, max_iter=100000),
        "conv4": dict(C=1e-4, tol=1e-10, max_iter=100000)}
rows = []
for d in sorted((LIB / "cavs").iterdir()):
    name = d.name
    stored = np.load(d / "concept_v1.npy").astype(np.float64)
    stored /= np.linalg.norm(stored)
    scaler, pca = joblib.load(d / "scaler_v1.pkl"), joblib.load(d / "pca_v1.pkl")
    pos, neg = np.load(d / "pos.npy"), np.load(d / "neg.npy")
    X = pca.transform(scaler.transform(np.vstack([pos, neg])))
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    r = {"cav": name, "n_pos": len(pos), "n_neg": len(neg)}
    vecs = {}
    for arm, kw in ARMS.items():
        clf = LogisticRegression(solver="lbfgs", class_weight="balanced", random_state=42, **kw)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", ConvergenceWarning)
            clf.fit(X, y)
        v = clf.coef_[0] / np.linalg.norm(clf.coef_[0])
        vecs[arm] = v
        r[f"n_iter_{arm}"] = int(clf.n_iter_[0])
        r[f"warn_{arm}"] = any(issubclass(x.category, ConvergenceWarning) for x in w)
        r[f"cos_stored_{arm}"] = float(v @ stored)
        r[f"auroc_{arm}"] = roc_auc_score(y, X @ v)
        out = Path(f"cavs_{arm}") / name
        out.mkdir(parents=True, exist_ok=True)
        np.save(out / "concept_v1.npy", v)
    r["cos_conv1_conv4"] = float(vecs["conv1"] @ vecs["conv4"])
    rows.append(r)
    print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k != "cav"}, flush=True)
df = pd.DataFrame(rows)
df.to_csv("refit_summary.tsv", sep="\t", index=False)
print(df.describe().T.to_string())
