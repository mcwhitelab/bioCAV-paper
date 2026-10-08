#!/usr/bin/env python3
"""No-PCA arms: same pos/neg, shared scaler frame (768-d), no PCA.
  np_orig  : C=1, tol=1e-4, max_iter=1000  (what the current scaler-only pipeline's old default would give)
  np_conv1 : C=1, tol=1e-10, max_iter=1e5
  np_conv4 : C=1e-4, tol=1e-10, max_iter=1e5   (= protein-library estimator)
Also stores the stored PCA CAV mapped into the scaler frame (pca.components_.T @ v) as 'stored_sc'
so cosines are comparable across frames."""
import warnings
from pathlib import Path
import numpy as np, pandas as pd, joblib
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning
LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
ARMS = {"np_orig": dict(C=1.0, tol=1e-4, max_iter=1000),
        "np_conv1": dict(C=1.0, tol=1e-10, max_iter=100000),
        "np_conv4": dict(C=1e-4, tol=1e-10, max_iter=100000)}
rows = []
for d in sorted((LIB / "cavs").iterdir()):
    name = d.name
    scaler, pca = joblib.load(d / "scaler_v1.pkl"), joblib.load(d / "pca_v1.pkl")
    pos, neg = np.load(d / "pos.npy"), np.load(d / "neg.npy")
    X = scaler.transform(np.vstack([pos, neg])); y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    st = pca.components_.T @ np.load(d / "concept_v1.npy").astype(float); st /= np.linalg.norm(st)
    c4 = np.load(f"cavs_conv4/{name}/concept_v1.npy"); c4 = pca.components_.T @ c4; c4 /= np.linalg.norm(c4)
    md = X[y == 1].mean(0) - X[y == 0].mean(0); md /= np.linalg.norm(md)
    r = {"cav": name, "n_pos": len(pos), "auroc_stored": roc_auc_score(y, X @ st)}
    vecs = {}
    for arm, kw in ARMS.items():
        clf = LogisticRegression(solver="lbfgs", class_weight="balanced", random_state=42, **kw)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", ConvergenceWarning); clf.fit(X, y)
        v = clf.coef_[0] / np.linalg.norm(clf.coef_[0]); vecs[arm] = v
        r[f"n_iter_{arm}"] = int(clf.n_iter_[0]); r[f"warn_{arm}"] = any(issubclass(x.category, ConvergenceWarning) for x in w)
        r[f"cos_stored_{arm}"] = float(v @ st); r[f"cos_pcaconv4_{arm}"] = float(v @ c4)
        r[f"cos_meandiff_{arm}"] = float(v @ md); r[f"auroc_{arm}"] = roc_auc_score(y, X @ v)
        out = Path(f"cavs_{arm}") / name; out.mkdir(parents=True, exist_ok=True); np.save(out / "concept_v1.npy", v)
    r["cos_np_conv1_np_conv4"] = float(vecs["np_conv1"] @ vecs["np_conv4"])
    rows.append(r); print(name, r["n_iter_np_conv1"], round(r["cos_stored_np_conv4"], 3), flush=True)
df = pd.DataFrame(rows); df.to_csv("refit_nopca_summary.tsv", sep="\t", index=False)
print(df.drop(columns="cav").describe(percentiles=[.05, .5, .95]).T.round(4).to_string())
