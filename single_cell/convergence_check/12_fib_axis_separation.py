#!/usr/bin/env python3
"""How well does the fibroblast/colorectum L2 axis separate normal from cancer cells?
Paired-donor fibroblasts, current (stored) CAVs, PCA-128 frame.
 - AUROC overall, per donor, and excluding the cells used as CAV training positives
   (matched by identical embedding rows of pos.npy)
 - null: AUROC of random unit directions orthogonalised against the same L0/L1 basis"""
import numpy as np, pandas as pd, joblib, anndata as ad, os, sys
from pathlib import Path
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.join(os.environ["BIOCAV_REPO"], "core"))
from src.utils.data_loader import load_sequence_embeddings
LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
ARM = sys.argv[1] if len(sys.argv) > 1 else "stored"   # conv1 = converged C=1 (Fig 6 since 2026-10-07)
ax = np.load(f"axes_{ARM}.npz")
L0, L1, L2 = ax["L0__fibroblast"], ax["L1__fibroblast__colorectum"], ax["L2__fibroblast__colorectum__colorectal_cancer"]
d = LIB / "cavs" / "fibroblast__colorectum__colorectal_cancer"
sc, pca = joblib.load(d / "scaler_v1.pkl"), joblib.load(d / "pca_v1.pkl")
emb, ids = load_sequence_embeddings(str(LIB / "embeddings" / "cells.pkl")); ids = np.array([str(i) for i in ids])
norm = lambda s: s.astype(str).str.lower().str.replace("_", " ")
a = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r")
obs = a.obs[["cell_type", "tissue", "disease", "donor_id"]].copy(); obs.index = obs.index.astype(str); a.file.close()
m = (norm(obs.cell_type) == "fibroblast") & (norm(obs.tissue) == "colorectum") & norm(obs.disease).isin(["normal", "colorectal cancer"])
sub = obs[m]; dc = sub.groupby("donor_id", observed=True).disease.nunique()
sub = sub[sub.donor_id.isin(dc[dc == 2].index)]
row = {c: i for i, c in enumerate(ids)}
E = emb[[row[c] for c in sub.index]]
X = pca.transform(sc.transform(E))
y = (norm(sub.disease) == "colorectal cancer").astype(int).values
s = X @ L2
# training positives: exact embedding-row matches to pos.npy (cancer CAV) and normal CAV positives
def trained(path):
    P = np.load(path); key = {r.tobytes() for r in P.astype(E.dtype)}
    return np.array([r.tobytes() in key for r in E])
tr_c = trained(d / "pos.npy"); tr_n = trained(LIB / "cavs" / "fibroblast__colorectum__normal" / "pos.npy")
print(f"paired-donor cells: {len(y)} ({y.sum()} cancer); used as cancer-CAV positives: {tr_c.sum()}, normal-CAV positives: {tr_n.sum()}")
print(f"AUROC all cells: {roc_auc_score(y, s):.3f}")
held = ~(tr_c | tr_n)
print(f"AUROC excluding all CAV training cells (n={held.sum()}, {y[held].sum()} cancer): {roc_auc_score(y[held], s[held]):.3f}")
print("mean L2 normal %.2f (sd %.2f), cancer %.2f (sd %.2f) -> shift = %.2f sd" % (
    s[y == 0].mean(), s[y == 0].std(), s[y == 1].mean(), s[y == 1].std(), (s[y == 1].mean() - s[y == 0].mean()) / s.std()))
print("per donor AUROC:", {dn: round(roc_auc_score(y[g], s[g]), 3) for dn in sub.donor_id.unique()
                           for g in [(sub.donor_id == dn).values] if len(set(y[g])) == 2})
rng = np.random.default_rng(0); null = []
for _ in range(2000):
    v = rng.standard_normal(len(L2)); v -= (v @ L0) * L0; v -= (v @ L1) * L1; v /= np.linalg.norm(v)
    au = roc_auc_score(y, X @ v); null.append(max(au, 1 - au))   # best orientation, to be fair to the null
null = np.array(null)
print(f"null (random L0/L1-orthogonal directions, best orientation): median {np.median(null):.3f}, 95th pct {np.quantile(null, .95):.3f}, 99th {np.quantile(null, .99):.3f}; "
      f"fraction >= observed: {(null >= roc_auc_score(y, s)).mean():.4f}")
# what the normal-vs-cancer separating direction in this space could reach (in-sample LDA-ish upper bound)
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, GroupKFold
pr = cross_val_predict(LogisticRegression(C=1, max_iter=5000, class_weight="balanced"), X, y, cv=GroupKFold(7), groups=sub.donor_id.values, method="predict_proba")[:, 1]
print(f"reference: a classifier trained directly on normal vs cancer fibroblasts, leave-one-donor-out AUROC = {roc_auc_score(y, pr):.3f}")
