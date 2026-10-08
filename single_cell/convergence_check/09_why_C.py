#!/usr/bin/env python3
"""Why C=1e-4 changes the single-cell L2 axes (fibroblast/colorectum focus).
In the PCA-128 frame (stored, conv1 = C=1 converged, conv4 = C=1e-4 converged):
 1. cos(CAV, mean difference) and cos(CAV, ridge-LDA Sigma_w^-1 dmu)
 2. background variance captured by the direction (w' Sigma_bg w, relative to the mean eigenvalue)
 3. cos(normal CAV, cancer CAV) and the norm of the L2 residual before normalisation
 4. lymphocyte-marker+ fraction among fibroblast positives, normal vs cancer"""
import numpy as np, pandas as pd, joblib, anndata as ad, scipy.sparse as sp
from pathlib import Path
LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
H = {a: np.load(f"axes_{a}.npz") for a in ["stored", "conv1", "conv4"]}
rows = []
for ct, t, cond in [("fibroblast", "colorectum", "colorectal_cancer"), ("epithelial_cell", "lung10x", "lung_cancer"),
                    ("neutrophil", "breast", "breast_cancer")]:
    dN, dC = LIB / "cavs" / f"{ct}__{t}__normal", LIB / "cavs" / f"{ct}__{t}__{cond}"
    sc, pca = joblib.load(dC / "scaler_v1.pkl"), joblib.load(dC / "pca_v1.pkl")
    f = lambda A: pca.transform(sc.transform(A))
    PN, PC, N = f(np.load(dN / "pos.npy")), f(np.load(dC / "pos.npy")), f(np.load(dC / "neg.npy"))
    Sbg = np.cov(N.T); ev = np.trace(Sbg) / len(Sbg)
    def lda(P):
        Sw = 0.5 * (np.cov(P.T) + Sbg); v = np.linalg.solve(Sw + 1e-3 * ev * np.eye(len(Sw)), P.mean(0) - N.mean(0)); return v / np.linalg.norm(v)
    md = lambda P: (P.mean(0) - N.mean(0)) / np.linalg.norm(P.mean(0) - N.mean(0))
    # direction of the actual normal->cancer shift among positives
    shift = PC.mean(0) - PN.mean(0); shift /= np.linalg.norm(shift)
    for arm in ["stored", "conv1", "conv4"]:
        u = lambda d: (lambda v: v / np.linalg.norm(v))(np.load(f"cavs_{arm}/{d.name}/concept_v1.npy").astype(float))
        wN, wC = u(dN), u(dC)
        L2 = H[arm][f"L2__{ct}__{t}__{cond}"]
        L0 = H[arm][f"L0__{ct}"]; L1 = H[arm][f"L1__{ct}__{t}"]
        resid = wC - (wC @ L0) * L0 - (wC @ L1) * L1
        rows.append({"pair": f"{ct}/{t}", "arm": arm,
                     "cos(wC, meandiff)": wC @ md(PC), "cos(wC, LDA)": wC @ lda(PC),
                     "bg var of wC (x mean eig)": wC @ Sbg @ wC / ev,
                     "bg var of L2 (x mean eig)": L2 @ Sbg @ L2 / ev,
                     "cos(wN, wC)": wN @ wC, "|L2 residual|": np.linalg.norm(resid),
                     "cos(L2, mean N->C shift)": L2 @ shift})
D = pd.DataFrame(rows); D.to_csv("why_C.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(D.round(3).to_string(index=False))

# lymphocyte contamination among paired-donor fibroblasts, by condition
a = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r")
obs = a.obs[["cell_type", "tissue", "disease", "donor_id"]].copy(); obs.index = obs.index.astype(str)
n = lambda s: s.astype(str).str.lower().str.replace("_", " ")
m = (n(obs.cell_type) == "fibroblast") & (n(obs.tissue) == "colorectum") & n(obs.disease).isin(["normal", "colorectal cancer"])
cells = obs.index[m]
name = pd.Series(a.var["feature_name"].astype(str).values, index=a.var_names.astype(str))
genes = ["PTPRC", "CD3E", "CD52", "MS4A1", "CD79A", "COL1A1", "DCN"]
idx = [list(name.index).index(name[name == g].index[0]) for g in genes]
X = a[cells, :].X; X = X.toarray() if sp.issparse(X) else np.asarray(X); a.file.close()
cnt = pd.DataFrame(X[:, idx] > 0, columns=genes, index=cells)
cnt["disease"] = obs.loc[cells, "disease"].astype(str).values
print("\nfraction of fibroblast cells (all donors) with >0 counts:")
print(cnt.groupby("disease").mean().round(3).to_string())
lym = (cnt[["PTPRC", "CD3E", "CD52"]].sum(1) >= 2)
print("\nPTPRC/CD3E/CD52 >=2 of 3 detected:", lym.groupby(cnt.disease).mean().round(3).to_dict())
