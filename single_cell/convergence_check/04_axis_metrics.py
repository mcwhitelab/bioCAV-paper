#!/usr/bin/env python3
"""Axis-level comparison for the 3 Fig 6 case pairs: cos of each arm's L2 axis
to the stored-library L2 axis, per-cell L2 score r on paired-donor cells,
in-sample normal-vs-cancer AUROC, and cos(raw CAV, balanced mean difference)."""
import numpy as np, pandas as pd, joblib, anndata as ad
from pathlib import Path
from sklearn.metrics import roc_auc_score
LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
CASES = [("neutrophil", "breast", "breast cancer", None, "L2__neutrophil__breast__breast_cancer"),
         ("epithelial cell", "lung", "lung cancer", "10x 3' v2", "L2__epithelial_cell__lung10x__lung_cancer"),
         ("fibroblast", "colorectum", "colorectal cancer", None, "L2__fibroblast__colorectum__colorectal_cancer")]
ARMS = ["stored", "orig", "conv1", "conv4", "np_orig", "np_conv1", "np_conv4"]
norm = lambda s: s.astype(str).str.lower().str.replace("_", " ").str.strip()
a = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r"); obs = a.obs[["cell_type", "tissue", "disease", "assay", "donor_id"]].copy(); a.file.close()
obs.index = obs.index.astype(str)
# all axes in the scaler frame (PCA arms mapped via pca.components_.T) so cosines compare across frames
axes = {arm: np.load(f"axes_{arm}.npz" if arm.startswith("np_") else f"axes_{arm}_sc.npz") for arm in ARMS}
coords = {arm: pd.read_csv(f"coords_{arm}.tsv", sep="\t", index_col=0, usecols=["cell_id"] + [c[4] for c in CASES]) for arm in ARMS}
for c in coords.values(): c.index = c.index.astype(str)
rows = []
for ct, t, cond, assay, col in CASES:
    m = (norm(obs.cell_type) == ct) & (norm(obs.tissue) == t) & norm(obs.disease).isin(["normal", cond])
    if assay: m &= obs.assay == assay
    sub = obs[m]; dc = sub.groupby("donor_id", observed=True).disease.nunique()
    sub = sub[sub.donor_id.isin(dc[dc == 2].index)]
    y = (norm(sub.disease) == cond).astype(int).values
    s0 = coords["stored"].loc[sub.index, col].values
    for arm in ARMS:
        s = coords[arm].loc[sub.index, col].values
        rows.append({"axis": col, "arm": arm, "n_cells": len(sub), "n_donors": sub.donor_id.nunique(),
                     "cos_L2_vs_stored": float(axes[arm][col] @ axes["stored"][col]),
                     "cos_L2_vs_conv4": float(axes[arm][col] @ axes["conv4"][col]),
                     "r_cellscore_vs_stored": np.corrcoef(s, s0)[0, 1],
                     "auroc_L2_paired_insample": roc_auc_score(y, s)})
df = pd.DataFrame(rows); df.to_csv("axis_metrics.tsv", sep="\t", index=False); print(df.round(4).to_string())

# cos(raw CAV, balanced class-mean difference) in the PCA frame
out = []
for name in ["neutrophil__breast__breast_cancer", "neutrophil__breast__normal", "epithelial_cell__lung10x__lung_cancer",
             "epithelial_cell__lung10x__normal", "fibroblast__colorectum__colorectal_cancer", "fibroblast__colorectum__normal"]:
    d = LIB / "cavs" / name
    sc, pca = joblib.load(d / "scaler_v1.pkl"), joblib.load(d / "pca_v1.pkl")
    P = pca.transform(sc.transform(np.load(d / "pos.npy"))); N = pca.transform(sc.transform(np.load(d / "neg.npy")))
    md = P.mean(0) - N.mean(0); md /= np.linalg.norm(md)
    r = {"cav": name}
    for arm in ["stored", "orig", "conv1", "conv4"]:
        v = np.load(f"cavs_{arm}/{name}/concept_v1.npy").astype(float); r[arm] = float(v @ md / np.linalg.norm(v))
    out.append(r)
print("\ncos(raw CAV, mean difference)"); print(pd.DataFrame(out).round(4).to_string())
pd.DataFrame(out).to_csv("cos_meandiff.tsv", sep="\t", index=False)
