#!/usr/bin/env python3
"""Case-selection preview: for every usable pair (static/paired_donor_pairs.csv), how well does
its L2 axis separate normal from cancer, and how does CAV gene r relate to pseudobulk DE?
Paired-donor cells, current (stored) CAVs, PCA-128 frame. Per pair:
  AUROC all cells / excluding the CAV training positives / per donor; mean shift (SD);
  null 95th pct of random L0/L1-orthogonal directions (best orientation);
  ceiling = leave-one-donor-out logistic regression on normal vs cancer embeddings;
  Spearman(CAV r, DESeq2 log2fc_mle) over the DE-tested genes; # DE genes (padj<0.05, |LFC|>=1).
Writes figures/figure_data_staging/singlecell_cases/{case_metrics.tsv, case_cells.csv, case_de_vs_cav.csv}"""
import os, sys
import numpy as np, pandas as pd, joblib, anndata as ad
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, GroupKFold
sys.path.insert(0, os.path.join(os.environ["BIOCAV_REPO"], "core"))
from src.utils.data_loader import load_sequence_embeddings
LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
OUT = Path("../../figures/figure_data_staging/singlecell_cases"); OUT.mkdir(parents=True, exist_ok=True)
DE = Path("de_pseudobulk_allpairs/min10")
pairs = pd.read_csv("../static/paired_donor_pairs.csv").fillna("")
ax = np.load("axes_stored.npz")
emb, ids = load_sequence_embeddings(str(LIB / "embeddings" / "cells.pkl")); ids = np.array([str(i) for i in ids])
row = {c: i for i, c in enumerate(ids)}
norm = lambda s: s.astype(str).str.lower().str.replace("_", " ").str.strip()
a = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r")
obs = a.obs[["cell_type", "tissue", "disease", "assay", "donor_id"]].copy(); obs.index = obs.index.astype(str); a.file.close()
ref = LIB / "cavs" / "epithelial_cell__lung__normal"
sc, pca = joblib.load(ref / "scaler_v1.pkl"), joblib.load(ref / "pca_v1.pkl")
rng = np.random.default_rng(0)
M, CELLS, DEV = [], [], []
for p in pairs.itertuples():
    ct, t, cond = norm(pd.Series([p.cell_type]))[0], norm(pd.Series([p.tissue_real]))[0], norm(pd.Series([p.condition]))[0]
    m = (norm(obs.cell_type) == ct) & (norm(obs.tissue) == t) & norm(obs.disease).isin(["normal", cond])
    if p.assay_filter: m &= obs.assay == p.assay_filter
    sub = obs[m]; dc = sub.groupby("donor_id", observed=True).disease.nunique()
    sub = sub[sub.donor_id.isin(dc[dc == 2].index)]
    E = emb[[row[c] for c in sub.index]]; X = pca.transform(sc.transform(E))
    y = (norm(sub.disease) == cond).astype(int).values
    _, ctk, tk, condk = p.score_col.split("__", 3)
    L0, L1, L2 = ax[f"L0__{ctk}"], ax[f"L1__{ctk}__{tk}"], ax[p.score_col]
    s = X @ L2
    def trained(name):
        P = np.load(LIB / "cavs" / name / "pos.npy").astype(E.dtype); key = {r.tobytes() for r in P}
        return np.array([r.tobytes() in key for r in E])
    tr = trained(f"{ctk}__{tk}__{condk}") | trained(f"{ctk}__{tk}__normal")
    held = ~tr
    null = []
    for _ in range(1000):
        v = rng.standard_normal(len(L2)); v -= (v @ L0) * L0; v -= (v @ L1) * L1; v /= np.linalg.norm(v)
        au = roc_auc_score(y, X @ v); null.append(max(au, 1 - au))
    donors = sub.donor_id.astype(str).values
    pd_au = [roc_auc_score(y[donors == d], s[donors == d]) for d in np.unique(donors) if len(set(y[donors == d])) == 2 and min((y[donors == d] == 0).sum(), (y[donors == d] == 1).sum()) >= 10]
    nd = len(np.unique(donors))
    ceil = roc_auc_score(y, cross_val_predict(LogisticRegression(C=1, max_iter=5000, class_weight="balanced"), X, y,
                                              cv=GroupKFold(nd), groups=donors, method="predict_proba")[:, 1])
    de = pd.read_csv(DE / f"{p.pair_label}.tsv", sep="\t")
    gc = pd.read_csv(LIB / "results" / "gene_correlation_paired_donors" / f"{p.pair_label}.tsv", sep="\t")
    mg = de[["gene", "gene_name", "log2fc_mle", "padj"]].merge(gc[["gene", "r"]], on="gene").dropna(subset=["log2fc_mle", "r"])
    mg["pair"] = p.pair_label; DEV.append(mg)
    M.append({"pair": p.pair_label, "n_donors": nd, "n_normal": int((y == 0).sum()), "n_cancer": int(y.sum()),
              "auroc_all": roc_auc_score(y, s), "auroc_heldout": roc_auc_score(y[held], s[held]) if len(set(y[held])) == 2 else np.nan,
              "n_heldout": int(held.sum()), "shift_sd": (s[y == 1].mean() - s[y == 0].mean()) / s.std(),
              "donor_auroc_min": min(pd_au) if pd_au else np.nan, "donor_auroc_max": max(pd_au) if pd_au else np.nan,
              "null95": np.quantile(null, .95), "ceiling_lodo": ceil,
              "spearman_cav_de": spearmanr(mg.r, mg.log2fc_mle)[0], "n_de_tested": len(mg),
              "n_de_sig": int(((mg.padj < 0.05) & (mg.log2fc_mle.abs() >= 1)).sum())})
    CELLS.append(pd.DataFrame({"pair": p.pair_label, "l2_score": s, "condition": np.where(y == 1, "cancer", "normal"),
                               "heldout": held, "donor": donors}))
    print(M[-1]["pair"], {k: round(v, 3) if isinstance(v, float) else v for k, v in M[-1].items() if k != "pair"}, flush=True)
R = pd.DataFrame(M).sort_values("auroc_heldout", ascending=False)
R.to_csv(OUT / "case_metrics.tsv", sep="\t", index=False)
pd.concat(CELLS).to_csv(OUT / "case_cells.csv", index=False); pd.concat(DEV).to_csv(OUT / "case_de_vs_cav.csv", index=False)
pd.set_option("display.width", 250); print(R.round(3).to_string(index=False))
