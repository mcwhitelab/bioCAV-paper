#!/usr/bin/env python3
"""
paired_donor_label_corr.py -- per-gene "label r": the Pearson correlation,
across a pair's single cells, between the binary condition label
(0 = normal, 1 = cancer) and log2(CP10K + 1) expression. This is the
CAV gene correlation (paired_donor_gene_corr.py) with the L2 CAV score
replaced by the condition label, i.e. what a per-cell correlation would
give if the CAV axis were nothing more than the normal/cancer call.

Same cells as every other Figure 6 statistic (paired donors only, lung
10x 3' v2 only), same normalisation, same gene set as the CAV table.

Output: results/label_correlation_paired_donors/<pair>.tsv
        (gene, gene_name, label_r, cav_r)

Usage: python scripts/paired_donor_label_corr.py   (run from single_cell/)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
CAV_DIR = LIB / "results" / "gene_correlation_paired_donors"
OUT_DIR = LIB / "results" / "label_correlation_paired_donors"

# pair_label, cell_type, tissue_real, condition, assay_filter
CASES = [
    ("neutrophil__breast__normal_vs_breast_cancer", "neutrophil", "breast", "breast cancer", None),
    ("epithelial_cell__lung10x__normal_vs_lung_cancer", "epithelial cell", "lung", "lung cancer", "10x 3' v2"),
    ("fibroblast__colorectum__normal_vs_colorectal_cancer", "fibroblast", "colorectum", "colorectal cancer", None),
]


def norm(v):
    return str(v).lower().replace("_", " ").strip()


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(LIB / "data" / "cells.h5ad")
    obs = adata.obs
    col = pd.Series(np.arange(adata.n_vars), index=adata.var_names)

    for pair_label, ct, tissue_real, cond, assay_filter in CASES:
        mask = ((obs["cell_type"].astype(str).map(norm) == norm(ct)) &
                (obs["tissue"].astype(str).map(norm) == norm(tissue_real)) &
                (obs["disease"].astype(str).map(norm).isin(["normal", norm(cond)])))
        if assay_filter:
            mask = mask & (obs["assay"] == assay_filter)
        sub = obs[mask]
        dc = sub.groupby(sub["donor_id"].astype(str))["disease"].nunique()
        mask = mask & obs["donor_id"].astype(str).isin(dc[dc == 2].index)

        X = adata.X[mask.values]
        X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
        lib = np.asarray(X.sum(axis=1)).ravel()          # all genes, before filtering
        y = (obs.loc[mask, "disease"].astype(str).map(norm) != "normal").values.astype(float)

        cav = pd.read_csv(CAV_DIR / f"{pair_label}.tsv", sep="\t")
        cav = cav[cav["gene"].isin(col.index)]
        G = X[:, col[cav["gene"]].values].toarray()
        G = np.log2(G / lib[:, None] * 1e4 + 1)

        yc = y - y.mean()
        Gc = G - G.mean(axis=0)
        sd = Gc.std(axis=0)
        label_r = np.where(sd > 0, (yc @ Gc) / (len(y) * y.std() * np.where(sd > 0, sd, 1)), np.nan)

        out = pd.DataFrame({"gene": cav["gene"].values, "gene_name": cav["gene_name"].values,
                            "label_r": label_r, "cav_r": cav["r"].values}).dropna()
        out.to_csv(OUT_DIR / f"{pair_label}.tsv", sep="\t", index=False)
        rho = out["label_r"].corr(out["cav_r"], method="spearman")
        print(f"{pair_label}: {int(mask.sum())} cells, {len(out)} genes, "
              f"Spearman(label r, CAV r) = {rho:.2f}", file=sys.stderr)


if __name__ == "__main__":
    main()
