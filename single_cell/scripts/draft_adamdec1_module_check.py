#!/usr/bin/env python3
"""
DRAFT / exploratory — checking a claimed fibroblast subpopulation module
(ADAMDEC1, CXCL14, EDNRB, PROCR) against the fibroblast/colorectum CAV data.
Not wired into the main figures.R pipeline. Mirrors export_subpop_supp_data.py's
pattern (paired-donor restriction, DE vs CAV scatter, co-expression matrix,
per-cell module score vs L2 axis).

Outputs (into --figure-data-dir, default ../../figures/figure_data):
  draft_adamdec1_gene_corr.csv    gene1, gene2, r
  draft_adamdec1_scatter.csv      pair, gene, gene_name, log2fc, cav_r
  draft_adamdec1_cell_scores.csv  cell_id, disease, module_score, l2_score
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
PAIR_LABEL = "fibroblast__colorectum__normal_vs_colorectal_cancer"
FIGURE_DATA_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../../figures/figure_data")

MODULE_GENES = {
    "ADAMDEC1": "ENSG00000134028",
    "CXCL14": "ENSG00000145824",
    "EDNRB": "ENSG00000136160",
    "PROCR": "ENSG00000101000",
}


def norm(v):
    return str(v).lower().replace("_", " ").strip()


def main():
    out_dir = FIGURE_DATA_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    coords = pd.read_csv(LIB / "results" / "hierarchy" / "cell_coordinates.tsv", sep="\t", index_col=0)
    coords.index = coords.index.astype(str)

    adata = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r")
    obs = adata.obs[["cell_type", "tissue", "disease", "donor_id"]].copy()
    obs.index = obs.index.astype(str)
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["disease"] = obs["disease"].astype(str)

    mask = ((obs["cell_type"].map(norm) == "fibroblast") &
            (obs["tissue"].map(norm) == "colorectum") &
            (obs["disease"].map(norm).isin(["normal", "colorectal cancer"])))
    sub = obs[mask].copy()
    dc = sub.groupby("donor_id")["disease"].nunique()
    paired = dc[dc == 2].index
    sub = sub[sub["donor_id"].isin(paired)]
    shared = sub.index.intersection(coords.index)
    sub = sub.loc[shared]
    scores = coords.loc[shared, "L2__fibroblast__colorectum__colorectal_cancer"]

    adata_sub = adata[shared, :]
    X_full = adata_sub.X
    if sp.issparse(X_full):
        X_full = X_full.toarray()
    X_full = np.asarray(X_full, dtype=np.float64)
    lib_size = X_full.sum(axis=1)

    var_names = list(adata.var_names)
    expr = {}
    for name, gid in MODULE_GENES.items():
        gidx = var_names.index(gid)
        expr[name] = np.log1p(X_full[:, gidx] / np.maximum(lib_size, 1) * 1e4)
    expr_df = pd.DataFrame(expr, index=shared)

    # --- DE vs CAV scatter table ---
    de = pd.read_csv(LIB / "results" / "de_mixedlm" / f"{PAIR_LABEL}.tsv", sep="\t")
    cav = pd.read_csv(LIB / "results" / "gene_correlation_paired_donors" / f"{PAIR_LABEL}.tsv", sep="\t")
    merged = de[["gene", "gene_name", "log2fc"]].merge(cav[["gene", "r"]], on="gene", how="inner")
    highlight = merged[merged["gene_name"].isin(MODULE_GENES)].copy()
    highlight["pair"] = PAIR_LABEL
    highlight[["pair", "gene", "gene_name", "log2fc", "r"]].rename(
        columns={"r": "cav_r"}).to_csv(out_dir / "draft_adamdec1_scatter.csv", index=False)

    # --- co-expression matrix ---
    corr = expr_df.corr()
    corr_long = corr.stack().reset_index()
    corr_long.columns = ["gene1", "gene2", "r"]
    corr_long.to_csv(out_dir / "draft_adamdec1_gene_corr.csv", index=False)

    # --- per-cell module score (z-mean across all 4 genes) vs L2 ---
    z = (expr_df - expr_df.mean()) / expr_df.std()
    module_score = z.mean(axis=1)
    cell_df = pd.DataFrame({
        "cell_id": shared, "disease": sub["disease"].values,
        "module_score": module_score.values, "l2_score": scores.values,
    })
    for name in MODULE_GENES:
        cell_df[f"expr_{name}"] = expr_df[name].values
    cell_df.to_csv(out_dir / "draft_adamdec1_cell_scores.csv", index=False)

    print("--- DE vs CAV for the 4 candidate module genes ---", file=sys.stderr)
    print(highlight[["gene_name", "log2fc", "r"]].to_string(index=False), file=sys.stderr)
    print("\n--- pairwise co-expression (single cell, log1p CPM) ---", file=sys.stderr)
    print(corr.round(2).to_string(), file=sys.stderr)
    print(f"\nWrote draft_adamdec1_{{scatter,gene_corr,cell_scores}}.csv to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
