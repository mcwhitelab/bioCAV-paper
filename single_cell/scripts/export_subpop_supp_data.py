#!/usr/bin/env python3
"""
export_subpop_supp_data.py — Export data for a supplemental figure arguing
CAV correlation can detect subpopulation-abundance shifts that population-
average DE (mixedlm) misses, using fibroblast/colorectum's quiescent-
fibroblast module (MGP, DCN, OGN, C3, CCDC80) as the worked example: all
five are DE-null (padj > 0.5) but strong, highly significant CAV
correlates (r=-0.33 to -0.37, p ~1e-53 to 1e-65), and are mutually
co-expressed (pairwise r=0.54-0.79), consistent with marking a real
subpopulation rather than being independently regulated per-cell.

Outputs (into --figure-data-dir):
  subpop_gene_scatter.csv   pair, gene, gene_name, log2fc, cav_r, module (quiescent/other/activated)
  subpop_gene_corr.csv      gene1, gene2, r   (long-format correlation matrix)
  subpop_cell_scores.csv    cell_id, disease, quiescent_score, l2_score
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
PAIR_LABEL = "fibroblast__colorectum__normal_vs_colorectal_cancer"
FIGURE_DATA_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../figures/figure_data")

QUIESCENT_GENES = {"MGP": "ENSG00000111341", "DCN": "ENSG00000011465", "OGN": "ENSG00000106809",
                   "C3": "ENSG00000125730", "CCDC80": "ENSG00000091986"}
OTHER_GENES = {"ADAMDEC1": "ENSG00000134028", "CXCR4": "ENSG00000121966",
               "CXCL14": "ENSG00000145824"}


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
    scores = coords.loc[shared, f"L2__fibroblast__colorectum__colorectal_cancer"]

    all_genes = {**QUIESCENT_GENES, **OTHER_GENES}
    adata_sub = adata[shared, :]
    X_full = adata_sub.X
    if sp.issparse(X_full):
        X_full = X_full.toarray()
    X_full = np.asarray(X_full, dtype=np.float64)
    lib_size = X_full.sum(axis=1)

    var_names = list(adata.var_names)
    expr = {}
    for name, gid in all_genes.items():
        gidx = var_names.index(gid)
        expr[name] = np.log1p(X_full[:, gidx] / np.maximum(lib_size, 1) * 1e4)
    expr_df = pd.DataFrame(expr, index=shared)

    # --- 1. scatter table (DE log2fc vs CAV r) ---
    de = pd.read_csv(LIB / "results" / "de_mixedlm" / f"{PAIR_LABEL}.tsv", sep="\t")
    cav = pd.read_csv(LIB / "results" / "gene_correlation_paired_donors" / f"{PAIR_LABEL}.tsv", sep="\t")
    merged = de[["gene", "gene_name", "log2fc"]].merge(cav[["gene", "r"]], on="gene", how="inner")
    module_map = {**{n: "quiescent" for n in QUIESCENT_GENES}, "CXCR4": "activated", "ADAMDEC1": "independent"}
    highlight = merged[merged["gene_name"].isin(module_map)].copy()
    highlight["module"] = highlight["gene_name"].map(module_map)
    highlight["pair"] = PAIR_LABEL
    highlight[["pair", "gene", "gene_name", "log2fc", "r", "module"]].to_csv(
        out_dir / "subpop_gene_scatter.csv", index=False)

    # background (all genes, for context cloud)
    merged["pair"] = PAIR_LABEL
    merged.rename(columns={"r": "cav_r"}).to_csv(out_dir / "subpop_gene_scatter_background.csv", index=False)

    # --- 2. correlation matrix ---
    corr = expr_df.corr()
    corr_long = corr.stack().reset_index()
    corr_long.columns = ["gene1", "gene2", "r"]
    corr_long.to_csv(out_dir / "subpop_gene_corr.csv", index=False)

    # --- 3. per-cell quiescent score + L2 score + composition flags ---
    z = (expr_df[list(QUIESCENT_GENES)] - expr_df[list(QUIESCENT_GENES)].mean()) / expr_df[list(QUIESCENT_GENES)].std()
    quiescent_score = z.mean(axis=1)
    module_negative = (expr_df[list(QUIESCENT_GENES)] == 0).all(axis=1)
    high_quiescent = quiescent_score > quiescent_score.quantile(0.75)
    cell_df = pd.DataFrame({
        "cell_id": shared, "disease": sub["disease"].values,
        "quiescent_score": quiescent_score.values, "l2_score": scores.values,
        "module_negative": module_negative.values, "high_quiescent": high_quiescent.values,
    })
    cell_df.to_csv(out_dir / "subpop_cell_scores.csv", index=False)

    print(f"Wrote subpop_gene_scatter.csv, subpop_gene_scatter_background.csv, "
          f"subpop_gene_corr.csv, subpop_cell_scores.csv to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
