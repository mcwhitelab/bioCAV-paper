#!/usr/bin/env python3
"""
lung10x_de_and_corr.py — DE (mixedlm) + gene correlation for the new
epithelial_cell__lung10x__normal_vs_lung_cancer pair, bypassing the usual
cell_type/tissue/disease metadata matching (there's no real "lung10x" tissue
value in obs -- it's a naming convention for the assay-restricted axes, not
a real category) by selecting cells explicitly (assay == "10x 3' v2") and
calling the same underlying functions the normal pipeline uses directly.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

_biocav_repo = os.environ.get("BIOCAV_REPO")
sys.path.insert(0, os.path.join(_biocav_repo, "specific_scripts"))
from de_pseudobulk import run_mixedlm_pair            # noqa: E402
from cav_gene_correlation import correlate_scores_with_genes  # noqa: E402

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
PAIR_LABEL = "epithelial_cell__lung10x__normal_vs_lung_cancer"


def norm(v):
    return str(v).lower().replace("_", " ").strip()


def main():
    print("Loading h5ad...", file=sys.stderr)
    adata = ad.read_h5ad(LIB / "data" / "cells.h5ad")
    obs = adata.obs

    mask = ((obs["cell_type"].map(norm) == "epithelial cell") &
            (obs["tissue"].map(norm) == "lung") &
            (obs["disease"].map(norm).isin(["normal", "lung cancer"])) &
            (obs["assay"] == "10x 3' v2"))
    print(f"Selected {mask.sum()} cells (10x 3' v2, epithelial, lung, normal/lung_cancer)", file=sys.stderr)
    adata_pair = adata[mask].copy()

    gene_name_map = dict(zip(adata.var_names, adata.var["feature_name"].astype(str)))

    # --- DE (mixedlm) ---
    print("Running mixedlm DE...", file=sys.stderr)
    de_df = run_mixedlm_pair(adata_pair, donor_col="donor_id", condition_col="disease",
                              baseline_value="normal")
    if de_df is None or de_df.empty:
        print("ERROR: DE returned no results", file=sys.stderr)
        sys.exit(1)
    de_df.insert(1, "gene_name", de_df["gene"].map(gene_name_map).fillna(""))
    de_df = (de_df.assign(abs_lfc=de_df["log2fc"].abs())
             .sort_values("abs_lfc", ascending=False).drop(columns="abs_lfc")
             .reset_index(drop=True))
    de_out = LIB / "results" / "de_mixedlm" / f"{PAIR_LABEL}.tsv"
    de_out.parent.mkdir(parents=True, exist_ok=True)
    de_df.to_csv(de_out, sep="\t", index=False)
    print(f"Wrote {de_out} ({len(de_df)} genes)", file=sys.stderr)

    # --- Gene correlation vs. L2 score ---
    print("Computing gene correlation vs. L2__epithelial_cell__lung10x__lung_cancer...", file=sys.stderr)
    coords = pd.read_csv(LIB / "results" / "hierarchy_celltype_tissue" / "cell_coordinates.tsv",
                         sep="\t", index_col=0)
    coords.index = coords.index.astype(str)
    cell_ids = adata_pair.obs_names.astype(str)
    scores = coords.reindex(cell_ids)["L2__epithelial_cell__lung10x__lung_cancer"].values

    import scipy.sparse as sp
    X = adata_pair.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    lib_size = X.sum(axis=1, keepdims=True)
    lib_size[lib_size == 0] = 1.0
    X_norm = np.log1p(X / lib_size * 1e4)

    corr_df = correlate_scores_with_genes(scores, X_norm, list(adata_pair.var_names), min_cells=10)
    corr_df.insert(1, "gene_name", corr_df["gene"].map(gene_name_map).fillna(""))
    corr_out = LIB / "results" / "gene_correlation_celltype_tissue" / f"{PAIR_LABEL}.tsv"
    corr_out.parent.mkdir(parents=True, exist_ok=True)
    corr_df.to_csv(corr_out, sep="\t", index=False)
    print(f"Wrote {corr_out} ({len(corr_df)} genes)", file=sys.stderr)


if __name__ == "__main__":
    main()
