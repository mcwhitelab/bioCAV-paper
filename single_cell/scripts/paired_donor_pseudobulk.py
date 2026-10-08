#!/usr/bin/env python3
"""
paired_donor_pseudobulk.py — Build donor x condition pseudobulk count matrices
for Figure 6's three case-study pairs, as input to the conventional pseudobulk
DE (scripts/paired_donor_pseudobulk_de.R: DESeq2 ~ donor + condition, apeglm
shrinkage; edgeR QL as a robustness check).

This is the second, conventional-log2FC version of the Figure 6 DE axis. The
first version (paired_donor_de_mixedlm.py) is a per-cell mixed model whose
"log2fc" is a difference in mean log-CP10K, not a log ratio.

Cell selection is identical to paired_donor_de_mixedlm.py /
export_fig5_case_studies.py: same mask, same lung 10x assay filter, same
restriction to donors that contribute BOTH conditions.

Raw integer counts (adata.X, verified integer) are SUMMED per donor x condition.
Note: this deliberately does not reuse de_methods.run_deseq2_pair, which
aggregates per donor only and labels each aggregate with its first cell's
condition -- with paired donors that merges each donor's normal and cancer
cells into one sample.

No min-cells filter is applied here; every sample is written with its n_cells
so the R step can apply (and sensitivity-test) the threshold.

Usage: python scripts/paired_donor_pseudobulk.py   (run from single_cell/)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
OUT_DIR = LIB / "results" / "de_pseudobulk_paired_donors" / "pseudobulk_counts"

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
    print("Loading h5ad...", file=sys.stderr)
    adata = ad.read_h5ad(LIB / "data" / "cells.h5ad")
    obs = adata.obs
    genes = pd.DataFrame({"gene": adata.var_names,
                          "gene_name": adata.var["feature_name"].astype(str).values})

    for pair_label, ct, tissue_real, cond, assay_filter in CASES:
        mask = ((obs["cell_type"].astype(str).map(norm) == norm(ct)) &
                (obs["tissue"].astype(str).map(norm) == norm(tissue_real)) &
                (obs["disease"].astype(str).map(norm).isin(["normal", norm(cond)])))
        if assay_filter:
            mask = mask & (obs["assay"] == assay_filter)

        sub_obs = obs[mask].copy()
        sub_obs["donor_id"] = sub_obs["donor_id"].astype(str)
        sub_obs["disease"] = sub_obs["disease"].astype(str)
        dc = sub_obs.groupby("donor_id")["disease"].nunique()
        paired = set(dc[dc == 2].index)
        keep = (mask & obs["donor_id"].astype(str).isin(paired)).values

        X = adata.X[keep]
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        o = obs[keep]
        donor = o["donor_id"].astype(str).values
        disease = o["disease"].astype(str).values
        condition = np.where(pd.Series(disease).map(norm).values == "normal", "normal", "cancer")

        rows, meta = [], []
        for d in sorted(set(donor)):
            for c in ("normal", "cancer"):
                sel = (donor == d) & (condition == c)
                rows.append(np.asarray(X[sel].sum(axis=0)).ravel())
                meta.append({"sample": f"{d}__{c}", "donor": d, "condition": c,
                             "n_cells": int(sel.sum())})
        counts = np.vstack(rows)
        assert np.allclose(counts, np.round(counts)), "non-integer counts"
        counts_df = pd.DataFrame(counts.astype(np.int64).T, index=adata.var_names,
                                 columns=[m["sample"] for m in meta])
        meta_df = pd.DataFrame(meta)

        counts_df.to_csv(OUT_DIR / f"{pair_label}.counts.tsv.gz", sep="\t")
        meta_df.to_csv(OUT_DIR / f"{pair_label}.samples.tsv", sep="\t", index=False)
        print(f"{pair_label}: {len(paired)}/{len(dc)} donors paired, "
              f"{int(keep.sum())} cells, {len(meta_df)} pseudobulk samples", file=sys.stderr)
        print(meta_df.to_string(index=False), file=sys.stderr)

    genes.to_csv(OUT_DIR / "genes.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
