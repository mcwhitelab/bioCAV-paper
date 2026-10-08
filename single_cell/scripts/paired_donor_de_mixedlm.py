#!/usr/bin/env python3
"""
paired_donor_de_mixedlm.py — Re-run the per-gene donor-random-intercept mixed
model DE for Figure 6's three case-study pairs, restricted to donors that
contribute BOTH conditions, so the DE axis of the DE-vs-CAV scatter is computed
on exactly the same cells as the CAV axis (results/gene_correlation_paired_donors/).

For neutrophil/breast (3/3 donors paired) and fibroblast/colorectum (7/7) this
reproduces the existing results/de_mixedlm/ tables; for epithelial/lung10x only
4 of 8 donors are paired, so this is a genuine re-fit on 658 rather than 1407 cells.

Cell selection is the same mask + donor-pairing logic as
export_fig5_case_studies.py, applied before the model sees any obs.

Usage: python paired_donor_de_mixedlm.py   (run from single_cell/)
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

_biocav_repo = os.environ.get("BIOCAV_REPO",
                              "/groups/clairemcwhite/cav_workspace/bioCAV")
sys.path.insert(0, os.path.join(_biocav_repo, "specific_scripts"))
from de_methods import run_mixedlm_pair  # noqa: E402

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
OUT_DIR = LIB / "results" / "de_mixedlm_paired_donors"

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
    gene_name_map = dict(zip(adata.var_names, adata.var["feature_name"].astype(str)))

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
        keep = mask & obs["donor_id"].astype(str).isin(paired)

        adata_pair = adata[keep.values].copy()
        adata_pair.obs["disease"] = adata_pair.obs["disease"].astype(str)
        adata_pair.obs["donor_id"] = adata_pair.obs["donor_id"].astype(str)
        vc = adata_pair.obs["disease"].value_counts()
        print(f"{pair_label}: {len(paired)}/{len(dc)} donors paired, "
              f"{adata_pair.n_obs} cells "
              f"({{{', '.join(f'{k}={v}' for k, v in vc.items() if v > 0)}}})",
              file=sys.stderr)

        de_df = run_mixedlm_pair(adata_pair, donor_col="donor_id",
                                 condition_col="disease", baseline_value="normal")
        if de_df is None or de_df.empty:
            print(f"  ERROR: no DE results for {pair_label}", file=sys.stderr)
            sys.exit(1)

        de_df.insert(1, "gene_name", de_df["gene"].map(gene_name_map).fillna(""))
        de_df = (de_df.assign(abs_lfc=de_df["log2fc"].abs())
                 .sort_values("abs_lfc", ascending=False).drop(columns="abs_lfc")
                 .reset_index(drop=True))
        out = OUT_DIR / f"{pair_label}.tsv"
        de_df.to_csv(out, sep="\t", index=False)
        print(f"  wrote {out} ({len(de_df)} genes, "
              f"{(~de_df['converged']).sum()} non-converged)", file=sys.stderr)


if __name__ == "__main__":
    main()
