#!/usr/bin/env python3
"""
paired_donor_gene_corr.py — Recompute CAV gene correlation restricted to
donors that have BOTH conditions present (a genuine paired normal-vs-cancer
comparison), instead of the pooled correlation across all cells regardless
of donor design.

The pooled correlation cav_gene_correlation.py computes treats every cell as
independent and is vulnerable to Simpson's-paradox-style confounds when a
pair mixes single-condition donors (whose L2-score variation reflects
something other than cancer status) with dual-condition donors -- see the
MT-CO1/lung10x investigation this script grew out of. Restricting to paired
donors removes that confound at the design level rather than statistically
adjusting for it.

For "lung10x" (not a real tissue value -- a naming convention for the
assay-10x-3'v2-restricted axes), pass --assay-filter to additionally
restrict to that assay.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

_biocav_repo_added = False


def norm(v):
    return str(v).lower().replace("_", " ").strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--coords", required=True)
    ap.add_argument("--lib-dir", required=True)
    ap.add_argument("--pairs-csv", required=True,
                    help="CSV with cell_type,tissue,condition,score_col,tissue_real,assay_filter columns")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    import os
    biocav_repo = os.environ.get("BIOCAV_REPO")
    sys.path.insert(0, os.path.join(biocav_repo, "specific_scripts"))
    from cav_gene_correlation import correlate_scores_with_genes  # noqa: E402

    pairs = pd.read_csv(args.pairs_csv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading h5ad...", file=sys.stderr)
    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs
    gene_name_map = dict(zip(adata.var_names, adata.var["feature_name"].astype(str)))

    coords = pd.read_csv(args.coords, sep="\t", index_col=0)
    coords.index = coords.index.astype(str)

    for _, row in pairs.iterrows():
        pair_label = row["pair_label"]
        ct, tissue_real, cond = norm(row["cell_type"]), norm(row["tissue_real"]), norm(row["condition"])
        mask = ((obs["cell_type"].astype(str).map(norm) == ct) &
                (obs["tissue"].astype(str).map(norm) == tissue_real) &
                (obs["disease"].astype(str).map(norm).isin(["normal", cond])))
        if pd.notna(row.get("assay_filter")):
            mask = mask & (obs["assay"] == row["assay_filter"])

        sub_obs = obs[mask].copy()
        sub_obs["disease"] = sub_obs["disease"].astype(str)
        sub_obs["donor_id"] = sub_obs["donor_id"].astype(str)

        dc = sub_obs.groupby("donor_id")["disease"].nunique()
        paired_donors = dc[dc == 2].index
        sub_obs = sub_obs[sub_obs["donor_id"].isin(paired_donors)]

        n_normal = (sub_obs["disease"].map(norm) == "normal").sum()
        n_cond = (sub_obs["disease"].map(norm) == cond).sum()
        print(f"{pair_label}: {len(paired_donors)} paired donors, "
              f"{n_normal} normal + {n_cond} condition = {len(sub_obs)} cells", file=sys.stderr)

        if len(sub_obs) < 10:
            print(f"  WARNING: too few cells — skipping", file=sys.stderr)
            continue

        cell_ids = sub_obs.index.astype(str)
        shared = cell_ids.intersection(coords.index)
        scores = coords.loc[shared, row["score_col"]].values

        adata_sub = adata[shared, :]
        X = adata_sub.X
        if sp.issparse(X):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float64)
        lib_size = X.sum(axis=1, keepdims=True)
        lib_size[lib_size == 0] = 1.0
        X_norm = np.log1p(X / lib_size * 1e4)

        corr_df = correlate_scores_with_genes(scores, X_norm, list(adata_sub.var_names), min_cells=10)
        corr_df.insert(1, "gene_name", corr_df["gene"].map(gene_name_map).fillna(""))
        out_path = out_dir / f"{pair_label}.tsv"
        corr_df.to_csv(out_path, sep="\t", index=False)
        print(f"  Wrote {out_path} ({len(corr_df)} genes)", file=sys.stderr)


if __name__ == "__main__":
    main()
