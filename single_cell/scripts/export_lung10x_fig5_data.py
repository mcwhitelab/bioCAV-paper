#!/usr/bin/env python3
"""
export_lung10x_fig5_data.py — Append the epithelial_cell__lung10x__normal_vs_
lung_cancer pair to Figure 5's sc_*.csv files (produced by 08_export_fig5_data.sh
for the other case-study pairs).

get_pair_data() (from cav_continuum_viz.py) selects cells by matching
cell_type/tissue/condition against real obs columns -- but "lung10x" isn't a
real tissue value, it's a naming convention for the assay-restricted axes.
So rather than special-casing cav_continuum_viz.py (shared infra used
elsewhere), this script patches a COPY of obs["tissue"] to "lung10x" for
exactly the 10x 3' v2 epithelial/lung cells, then calls get_pair_data()
normally -- it takes obs as a plain argument, so this works transparently.

Run AFTER 08_export_fig5_data.sh (with epithelial_cell__colorectum swapped
out of its --pairs list) -- this script reads and appends to its output.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from cav_continuum_viz import load_coords, load_obs, load_expression, get_pair_data, minmax  # noqa: E402

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
PAIR_LABEL = "epithelial_cell__lung10x__normal_vs_lung_cancer"
FIGURE_DATA_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../figures/figure_data")


def norm(v):
    return str(v).lower().replace("_", " ").strip()


def main():
    coords = load_coords(str(LIB / "results" / "hierarchy" / "cell_coordinates.tsv"), "L2")
    obs = load_obs(str(LIB / "data" / "cells.h5ad"), ["cell_type", "tissue", "disease", "assay"])
    X_expr, gene_ids, expr_cell_ids, gene_name_map = load_expression(str(LIB / "data" / "cells.h5ad"))

    lib_size = X_expr.sum(axis=1, keepdims=True)
    lib_size[lib_size == 0] = 1.0
    X_expr = np.log1p(X_expr / lib_size * 1e4)

    mask = ((obs["cell_type"].map(norm) == "epithelial cell") &
            (obs["tissue"].map(norm) == "lung") &
            (obs["disease"].map(norm).isin(["normal", "lung cancer"])) &
            (obs["assay"] == "10x 3' v2"))
    print(f"Patching tissue -> 'lung10x' for {mask.sum()} cells", file=sys.stderr)
    obs_patched = obs.copy()
    if hasattr(obs_patched["tissue"], "cat"):
        obs_patched["tissue"] = obs_patched["tissue"].astype(str)
    obs_patched.loc[mask, "tissue"] = "lung10x"

    pd_ = get_pair_data(
        pair_label=PAIR_LABEL, coords=coords, X_expr=X_expr,
        gene_ids=gene_ids, gene_name_map=gene_name_map, obs=obs_patched,
        group_col="cell_type", context_col="tissue", condition_col="disease",
        level="L2", gene_corr_dir=LIB / "results" / "gene_correlation_celltype_tissue",
        top_n=8, max_cells=1000, baseline_value="normal", top_n_each=5,
    )
    if pd_ is None:
        print("ERROR: get_pair_data returned None", file=sys.stderr)
        sys.exit(1)

    meta_row = {"pair": PAIR_LABEL, "group": pd_["group"], "context": pd_["context"],
                "baseline": pd_["baseline"], "condition": pd_["condition"],
                "n_normal": int(pd_["n_normal"]), "n_cancer": int(pd_["n_cancer"])}

    is_base = pd_["labels"] == pd_["baseline"]
    cell_df = pd.DataFrame({"pair": PAIR_LABEL, "l2_score": pd_["scores"], "is_baseline": is_base})

    pos = pd_["gene_df"][pd_["gene_df"]["r"] >= 0].sort_values("r", ascending=False)
    neg = pd_["gene_df"][pd_["gene_df"]["r"] < 0].sort_values("r", ascending=True)
    gene_df = pd.concat([pos, neg]).reset_index(drop=True)
    gene_rows = []
    for rank, row in gene_df.iterrows():
        gid = row["gene"]
        if gid not in pd_["gene_expr"]:
            continue
        expr_sc = minmax(pd_["gene_expr"][gid].astype(np.float32), clip_pct=99.0)
        gene_rows.append(pd.DataFrame({"pair": PAIR_LABEL, "gene": gid, "gene_name": row["display"],
                                       "r": row["r"], "rank": rank, "l2_score": pd_["scores"],
                                       "expr_scaled": expr_sc}))
    gene_df_out = pd.concat(gene_rows, ignore_index=True)

    de = pd.read_csv(LIB / "results" / "de_mixedlm" / f"{PAIR_LABEL}.tsv", sep="\t")
    cav = pd.read_csv(LIB / "results" / "gene_correlation_celltype_tissue" / f"{PAIR_LABEL}.tsv", sep="\t")
    merged = (de[["gene", "log2fc", "gene_name"]].rename(columns={"log2fc": "lfc"})
              .merge(cav[["gene", "r"]].rename(columns={"r": "cav_r"}), on="gene", how="inner")
              .dropna(subset=["lfc", "cav_r"]))
    merged["pair"] = PAIR_LABEL
    de_df_out = merged[["pair", "gene", "gene_name", "lfc", "cav_r"]].rename(columns={"lfc": "log2fc"})

    # Append to existing CSVs from 08_export_fig5_data.sh
    for name, new_df in [("sc_pair_meta.csv", pd.DataFrame([meta_row])),
                         ("sc_continuum_cells.csv", cell_df),
                         ("sc_continuum_genes.csv", gene_df_out),
                         ("sc_de_vs_cav.csv", de_df_out)]:
        path = FIGURE_DATA_DIR / name
        existing = pd.read_csv(path)
        existing = existing[existing["pair"] != PAIR_LABEL]  # idempotent re-run
        combined = pd.concat([existing, new_df], ignore_index=True)
        combined.to_csv(path, index=False)
        print(f"  Updated {path} (+{len(new_df)} rows for {PAIR_LABEL})", file=sys.stderr)


if __name__ == "__main__":
    main()
