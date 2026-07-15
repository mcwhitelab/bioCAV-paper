#!/usr/bin/env python3
"""
export_fig5_case_studies.py — Export sc_*.csv for Figure 5's three
case-study pairs, all restricted to donors with BOTH conditions present (a
genuine paired comparison -- see the MT-CO1/epithelial-lung Simpson's-
paradox investigation this grew out of). Replaces the earlier two-script
split (08_export_fig5_data.sh + export_lung10x_fig5_data.py) with one
consistent path so the continuum strips and the DE-vs-CAV scatter both
reflect the same paired-donor cell population the gene correlations were
computed on.

For "lung10x" (not a real tissue value -- a naming convention for the
assay-10x-3'v2-restricted axes), the tissue column is additionally patched
for the assay-matched subset before donor-pairing is applied.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from cav_continuum_viz import load_coords, load_obs, load_expression, get_pair_data, minmax  # noqa: E402

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
GENE_CORR_DIR = LIB / "results" / "gene_correlation_paired_donors"
DE_DIR = LIB / "results" / "de_mixedlm"
FIGURE_DATA_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../figures/figure_data")

CASES = [
    # pair_label, cell_type, tissue_real, condition, assay_filter (None if not needed)
    ("neutrophil__breast__normal_vs_breast_cancer", "neutrophil", "breast", "breast cancer", None),
    ("epithelial_cell__lung10x__normal_vs_lung_cancer", "epithelial cell", "lung", "lung cancer", "10x 3' v2"),
    # fibroblast/colorectum: 7/7 donors paired (no cells lost), and CAV
    # recovers ADAMDEC1 (r=0.21, p=7e-21) as a top correlate despite near-
    # zero DE log2fc (0.016, p=0.77) -- ADAMDEC1 marks a specific mucosal
    # fibroblast subpopulation in gut atlases (e.g. Kinchen et al.), so this
    # likely reflects a subpopulation-abundance shift that a population-
    # average DE test misses but the per-cell CAV score picks up.
    ("fibroblast__colorectum__normal_vs_colorectal_cancer", "fibroblast", "colorectum", "colorectal cancer", None),
]


def norm(v):
    return str(v).lower().replace("_", " ").strip()


def main():
    coords = load_coords(str(LIB / "results" / "hierarchy" / "cell_coordinates.tsv"), "L2")
    obs = load_obs(str(LIB / "data" / "cells.h5ad"),
                   ["cell_type", "tissue", "disease", "assay", "donor_id"])
    obs["donor_id"] = obs["donor_id"].astype(str)
    X_expr, gene_ids, expr_cell_ids, gene_name_map = load_expression(str(LIB / "data" / "cells.h5ad"))

    lib_size = X_expr.sum(axis=1, keepdims=True)
    lib_size[lib_size == 0] = 1.0
    X_expr = np.log1p(X_expr / lib_size * 1e4)

    if hasattr(obs["tissue"], "cat"):
        obs["tissue"] = obs["tissue"].astype(str)
    if hasattr(obs["disease"], "cat"):
        obs["disease"] = obs["disease"].astype(str)

    meta_rows, cell_rows, gene_rows, de_rows = [], [], [], []

    for pair_label, ct, tissue_real, cond, assay_filter in CASES:
        mask = ((obs["cell_type"].map(norm) == norm(ct)) &
                (obs["tissue"].map(norm) == norm(tissue_real)) &
                (obs["disease"].map(norm).isin(["normal", norm(cond)])))
        if assay_filter:
            mask = mask & (obs["assay"] == assay_filter)

        sub = obs[mask]
        dc = sub.groupby("donor_id")["disease"].nunique()
        paired_donors = set(dc[dc == 2].index)
        excluded = obs.index[mask & ~obs["donor_id"].isin(paired_donors)]
        print(f"{pair_label}: excluding {len(excluded)} cells from non-paired donors "
              f"({len(paired_donors)} paired donors kept)", file=sys.stderr)

        obs_patched = obs.copy()
        tissue_label = "lung10x" if assay_filter else tissue_real
        if assay_filter:
            obs_patched.loc[mask & obs["donor_id"].isin(paired_donors), "tissue"] = tissue_label
        # Exclude non-paired-donor cells from matching by giving them a
        # disease value that won't match "normal" or the real condition.
        obs_patched.loc[excluded, "disease"] = "excluded_unpaired_donor"

        pd_ = get_pair_data(
            pair_label=pair_label, coords=coords, X_expr=X_expr,
            gene_ids=gene_ids, gene_name_map=gene_name_map, obs=obs_patched,
            group_col="cell_type", context_col="tissue", condition_col="disease",
            level="L2", gene_corr_dir=GENE_CORR_DIR, top_n=8, max_cells=1000,
            baseline_value="normal", top_n_each=5,
        )
        if pd_ is None:
            print(f"  ERROR: get_pair_data returned None for {pair_label}", file=sys.stderr)
            continue

        meta_rows.append({"pair": pair_label, "group": pd_["group"], "context": pd_["context"],
                          "baseline": pd_["baseline"], "condition": pd_["condition"],
                          "n_normal": int(pd_["n_normal"]), "n_cancer": int(pd_["n_cancer"])})

        is_base = pd_["labels"] == pd_["baseline"]
        cell_rows.append(pd.DataFrame({"pair": pair_label, "l2_score": pd_["scores"], "is_baseline": is_base}))

        pos = pd_["gene_df"][pd_["gene_df"]["r"] >= 0].sort_values("r", ascending=False)
        neg = pd_["gene_df"][pd_["gene_df"]["r"] < 0].sort_values("r", ascending=True)
        gene_df = pd.concat([pos, neg]).reset_index(drop=True)
        for rank, row in gene_df.iterrows():
            gid = row["gene"]
            if gid not in pd_["gene_expr"]:
                continue
            expr_sc = minmax(pd_["gene_expr"][gid].astype(np.float32), clip_pct=99.0)
            gene_rows.append(pd.DataFrame({"pair": pair_label, "gene": gid, "gene_name": row["display"],
                                           "r": row["r"], "rank": rank, "l2_score": pd_["scores"],
                                           "expr_scaled": expr_sc}))

        de = pd.read_csv(DE_DIR / f"{pair_label}.tsv", sep="\t")
        cav = pd.read_csv(GENE_CORR_DIR / f"{pair_label}.tsv", sep="\t")
        merged = (de[["gene", "log2fc", "gene_name"]].rename(columns={"log2fc": "lfc"})
                  .merge(cav[["gene", "r"]].rename(columns={"r": "cav_r"}), on="gene", how="inner")
                  .dropna(subset=["lfc", "cav_r"]))
        merged["pair"] = pair_label
        de_rows.append(merged[["pair", "gene", "gene_name", "lfc", "cav_r"]].rename(columns={"lfc": "log2fc"}))

    out_dir = FIGURE_DATA_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(meta_rows).to_csv(out_dir / "sc_pair_meta.csv", index=False)
    pd.concat(cell_rows, ignore_index=True).to_csv(out_dir / "sc_continuum_cells.csv", index=False)
    pd.concat(gene_rows, ignore_index=True).to_csv(out_dir / "sc_continuum_genes.csv", index=False)
    pd.concat(de_rows, ignore_index=True).to_csv(out_dir / "sc_de_vs_cav.csv", index=False)
    print(f"Wrote sc_pair_meta.csv, sc_continuum_cells.csv, sc_continuum_genes.csv, "
          f"sc_de_vs_cav.csv to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
