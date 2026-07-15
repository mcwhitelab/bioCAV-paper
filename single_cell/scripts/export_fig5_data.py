#!/usr/bin/env python3
"""
export_fig5_data.py — Export flat CSVs for Figure 5 (single-cell case studies)
so figures.R can render native ggplot panels instead of embedding the
matplotlib PNGs from cav_continuum_viz.py / de_pseudobulk.py.

Reuses cav_continuum_viz.py's data-assembly (get_pair_data) so the cell/gene
selection and expression scaling exactly match the existing exploratory
figures — this script only reshapes that same data into long CSVs plus the
DE-vs-CAV scatter table.

Outputs (into --figure-data-dir):
  sc_pair_meta.csv       pair, group, context, baseline, condition, n_normal, n_cancer
  sc_continuum_cells.csv pair, l2_score, is_baseline   (top strip)
  sc_continuum_genes.csv pair, gene, gene_name, r, rank, l2_score, expr_scaled  (gene strips, long: one row per cell x gene)
  sc_de_vs_cav.csv       pair, gene, gene_name, log2fc, cav_r   (unfiltered — R applies the lfc threshold)
"""

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from cav_continuum_viz import load_coords, load_obs, load_expression, get_pair_data, minmax  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coords", required=True)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--gene-corr-dir", required=True)
    ap.add_argument("--de-dir", required=True,
                     help="Directory of per-pair DE TSVs (gene, gene_name, log2fc, ...)")
    ap.add_argument("--pairs", required=True, nargs="+")
    ap.add_argument("--group-col", default="cell_type")
    ap.add_argument("--context-col", default="tissue")
    ap.add_argument("--condition-col", default="disease")
    ap.add_argument("--level", default="L2")
    ap.add_argument("--top-n-each", type=int, default=5)
    ap.add_argument("--max-cells", type=int, default=1000)
    ap.add_argument("--clip-pct", type=float, default=99.0)
    ap.add_argument("--figure-data-dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.figure_data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    coords = load_coords(args.coords, args.level)
    obs_cols = list({args.group_col, args.context_col, args.condition_col})
    obs = load_obs(args.h5ad, obs_cols)
    X_expr, gene_ids, expr_cell_ids, gene_name_map = load_expression(args.h5ad)

    # cav_continuum_viz.load_expression() returns raw counts. Cells vary
    # hugely in total library size, so raw counts confound a gene's true
    # relative expression with how deeply that particular cell was
    # sequenced -- e.g. a handful of very-high-depth cancer cells can show
    # a large raw count for a gene that is actually a *smaller* share of
    # their transcriptome than in normal cells. cav_gene_correlation.py
    # (which computes the "r" values these continuum strips are labeled
    # with) normalizes for exactly this reason -- match its convention
    # here so the strip coloring is consistent with the displayed r.
    lib_size = X_expr.sum(axis=1, keepdims=True)
    lib_size[lib_size == 0] = 1.0
    X_expr = np.log1p(X_expr / lib_size * 1e4)

    gene_corr_dir = Path(args.gene_corr_dir)
    de_dir = Path(args.de_dir)

    meta_rows = []
    cell_rows = []
    gene_rows = []
    de_rows = []

    for pair_label in args.pairs:
        pd_ = get_pair_data(
            pair_label=pair_label, coords=coords, X_expr=X_expr,
            gene_ids=gene_ids, gene_name_map=gene_name_map, obs=obs,
            group_col=args.group_col, context_col=args.context_col,
            condition_col=args.condition_col, level=args.level,
            gene_corr_dir=gene_corr_dir, top_n=8, max_cells=args.max_cells,
            baseline_value="normal", top_n_each=args.top_n_each,
        )
        if pd_ is None:
            print(f"WARNING: skipping {pair_label} (no data)", file=sys.stderr)
            continue

        meta_rows.append({
            "pair": pair_label, "group": pd_["group"], "context": pd_["context"],
            "baseline": pd_["baseline"], "condition": pd_["condition"],
            "n_normal": int(pd_["n_normal"]), "n_cancer": int(pd_["n_cancer"]),
        })

        is_base = pd_["labels"] == pd_["baseline"]
        cell_rows.append(pd.DataFrame({
            "pair": pair_label,
            "l2_score": pd_["scores"],
            "is_baseline": is_base,
        }))

        # gene_df comes back as [best-positive..weakest-positive, most-negative..
        # least-negative] (matching draw_gene_strips' row order in the PNG:
        # positives ranked best-first, then negatives ranked strongest-first).
        # A flat sort by r would interleave the two blocks, so re-derive it
        # explicitly rather than trusting a single sort key.
        pos = pd_["gene_df"][pd_["gene_df"]["r"] >= 0].sort_values("r", ascending=False)
        neg = pd_["gene_df"][pd_["gene_df"]["r"] < 0].sort_values("r", ascending=True)
        gene_df = pd.concat([pos, neg]).reset_index(drop=True)
        for rank, row in gene_df.iterrows():
            gid = row["gene"]
            if gid not in pd_["gene_expr"]:
                continue
            expr_sc = minmax(pd_["gene_expr"][gid].astype(np.float32), clip_pct=args.clip_pct)
            gene_rows.append(pd.DataFrame({
                "pair": pair_label,
                "gene": gid,
                "gene_name": row["display"],
                "r": row["r"],
                "rank": rank,
                "l2_score": pd_["scores"],
                "expr_scaled": expr_sc,
            }))

        # DE vs CAV scatter table: merge DE log2fc with CAV r for every shared gene
        de = pd.read_csv(de_dir / f"{pair_label}.tsv", sep="\t")
        cav = pd.read_csv(gene_corr_dir / f"{pair_label}.tsv", sep="\t")
        de_cols = ["gene", "log2fc"] + (["gene_name"] if "gene_name" in de.columns else [])
        cav_cols = ["gene", "r"] + (["gene_name"] if "gene_name" in cav.columns else [])
        merged = (de[de_cols].rename(columns={"log2fc": "lfc"})
                  .merge(cav[cav_cols].rename(columns={"r": "cav_r", "gene_name": "gene_name_cav"}),
                         on="gene", how="inner")
                  .dropna(subset=["lfc", "cav_r"]))
        if "gene_name" in merged.columns:
            merged["gene_name"] = merged["gene_name"].fillna(merged.get("gene_name_cav", merged["gene"]))
        elif "gene_name_cav" in merged.columns:
            merged["gene_name"] = merged["gene_name_cav"].fillna(merged["gene"])
        else:
            merged["gene_name"] = merged["gene"]
        merged["pair"] = pair_label
        de_rows.append(merged[["pair", "gene", "gene_name", "lfc", "cav_r"]]
                        .rename(columns={"lfc": "log2fc"}))

    pd.DataFrame(meta_rows).to_csv(out_dir / "sc_pair_meta.csv", index=False)
    pd.concat(cell_rows, ignore_index=True).to_csv(out_dir / "sc_continuum_cells.csv", index=False)
    pd.concat(gene_rows, ignore_index=True).to_csv(out_dir / "sc_continuum_genes.csv", index=False)
    pd.concat(de_rows, ignore_index=True).to_csv(out_dir / "sc_de_vs_cav.csv", index=False)
    print(f"Wrote sc_pair_meta.csv, sc_continuum_cells.csv, sc_continuum_genes.csv, "
          f"sc_de_vs_cav.csv to {out_dir}")


if __name__ == "__main__":
    main()
