#!/usr/bin/env python3
"""
export_projection_comparison_supp_data.py — Candidate supplemental figure
data: does orthogonalizing the condition CAV against tissue/cell-type
structure actually change the DE-vs-CAV relationship, compared to using the
raw (unprocessed) condition CAV directly?

For each of Figure 5's three case-study pairs (neutrophil/breast,
epithelial_cell/lung10x, fibroblast/colorectum), the same paired-donor cell
population used for the sc_de_vs_cav.csv scatter (panels E-G) is scored at
three stages of the SAME orthogonalization chain that ultimately produces L2:

  raw         — condition CAV direction, no orthogonalization at all
  L0_removed  — condition CAV with only the cell-type baseline (L0) removed
  L2          — condition CAV with cell-type (L0) AND tissue-within-cell-type
                (L1) removed — identical to the L2 axis already used
                everywhere else in Figure 5

Note there are only three mathematically distinct stages here (not four):
L0 and L1 are the only two things this hierarchy has to remove from a
condition CAV, so "L0-removed-only" and "L0+L1-removed" (=L2) are the full
set of intermediate points between raw and fully orthogonalized.

Uses the DEFAULT group-context-condition hierarchy (group=cell_type,
context=tissue) — same one export_fig5_case_studies.py and
paired_donor_gene_corr.py already use for L2 — so the "L2" row here should
reproduce results/gene_correlation_paired_donors/*.tsv and the existing
sc_de_vs_cav.csv rows exactly, which doubles as a sanity check on the
raw/L0_removed computation.

Output (into --figure-data-dir):
  proj_de_vs_cav.csv   pair, level, gene, gene_name, log2fc, cav_r
"""

import argparse
import os
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

_biocav_repo = os.environ.get("BIOCAV_REPO")
if not _biocav_repo:
    raise RuntimeError("BIOCAV_REPO is not set — source config/paths.sh before running this script.")
sys.path.insert(0, os.path.join(_biocav_repo, "core"))
sys.path.insert(0, os.path.join(_biocav_repo, "specific_scripts"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402
from src.utils.preprocessing import preprocess_embeddings   # noqa: E402
from cav_gene_correlation import correlate_scores_with_genes  # noqa: E402

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
CAVS = LIB / "cavs"
DE_DIR = LIB / "results" / "de_mixedlm"
COORDS_PATH = LIB / "results" / "hierarchy" / "cell_coordinates.tsv"

# pair_label, cell_type (obs value), tissue_real (obs value), condition (obs
# value), assay_filter (None if not needed), cav_tissue_token (second part
# of the CAV directory name — differs from tissue_real only for lung10x)
CASES = [
    ("neutrophil__breast__normal_vs_breast_cancer",
     "neutrophil", "breast", "breast cancer", None, "breast"),
    ("epithelial_cell__lung10x__normal_vs_lung_cancer",
     "epithelial cell", "lung", "lung cancer", "10x 3' v2", "lung10x"),
    ("fibroblast__colorectum__normal_vs_colorectal_cancer",
     "fibroblast", "colorectum", "colorectal cancer", None, "colorectum"),
]

# L0[group] baseline sets, matching exactly what the original hierarchy
# build used (see extend_hierarchy_lung10x.py for the epithelial_cell case —
# lung10x is deliberately excluded so L0 is reproduced unchanged).
L0_BASELINE_SETS = {
    "neutrophil": ["neutrophil__breast__normal", "neutrophil__lung__normal"],
    "epithelial_cell": ["epithelial_cell__colorectum__normal", "epithelial_cell__lung__normal",
                        "epithelial_cell__ovary__normal", "epithelial_cell__skin_epidermis__normal"],
    "fibroblast": ["fibroblast__colorectum__normal", "fibroblast__lung__normal",
                   "fibroblast__ovary__normal", "fibroblast__skin_epidermis__normal"],
}


def norm(v):
    return str(v).lower().replace("_", " ").strip()


def load_unit_cav(name):
    v = np.load(CAVS / name / f"concept_v1.npy").astype(np.float64)
    return v / np.linalg.norm(v)


def orthogonalize(v, basis):
    v = v.copy()
    for u in basis:
        v = v - np.dot(v, u) * u
    n = np.linalg.norm(v)
    return v / n if n > 1e-10 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure-data-dir", default="../figures/figure_data")
    args = ap.parse_args()
    out_dir = Path(args.figure_data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading shared scaler/PCA + all-cell embeddings...", file=sys.stderr)
    scaler = joblib.load(CAVS / "fibroblast__colorectum__normal" / "scaler_v1.pkl")
    pca = joblib.load(CAVS / "fibroblast__colorectum__normal" / "pca_v1.pkl")
    embs_raw, cell_ids = load_sequence_embeddings(str(LIB / "embeddings" / "cells.pkl"))
    cell_ids = np.array([str(c) for c in cell_ids])
    X_all = pca.transform(preprocess_embeddings(embs_raw, scaler))
    id_to_row = {c: i for i, c in enumerate(cell_ids)}

    print("Loading obs + coords (for sanity check)...", file=sys.stderr)
    adata = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r")
    obs = adata.obs[["cell_type", "tissue", "disease", "assay", "donor_id"]].copy()
    obs.index = obs.index.astype(str)
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["disease"] = obs["disease"].astype(str)
    gene_name_map = dict(zip(adata.var_names, adata.var["feature_name"].astype(str)))

    coords = pd.read_csv(COORDS_PATH, sep="\t", index_col=0)
    coords.index = coords.index.astype(str)

    rows = []

    for pair_label, ct, tissue_real, cond, assay_filter, cav_tissue in CASES:
        print(f"\n=== {pair_label} ===", file=sys.stderr)
        mask = ((obs["cell_type"].map(norm) == norm(ct)) &
                (obs["tissue"].map(norm) == norm(tissue_real)) &
                (obs["disease"].map(norm).isin(["normal", norm(cond)])))
        if assay_filter:
            mask = mask & (obs["assay"] == assay_filter)

        sub = obs[mask].copy()
        dc = sub.groupby("donor_id")["disease"].nunique()
        paired_donors = set(dc[dc == 2].index)
        sub = sub[sub["donor_id"].isin(paired_donors)]
        print(f"  {len(paired_donors)} paired donors, {len(sub)} cells", file=sys.stderr)

        cell_ids_sub = np.array([c for c in sub.index if c in id_to_row])
        rows_idx = [id_to_row[c] for c in cell_ids_sub]
        X_sub = X_all[rows_idx]

        group_us = ct.replace(" ", "_")
        cond_us = cond.replace(" ", "_")
        baseline_name = f"{group_us}__{cav_tissue}__normal"
        condition_name = f"{group_us}__{cav_tissue}__{cond_us}"

        baseline_dir = load_unit_cav(baseline_name)
        condition_dir = load_unit_cav(condition_name)

        L0_vecs = [load_unit_cav(n) for n in L0_BASELINE_SETS[group_us]]
        L0 = np.mean(L0_vecs, axis=0)
        L0 = L0 / np.linalg.norm(L0)

        L1 = orthogonalize(baseline_dir, [L0])

        raw_dir = condition_dir
        l0_removed_dir = orthogonalize(condition_dir, [L0])
        l2_dir = orthogonalize(condition_dir, [L0, L1])

        scores = {
            "raw": X_sub @ raw_dir,
            "L0_removed": X_sub @ l0_removed_dir,
            "L2": X_sub @ l2_dir,
        }

        # Sanity check: recomputed L2 should match the existing coords column.
        l2_col = f"L2__{group_us}__{tissue_real.replace(' ', '_')}__{cond_us}"
        if assay_filter:
            l2_col = f"L2__{group_us}__lung10x__{cond_us}"
        if l2_col in coords.columns:
            existing = coords.reindex(cell_ids_sub)[l2_col].values
            diff = np.nanmax(np.abs(scores["L2"] - existing))
            print(f"  Sanity check ({l2_col}): max |recomputed L2 - existing L2| = {diff:.6g}",
                  file=sys.stderr)
            if diff > 1e-3:
                print("  WARNING: mismatch — check L0 baseline set / CAV naming.", file=sys.stderr)
        else:
            print(f"  Sanity check skipped: {l2_col} not in coords columns", file=sys.stderr)

        # Gene expression for this cell subset (same lib-size + log1p convention
        # used throughout Figure 5's gene correlation).
        adata_sub = adata[cell_ids_sub, :]
        X_expr = adata_sub.X
        if sp.issparse(X_expr):
            X_expr = X_expr.toarray()
        X_expr = np.asarray(X_expr, dtype=np.float64)
        lib_size = X_expr.sum(axis=1, keepdims=True)
        lib_size[lib_size == 0] = 1.0
        X_expr = np.log1p(X_expr / lib_size * 1e4)
        var_names = list(adata_sub.var_names)

        de = pd.read_csv(DE_DIR / f"{pair_label}.tsv", sep="\t")

        for level, s in scores.items():
            corr_df = correlate_scores_with_genes(s, X_expr, var_names, min_cells=10)
            merged = (de[["gene", "log2fc", "gene_name"]]
                      .merge(corr_df[["gene", "r"]].rename(columns={"r": "cav_r"}),
                             on="gene", how="inner")
                      .dropna(subset=["log2fc", "cav_r"]))
            merged["pair"] = pair_label
            merged["level"] = level
            rows.append(merged[["pair", "level", "gene", "gene_name", "log2fc", "cav_r"]])
            print(f"  {level}: {len(merged)} genes", file=sys.stderr)

    out = pd.concat(rows, ignore_index=True)
    out_path = out_dir / "proj_de_vs_cav.csv"
    out.to_csv(out_path, index=False)
    print(f"\nWrote {out_path} ({len(out)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
