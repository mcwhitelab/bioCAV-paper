#!/usr/bin/env python3
"""
export_projection_comparison_supp_data.py — Candidate supplemental figure
data: does orthogonalizing the condition CAV against tissue/cell-type
structure actually change the DE-vs-CAV relationship, compared to using the
raw (unprocessed) condition CAV directly?

For each of Figure 6's three case-study pairs (neutrophil/breast,
epithelial_cell/lung10x, neutrophil/lung since 2026-10-08), the same paired-donor cell
population used for the sc_de_vs_cav.csv scatter is scored at the stages of the
orthogonalization chain that produces L2, plus the axis Fig 6 now uses:

  raw         — condition CAV direction, no orthogonalization at all
  L0_removed  — condition CAV with only the cell-type baseline (L0) removed
  L2          — condition CAV with cell-type (L0) AND tissue-within-cell-type
                (L1) removed — identical to the L2 axis already used
                in the library hierarchy
  diff        — unit(condition CAV) - unit(matched normal CAV), no projection: the Fig 6
                axis since 2026-10-08 (single_cell/official_tokenizer/16a_build_diff_coords.py)

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

--run-dir <official_tokenizer run>  takes scaler.pkl, pca.pkl, cavs/ and the library
embeddings (--readout, default cls_penult) from an official-tokenizer run instead of the
library; the L2 / diff rows are checked against <run-dir>/coords.tsv / coords_diff.tsv.

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
# DE axis = the same pseudobulk DESeq2 table as the Figure 6 scatter
# (export_fig5_case_studies.py): paired donors, ~ donor + condition, samples
# < 10 cells dropped, unshrunk log2FoldChange (column log2fc_mle). With the
# conv1 CAVs the L2 row therefore reproduces Fig 6 exactly. (Before
# 2026-10-07: per-cell mixed model, results/de_mixedlm_paired_donors.)
DE_DIR = LIB / "results" / "de_pseudobulk_paired_donors" / "min10"
# all ten paired-donor pairs, same recipe (identical for the shared pairs); has neutrophil/lung
DE_DIR_ALL = Path("convergence_check/de_pseudobulk_allpairs/min10")
XD_EMB = Path("/xdisk/clairemcwhite/clairemcwhite/sc_geneformer_official/emb/library")
COORDS_PATH = LIB / "results" / "hierarchy" / "cell_coordinates.tsv"

# pair_label, cell_type (obs value), tissue_real (obs value), condition (obs
# value), assay_filter (None if not needed), cav_tissue_token (second part
# of the CAV directory name — differs from tissue_real only for lung10x)
CASES = [
    ("neutrophil__breast__normal_vs_breast_cancer",
     "neutrophil", "breast", "breast cancer", None, "breast"),
    ("epithelial_cell__lung10x__normal_vs_lung_cancer",
     "epithelial cell", "lung", "lung cancer", "10x 3' v2", "lung10x"),
    ("neutrophil__lung__normal_vs_lung_cancer",
     "neutrophil", "lung", "lung cancer", None, "lung"),
]
# previous third case (until 2026-10-08):
# ("fibroblast__colorectum__normal_vs_colorectal_cancer", "fibroblast", "colorectum", "colorectal cancer", None, "colorectum")

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


# Concept vectors are read from CAV_VEC_DIR (default: the library's stored CAVs). Since
# 2026-10-07 the supplement uses convergence_check/cavs_conv1 (C=1, PCA-128, refit to lbfgs
# convergence), matching Fig 6. Scaler/PCA always come from the library (shared, unchanged).
CAV_VEC_DIR = CAVS


def load_unit_cav(name):
    v = np.load(CAV_VEC_DIR / name / f"concept_v1.npy").astype(np.float64)
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
    ap.add_argument("--cav-dir", default=None,
                    help="directory of <cav>/concept_v1.npy to use instead of the library's stored CAVs")
    ap.add_argument("--coords", default=None,
                    help="hierarchy coordinates for the L2 sanity check (default: library cell_coordinates.tsv)")
    ap.add_argument("--run-dir", default=None, help="official-tokenizer run (scaler, pca, cavs, coords)")
    ap.add_argument("--readout", default="cls_penult", help="embedding readout for --run-dir")
    args = ap.parse_args()
    global CAV_VEC_DIR
    run = Path(args.run_dir) if args.run_dir else None
    if args.cav_dir:
        CAV_VEC_DIR = Path(args.cav_dir)
    elif run:
        CAV_VEC_DIR = run / "cavs"
    out_dir = Path(args.figure_data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading shared scaler/PCA + all-cell embeddings...", file=sys.stderr)
    if run:
        scaler, pca = joblib.load(run / "scaler.pkl"), joblib.load(run / "pca.pkl")
        p = XD_EMB / f"{args.readout}.pkl"
        paths = [p] if p.exists() else sorted(XD_EMB.glob(f"shard*/{args.readout}.pkl"))
        parts = [load_sequence_embeddings(str(q)) for q in paths]
        embs_raw = np.vstack([e for e, _ in parts]); cell_ids = [c for _, i in parts for c in i]
    else:
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

    coords = pd.read_csv(args.coords or (run / "coords.tsv" if run else COORDS_PATH), sep="\t", index_col=0)
    coords.index = coords.index.astype(str)
    coords_diff = None
    if run and (run / "coords_diff.tsv").exists():
        coords_diff = pd.read_csv(run / "coords_diff.tsv", sep="\t", index_col=0)
        coords_diff.index = coords_diff.index.astype(str)

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

        diff_dir = condition_dir - baseline_dir
        diff_dir = diff_dir / np.linalg.norm(diff_dir)
        scores = {
            "raw": X_sub @ raw_dir,
            "L0_removed": X_sub @ l0_removed_dir,
            "L2": X_sub @ l2_dir,
            "diff": X_sub @ diff_dir,
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
        if coords_diff is not None and l2_col in coords_diff.columns:
            dd = np.nanmax(np.abs(scores["diff"] - coords_diff.reindex(cell_ids_sub)[l2_col].values))
            print(f"  Sanity check (diff vs coords_diff.tsv): max |diff| = {dd:.6g}", file=sys.stderr)
            if dd > 1e-3:
                print("  WARNING: diff axis mismatch.", file=sys.stderr)

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

        de_dir = DE_DIR if (DE_DIR / f"{pair_label}.tsv").exists() else DE_DIR_ALL
        de = (pd.read_csv(de_dir / f"{pair_label}.tsv", sep="\t")
              .drop(columns="log2fc").rename(columns={"log2fc_mle": "log2fc"}))

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
