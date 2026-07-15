#!/usr/bin/env python3
"""
extend_hierarchy_lung10x.py — Surgically add the new lung10x pair to
results/hierarchy/cell_coordinates.tsv without touching any other axis.

L0[epithelial_cell] is recomputed exactly as the original build did (mean
of the 4 original epithelial baseline CAVs: colorectum, lung, ovary,
skin_epidermis -- NOT including lung10x), so every existing column in
cell_coordinates.tsv is reproduced unchanged. Only three new columns are
added:
  L1__epithelial_cell__lung10x
  L2__epithelial_cell__lung10x__normal
  L2__epithelial_cell__lung10x__lung_cancer
"""

import argparse
import os
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

_biocav_repo = os.environ.get("BIOCAV_REPO")
sys.path.insert(0, os.path.join(_biocav_repo, "core"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402
from src.utils.preprocessing import preprocess_embeddings   # noqa: E402

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
CAVS = LIB / "cavs"
ORIGINAL_EPITHELIAL_BASELINES = [
    "epithelial_cell__colorectum__normal",
    "epithelial_cell__lung__normal",
    "epithelial_cell__ovary__normal",
    "epithelial_cell__skin_epidermis__normal",
]
TARGET_COORDS = [
    LIB / "results" / "hierarchy" / "cell_coordinates.tsv",
    LIB / "results" / "hierarchy_celltype_tissue" / "cell_coordinates.tsv",
]


def load_unit_cav(name):
    v = np.load(CAVS / name / "concept_v1.npy").astype(np.float64)
    return v / np.linalg.norm(v)


def orthogonalize(v, basis):
    v = v.copy()
    for u in basis:
        v -= np.dot(v, u) * u
    norm = np.linalg.norm(v)
    return v / norm if norm > 1e-10 else None


def main():
    print("Recomputing L0[epithelial_cell] (unchanged from original build)...", file=sys.stderr)
    baseline_vecs = [load_unit_cav(n) for n in ORIGINAL_EPITHELIAL_BASELINES]
    L0_epithelial = np.mean(baseline_vecs, axis=0)
    L0_epithelial = L0_epithelial / np.linalg.norm(L0_epithelial)

    print("Computing new lung10x L1/L2 axes...", file=sys.stderr)
    lung10x_baseline = load_unit_cav("epithelial_cell__lung10x__normal")
    lung10x_condition = load_unit_cav("epithelial_cell__lung10x__lung_cancer")

    L1_lung10x = orthogonalize(lung10x_baseline, [L0_epithelial])
    # The baseline's own L2 residual is always exactly zero by construction
    # (L1 IS baseline's residual after removing L0, so removing [L0, L1]
    # from baseline leaves nothing) -- matches the original pipeline, which
    # skips this case rather than emitting a degenerate axis.
    L2_cancer = orthogonalize(lung10x_condition, [L0_epithelial, L1_lung10x])

    print("Loading scaler/PCA and all-cell embeddings...", file=sys.stderr)
    scaler = joblib.load(LIB / "cavs" / "epithelial_cell__lung__normal" / "scaler_v1.pkl")
    pca = joblib.load(LIB / "cavs" / "epithelial_cell__lung__normal" / "pca_v1.pkl")
    embs_raw, cell_ids = load_sequence_embeddings(str(LIB / "embeddings" / "cells.pkl"))
    cell_ids = [str(c) for c in cell_ids]

    X = pca.transform(preprocess_embeddings(embs_raw, scaler))

    new_cols = pd.DataFrame({
        "cell_id": cell_ids,
        "L1__epithelial_cell__lung10x": X @ L1_lung10x,
        "L2__epithelial_cell__lung10x__lung_cancer": X @ L2_cancer,
    }).set_index("cell_id")

    for coords_path in TARGET_COORDS:
        backup_path = coords_path.with_suffix(".tsv.bak")
        source_path = backup_path if backup_path.exists() else coords_path
        print(f"\nLoading {source_path}...", file=sys.stderr)
        coords = pd.read_csv(source_path, sep="\t", index_col=0)
        coords.index = coords.index.astype(str)

        # Sanity check: recomputed L0[epithelial_cell] should reproduce the
        # existing column exactly (up to float precision), confirming we used
        # the right baseline set / scaler / pca.
        if "L0__epithelial_cell" in coords.columns:
            check = X @ L0_epithelial
            existing = coords.reindex(cell_ids)["L0__epithelial_cell"].values
            max_diff = np.nanmax(np.abs(check - existing))
            print(f"  Sanity check: max |recomputed L0 - existing L0| = {max_diff:.6g}", file=sys.stderr)
            if max_diff > 1e-3:
                print("  WARNING: recomputed L0 does not match existing column -- "
                      "investigate before trusting the new axes.", file=sys.stderr)

        merged = coords.join(new_cols, how="left")

        if not backup_path.exists():
            coords_path.rename(backup_path)
            print(f"  Backed up original to {backup_path}", file=sys.stderr)
        merged.to_csv(coords_path, sep="\t")
        print(f"  Wrote extended {coords_path} ({merged.shape[0]} cells x {merged.shape[1]} axes)", file=sys.stderr)


if __name__ == "__main__":
    main()
