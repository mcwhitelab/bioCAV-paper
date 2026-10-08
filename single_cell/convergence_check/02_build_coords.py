#!/usr/bin/env python3
"""Rebuild L0/L1/L2 hierarchy coordinates for each refit arm, in the PCA-128
frame, using cav_hierarchy.build_hierarchy for the 91 original CAVs and the
extend_hierarchy_lung10x.py rule for the lung10x pair (L0[epithelial] from
the 4 original baselines only). Gate: the 'orig' arm must reproduce
results/hierarchy/cell_coordinates.tsv.
Output: coords_<arm>.tsv (same column names as the library file).
"""
import os, sys, logging
from pathlib import Path
import numpy as np, pandas as pd, joblib
sys.path.insert(0, os.path.join(os.environ["BIOCAV_REPO"], "core"))
sys.path.insert(0, os.path.join(os.environ["BIOCAV_REPO"], "specific_scripts"))
from src.utils.data_loader import load_sequence_embeddings
from cav_hierarchy import build_hierarchy, orthogonalize
logging.basicConfig(level=logging.WARNING)

LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
ref = LIB / "cavs" / "epithelial_cell__lung__normal"
scaler, pca = joblib.load(ref / "scaler_v1.pkl"), joblib.load(ref / "pca_v1.pkl")
embs, cell_ids = load_sequence_embeddings(str(LIB / "embeddings" / "cells.pkl"))
cell_ids = [str(c) for c in cell_ids]
Xs = scaler.transform(embs)
X_pca = pca.transform(Xs)
# arms named np_* live in the scaler frame (no PCA); others in PCA-128
frame = lambda arm: Xs if arm.startswith("np_") else X_pca

existing = pd.read_csv(LIB / "results" / "hierarchy" / "cell_coordinates.tsv", sep="\t", index_col=0)
existing.index = existing.index.astype(str)
existing = existing.reindex(cell_ids)

def unit(p):
    v = np.load(p).astype(np.float64); return v / np.linalg.norm(v)

for arm in sys.argv[1:] or ["orig", "conv1", "conv4"]:
    cdir = Path(f"cavs_{arm}")
    names = {d.name: tuple(d.name.split("__")) for d in sorted(cdir.iterdir()) if "lung10x" not in d.name}
    dirs = {n: unit(cdir / n / "concept_v1.npy") for n in names}
    H = build_hierarchy(names, dirs, lib_dir=LIB)
    cols = {}
    for ct, v in sorted(H["L0"].items()): cols[f"L0__{ct}"] = v
    for (ct, t), v in sorted(H["L1"].items()): cols[f"L1__{ct}__{t}"] = v
    for (ct, t, d), v in sorted(H["L2"].items()): cols[f"L2__{ct}__{t}__{d}"] = v
    # lung10x, as extend_hierarchy_lung10x.py
    L0e = H["L0"]["epithelial_cell"]
    L1 = orthogonalize(unit(cdir / "epithelial_cell__lung10x__normal" / "concept_v1.npy"), [L0e])
    cols["L1__epithelial_cell__lung10x"] = L1
    cols["L2__epithelial_cell__lung10x__lung_cancer"] = orthogonalize(
        unit(cdir / "epithelial_cell__lung10x__lung_cancer" / "concept_v1.npy"), [L0e, L1])
    np.savez(f"axes_{arm}.npz", **cols)
    if not arm.startswith("np_"):  # same axes mapped into the scaler frame, for cross-frame cosines
        np.savez(f"axes_{arm}_sc.npz", **{k: pca.components_.T @ v for k, v in cols.items()})
    X = frame(arm)
    out = pd.DataFrame({k: X @ v for k, v in cols.items()}, index=pd.Index(cell_ids, name="cell_id"))
    if arm in ("orig", "stored"):
        common = [c for c in out.columns if c in existing.columns]
        missing = sorted(set(existing.columns) - set(out.columns)); extra = sorted(set(out.columns) - set(existing.columns))
        diffs = {c: np.nanmax(np.abs(out[c].values - existing[c].values)) for c in common}
        worst = max(diffs, key=diffs.get)
        print(f"GATE orig: {len(common)} common cols, missing={missing}, extra={extra}, "
              f"max|diff|={diffs[worst]:.3g} ({worst})", flush=True)
    out.to_csv(f"coords_{arm}.tsv", sep="\t")
    print(f"wrote coords_{arm}.tsv {out.shape}", flush=True)
