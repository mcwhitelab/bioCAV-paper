#!/usr/bin/env python3
"""
train_10x_lung_cavs.py — Train epithelial_cell__lung10x__{normal,lung_cancer}
CAVs restricted to assay == "10x 3' v2" cells, to avoid the Smart-seq2/10x
sequencing-technology confound in the original epithelial_cell__lung__*
CAVs (which mixed multiple assays, ~91% 10x on normal vs ~50% on cancer).

Replicates the pre-1489e3a training convention (scale -> shared 128-dim
global PCA -> logistic regression, C=1.0, class_weight=balanced, 5-fold CV)
so the new CAVs are dimensionally and spatially compatible with the rest of
this library's baseline CAVs -- reuses the *existing* shared scaler_v1.pkl,
pca_v1.pkl, and reference-negative pool (neg.npy/ref_neg_ids.txt) rather
than refitting anything, for exact consistency.

New tissue label "lung10x" (distinct from "lung") keeps this from
clobbering or silently altering the original mixed-assay lung CAVs.
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, average_precision_score

_biocav_repo = os.environ.get("BIOCAV_REPO")
if not _biocav_repo:
    raise RuntimeError("BIOCAV_REPO is not set — source config/paths.sh before running this script.")
sys.path.insert(0, os.path.join(_biocav_repo, "core"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
REFERENCE_CAV_DIR = LIB / "cavs" / "epithelial_cell__lung__normal"  # source of shared scaler/pca/neg pool
N_POS_CAP = 500
SEED = 42


def norm(v):
    return str(v).lower().replace("_", " ").strip()


def load_pos_embeddings(cell_ids_wanted, embs_raw, all_cell_ids):
    id_to_row = {cid: i for i, cid in enumerate(all_cell_ids)}
    rows = [id_to_row[c] for c in cell_ids_wanted if c in id_to_row]
    missing = len(cell_ids_wanted) - len(rows)
    if missing:
        print(f"  WARNING: {missing} requested cell ids not found in embeddings", file=sys.stderr)
    return embs_raw[rows]


def train_and_save(pos_emb, neg_emb, scaler, pca, out_dir: Path, config: dict):
    out_dir.mkdir(parents=True, exist_ok=True)

    X_raw = np.vstack([pos_emb, neg_emb])
    y = np.hstack([np.ones(len(pos_emb)), np.zeros(len(neg_emb))])

    X_scaled = scaler.transform(X_raw)
    X_pca = pca.transform(X_scaled)

    clf = LogisticRegression(C=config["regularization_C"], max_iter=1000,
                             random_state=SEED, solver="lbfgs",
                             class_weight=config["class_weight"])
    cv = StratifiedKFold(n_splits=config["cv_folds"], shuffle=True, random_state=SEED)
    cv_scores = cross_val_score(clf, X_pca, y, cv=cv, scoring="roc_auc")
    clf.fit(X_pca, y)

    y_pred = clf.predict_proba(X_pca)[:, 1]
    train_auroc = roc_auc_score(y, y_pred)
    train_auprc = average_precision_score(y, y_pred)

    concept = clf.coef_[0]
    concept = concept / (np.linalg.norm(concept) + 1e-10)

    np.save(out_dir / "concept_v1.npy", concept)
    np.save(out_dir / "pos.npy", pos_emb)
    np.save(out_dir / "neg.npy", neg_emb)
    joblib.dump(scaler, out_dir / "scaler_v1.pkl")
    joblib.dump(pca, out_dir / "pca_v1.pkl")

    report = {
        "artifact_version": "v1",
        "timestamp": datetime.now().isoformat(),
        "concept_metrics": {
            "cv_auroc_mean": float(cv_scores.mean()),
            "cv_auroc_std": float(cv_scores.std()),
            "cv_auroc_scores": cv_scores.tolist(),
            "train_auroc": float(train_auroc),
            "train_auprc": float(train_auprc),
            "n_train": len(X_pca),
            "n_positive": int(y.sum()),
            "n_negative": int((1 - y).sum()),
            "global_pca": True,
            "pca_dim": pca.n_components,
        },
        "config": config,
        "note": "assay-restricted retrain (10x 3' v2 only), reusing shared scaler/pca/ref-neg pool from "
                + str(REFERENCE_CAV_DIR),
    }
    (out_dir / "report_v1.json").write_text(json.dumps(report, indent=2))

    manifest = {
        "version": "v1", "created": datetime.now().isoformat(),
        "files": {k: {"filename": f, "size_bytes": (out_dir / f).stat().st_size, "exists": True}
                 for k, f in [("concept_cav", "concept_v1.npy"), ("scaler", "scaler_v1.pkl"),
                              ("pca", "pca_v1.pkl"), ("report", "report_v1.json")]},
    }
    (out_dir / "manifest_v1.json").write_text(json.dumps(manifest, indent=2))

    print(f"  Saved {out_dir}  CV AUROC={cv_scores.mean():.3f}±{cv_scores.std():.3f}  "
          f"train AUROC={train_auroc:.3f}  n_pos={int(y.sum())} n_neg={int((1-y).sum())}")


def main():
    print("Loading shared scaler/PCA/reference-negative pool...", file=sys.stderr)
    scaler = joblib.load(REFERENCE_CAV_DIR / "scaler_v1.pkl")
    pca = joblib.load(REFERENCE_CAV_DIR / "pca_v1.pkl")
    neg_emb = np.load(REFERENCE_CAV_DIR / "neg.npy")
    ref_neg_ids_text = (REFERENCE_CAV_DIR / "ref_neg_ids.txt").read_text()

    print("Loading h5ad obs + embeddings...", file=sys.stderr)
    adata = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r")
    obs = adata.obs[["cell_type", "tissue", "disease", "assay"]].copy()
    obs.index = obs.index.astype(str)
    adata.file.close()

    embs_raw, cell_ids = load_sequence_embeddings(str(LIB / "embeddings" / "cells.pkl"))
    cell_ids = [str(c) for c in cell_ids]

    for disease_val, out_name in [("normal", "epithelial_cell__lung10x__normal"),
                                   ("lung cancer", "epithelial_cell__lung10x__lung_cancer")]:
        mask = ((obs["cell_type"].map(norm) == "epithelial cell") &
                (obs["tissue"].map(norm) == "lung") &
                (obs["disease"].map(norm) == norm(disease_val)) &
                (obs["assay"] == "10x 3' v2"))
        pos_ids = obs.index[mask].tolist()
        rng = np.random.default_rng(SEED)
        if len(pos_ids) > N_POS_CAP:
            pos_ids = list(rng.choice(pos_ids, N_POS_CAP, replace=False))
        print(f"\n{out_name}: {len(pos_ids)} positive (10x 3' v2) cells", file=sys.stderr)

        pos_emb = load_pos_embeddings(pos_ids, embs_raw, cell_ids)

        out_dir = LIB / "cavs" / out_name
        train_and_save(pos_emb, neg_emb, scaler, pca, out_dir,
                       config={"regularization_C": 1.0, "cv_folds": 5,
                               "random_seed": SEED, "class_weight": "balanced"})
        (out_dir / "ref_neg_ids.txt").write_text(ref_neg_ids_text)


if __name__ == "__main__":
    main()
