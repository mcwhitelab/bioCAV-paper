#!/usr/bin/env python3
"""
export_fig5_umap_data.py — Export a flat CSV for the overview row at the top
of Figure 5: t-SNE of cells at four stages of CAV processing --
  raw — every individual baseline CAV direction, unmodified, no averaging
        or orthogonalization at all (one feature per group__tissue__normal
        CAV, e.g. 31 of them). The rawest possible starting point.
  L0  — primary/group axes      (tissue, per cav_hierarchy.py's
                                  context-group-condition level order --
                                  averages the raw baseline CAVs within
                                  each tissue)
  L1  — secondary/context axes  (cell type, tissue projected away)
  L2  — condition-residual axes (normal vs. cancer, tissue + cell type
                                  projected away)

Each level gets its own t-SNE fit (same reducer/params cav_viz.py's
--reducer tsne uses) on the same subsampled cells, colored by whatever that
level's axes are actually built to separate.

Output (into --figure-data-dir):
  sc_umap_overview.csv   level, cell_id, x, y, cell_type, tissue, is_baseline
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE

_biocav_repo = os.environ.get("BIOCAV_REPO")
if not _biocav_repo:
    raise RuntimeError("BIOCAV_REPO is not set — source config/paths.sh before running this script.")
sys.path.insert(0, os.path.join(_biocav_repo, "core"))
sys.path.insert(0, os.path.join(_biocav_repo, "specific_scripts"))
from src.utils.data_loader import load_sequence_embeddings  # noqa: E402
from src.utils.preprocessing import preprocess_embeddings   # noqa: E402
from cav_hierarchy import load_cav_direction                # noqa: E402


def make_tsne(seed: int):
    return TSNE(n_components=2, random_state=seed, perplexity=30,
                max_iter=1000, init="pca", learning_rate="auto")


def project_onto_raw_cavs(lib_dir: Path, cell_ids: list, X_scaled_by_id: dict,
                          baseline_value: str) -> pd.DataFrame:
    """Every individual baseline CAV, dot-producted against scaled embeddings.
    No averaging, no orthogonalization -- one column per group__tissue CAV.

    These CAVs were trained in a PCA-reduced space (concept_v1.npy is
    128-dim, not the 768-dim scaler space cav_hierarchy.py's docstring
    assumes) -- confirmed by report_v1.json's "global_pca": true and by
    the pca_v1.pkl objects being byte-identical across different CAV
    dirs. So there's one more step than L0/L1/L2's own scoring: scale ->
    (shared) PCA -> dot with the 128-dim concept vector.
    """
    import joblib

    cavs_dir = lib_dir / "cavs"
    baseline_dirs = sorted(d.name for d in cavs_dir.iterdir()
                           if d.is_dir() and d.name.endswith(f"__{baseline_value}"))
    print(f"  {len(baseline_dirs)} individual baseline CAVs found", file=sys.stderr)

    shared_pca = joblib.load(cavs_dir / baseline_dirs[0] / "pca_v1.pkl")
    X_scaled = np.stack([X_scaled_by_id[cid] for cid in cell_ids])
    X = shared_pca.transform(X_scaled)

    feats = {}
    for name in baseline_dirs:
        try:
            unit_vec = load_cav_direction(cavs_dir / name)
        except Exception as e:
            print(f"  WARNING: could not load {name}: {e}", file=sys.stderr)
            continue
        feats[name] = X @ unit_vec
    return pd.DataFrame(feats, index=cell_ids)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coords", required=True, help="cell_coordinates.tsv from cav_hierarchy.py")
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--lib-dir", required=True,
                    help="CAV library dir (contains cavs/ and global_scaler_v1.pkl)")
    ap.add_argument("--raw-pkl", required=True, help="Raw pre-CAV embeddings .pkl")
    ap.add_argument("--scaler-pkl", default=None,
                    help="Defaults to <lib-dir>/global_scaler_v1.pkl")
    ap.add_argument("--group-col", default="cell_type")
    ap.add_argument("--context-col", default="tissue")
    ap.add_argument("--condition-col", default="disease")
    ap.add_argument("--baseline-value", default="normal")
    ap.add_argument("--max-cells", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--figure-data-dir", required=True)
    ap.add_argument("--clean-pairs-csv", default=None,
                    help="Optional CSV with cell_type,tissue columns -- if given, "
                         "restrict cells to only these (assay-matched) combos before "
                         "subsampling, so condition separation in the t-SNE can't be "
                         "a sequencing-technology batch effect.")
    ap.add_argument("--paired-donors-only", action="store_true",
                    help="Within each (cell_type, tissue) group, keep only donors that "
                         "have both the baseline and non-baseline condition present -- "
                         "same paired-donor restriction used for panels E-G's gene "
                         "correlations, so the overview t-SNE reflects the same cell "
                         "population rather than a broader any-donor pool.")
    ap.add_argument("--donor-col", default="donor_id")
    args = ap.parse_args()

    import anndata as ad
    import joblib

    lib_dir = Path(args.lib_dir)
    scaler_pkl = Path(args.scaler_pkl) if args.scaler_pkl else lib_dir / "global_scaler_v1.pkl"

    print("Loading CAV-projected coordinates...", file=sys.stderr)
    coords = pd.read_csv(args.coords, sep="\t", index_col=0)
    coords.index = coords.index.astype(str)

    print("Loading obs metadata...", file=sys.stderr)
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs_cols = [args.group_col, args.context_col, args.condition_col]
    if args.paired_donors_only:
        obs_cols.append(args.donor_col)
    obs = adata.obs[obs_cols].copy()
    adata.file.close()
    obs.index = obs.index.astype(str)

    print("Loading raw embeddings...", file=sys.stderr)
    embs_raw, cell_ids = load_sequence_embeddings(args.raw_pkl)
    cell_ids = [str(c) for c in cell_ids]

    shared = sorted(set(coords.index) & set(obs.index) & set(cell_ids))
    print(f"{len(shared)} cells shared across coords, obs, and raw embeddings", file=sys.stderr)

    if args.clean_pairs_csv:
        clean = pd.read_csv(args.clean_pairs_csv)
        allowed = set(zip(clean["cell_type"].str.lower().str.strip(),
                          clean["tissue"].str.lower().str.strip()))
        obs_key = list(zip(obs[args.group_col].astype(str).str.lower().str.strip(),
                           obs[args.context_col].astype(str).str.lower().str.strip()))
        obs_clean_mask = pd.Series([k in allowed for k in obs_key], index=obs.index)
        shared = sorted(set(shared) & set(obs.index[obs_clean_mask]))
        print(f"  Restricted to {len(allowed)} assay-matched (cell_type, tissue) combos: "
              f"{len(shared)} cells remain", file=sys.stderr)

    if args.paired_donors_only:
        is_baseline_col = obs[args.condition_col].astype(str).str.lower().str.contains(
            args.baseline_value.lower())
        donor = obs[args.donor_col].astype(str)
        group_key = list(zip(obs[args.group_col].astype(str), obs[args.context_col].astype(str), donor))
        cond_by_key = {}
        for k, b in zip(group_key, is_baseline_col):
            cond_by_key.setdefault(k, set()).add(b)
        paired_keys = {k for k, v in cond_by_key.items() if len(v) == 2}
        paired_mask = pd.Series([k in paired_keys for k in group_key], index=obs.index)
        shared = sorted(set(shared) & set(obs.index[paired_mask]))
        print(f"  Restricted to paired donors within each (cell_type, tissue) group: "
              f"{len(shared)} cells remain", file=sys.stderr)

    rng = np.random.default_rng(args.seed)
    n = min(args.max_cells, len(shared))
    chosen = list(rng.choice(shared, n, replace=False))

    meta = obs.loc[chosen]
    is_baseline = meta[args.condition_col].astype(str).str.lower().str.contains(
        args.baseline_value.lower())

    rows = []

    # --- raw level: every individual baseline CAV, unmodified ---
    print("Loading global scaler...", file=sys.stderr)
    scaler = joblib.load(scaler_pkl)["scaler"]
    id_to_row = {cid: i for i, cid in enumerate(cell_ids)}
    chosen_rows = np.array([id_to_row[c] for c in chosen])
    X_scaled = preprocess_embeddings(embs_raw[chosen_rows], scaler)
    X_scaled_by_id = {cid: X_scaled[i] for i, cid in enumerate(chosen)}

    raw_feats = project_onto_raw_cavs(lib_dir, chosen, X_scaled_by_id, args.baseline_value)
    print(f"Running t-SNE on raw ({raw_feats.shape})...", file=sys.stderr)
    emb = make_tsne(args.seed).fit_transform(raw_feats.values)
    rows.append(pd.DataFrame({
        "level": "raw", "cell_id": chosen,
        "x": emb[:, 0], "y": emb[:, 1],
        "cell_type": meta[args.group_col].astype(str).values,
        "tissue": meta[args.context_col].astype(str).values,
        "is_baseline": is_baseline.values,
    }))

    # --- L0 / L1 / L2 from the precomputed hierarchy ---
    for level in ["L0", "L1", "L2"]:
        level_cols = [c for c in coords.columns if c.startswith(level + "__")]
        if not level_cols:
            print(f"WARNING: no '{level}__' columns found — skipping {level}", file=sys.stderr)
            continue
        X = coords.loc[chosen, level_cols].values
        print(f"Running t-SNE on {level} ({X.shape})...", file=sys.stderr)
        emb = make_tsne(args.seed).fit_transform(X)
        rows.append(pd.DataFrame({
            "level": level, "cell_id": chosen,
            "x": emb[:, 0], "y": emb[:, 1],
            "cell_type": meta[args.group_col].astype(str).values,
            "tissue": meta[args.context_col].astype(str).values,
            "is_baseline": is_baseline.values,
        }))

    out_dir = Path(args.figure_data_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(out_dir / "sc_umap_overview.csv", index=False)
    print(f"Wrote sc_umap_overview.csv to {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
