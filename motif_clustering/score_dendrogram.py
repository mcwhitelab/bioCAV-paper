#!/usr/bin/env python3
"""
score_dendrogram.py — Evaluate a CAV dendrogram against Pfam clans.

Reports three families of metric, because "a good dendrogram" means different
things depending on what the tree is used for:

  1. Clan contiguity in leaf order.  For each clan, the span of its members'
     positions in the dendrogram leaf order, divided by clan size (1.0 = the
     clan's members are perfectly adjacent).  This scores the heatmap
     row-ordering use of the tree, which no cut-based metric measures.  A
     permutation baseline is reported alongside, since a large clan is
     contiguous by chance more easily than a small one.

  2. Cut quality vs clans.  Per cut column: ARI, AMI, and median per-clan
     best-match F1, computed on the subset of CAVs that actually have a clan.
     ARI is pair-based and clan sizes here span 5 to 699, so a single ARI is
     dominated by the handful of giant clans; the median per-clan F1 weights
     every clan equally and is the number to check for small-clan recovery.
     Read the curve across cut levels, not just its maximum — a good tree
     scores well over a broad plateau of resolutions.

  3. Cophenetic correlation against the raw cosine distances, estimated on a
     random sample of pairs (the full condensed matrix is 377M entries at 27k).
     Note this measures fidelity to the *input* distances; when the input is
     noise-dominated, higher is not automatically better.

Usage
-----
python score_dendrogram.py \\
    --clusters results/figures/pfam_clusters_diffusion.csv \\
    --cav-pattern "/path/to/cavs/*/L25_concept_v1.npy" \\
    --out results/figures/scores_diffusion.csv
"""

import re
import argparse
import logging
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.metrics.cluster import contingency_matrix

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

_VERSION_SUFFIX = re.compile(r'_v\d+$')


# ---------------------------------------------------------------------------
# 1. Leaf-order contiguity
# ---------------------------------------------------------------------------

def contiguity(order, clan, min_size=2, n_perm=20, seed=0):
    """
    order : 1-based dendrogram position per row
    clan  : clan accession per row ('' when unassigned)

    Returns (per_clan_df, summary_dict).  Normalised span of 1.0 means the
    clan occupies a contiguous block of the leaf order.
    """
    n = len(order)
    df = pd.DataFrame({"order": np.asarray(order), "clan": np.asarray(clan)})
    df = df[df["clan"].astype(str).str.strip() != ""]

    g = df.groupby("clan")["order"]
    stats = pd.DataFrame({"size": g.size(), "lo": g.min(), "hi": g.max()})
    stats = stats[stats["size"] >= min_size]
    stats["span"] = stats["hi"] - stats["lo"] + 1
    stats["norm_span"] = stats["span"] / stats["size"]
    stats["span_frac"] = stats["span"] / n

    rng = np.random.default_rng(seed)
    perm_medians = []
    vals = df["order"].to_numpy()
    for _ in range(n_perm):
        # Null: keep the observed set of leaf positions and each clan's size,
        # but scramble which clan sits at which position.
        tmp = pd.DataFrame({"order": rng.permutation(vals), "clan": df["clan"].to_numpy()})
        gg = tmp.groupby("clan")["order"]
        s = pd.DataFrame({"size": gg.size(), "lo": gg.min(), "hi": gg.max()})
        s = s[s["size"] >= min_size]
        perm_medians.append(float(np.median((s["hi"] - s["lo"] + 1) / s["size"])))

    summary = {
        "n_clans_scored":        int(len(stats)),
        "median_norm_span":      float(stats["norm_span"].median()),
        "mean_norm_span":        float(stats["norm_span"].mean()),
        "median_span_frac":      float(stats["span_frac"].median()),
        "frac_clans_norm_span_lt2": float((stats["norm_span"] < 2).mean()),
        "perm_median_norm_span": float(np.mean(perm_medians)),
    }
    summary["contiguity_ratio"] = (summary["perm_median_norm_span"]
                                   / summary["median_norm_span"])
    return stats.sort_values("norm_span"), summary


# ---------------------------------------------------------------------------
# 2. Cut quality vs clans
# ---------------------------------------------------------------------------

def median_per_clan_f1(true_labels, pred_labels):
    """Median over clans of the best-matching cluster's F1. Weights clans equally."""
    C = contingency_matrix(true_labels, pred_labels, sparse=True).tocsr().astype(np.float64)
    clan_sizes = np.asarray(C.sum(axis=1)).ravel()
    clus_sizes = np.asarray(C.sum(axis=0)).ravel()

    best = np.zeros(C.shape[0])
    for i in range(C.shape[0]):
        s, e = C.indptr[i], C.indptr[i + 1]
        if s == e:
            continue
        inter = C.data[s:e]
        f1 = 2.0 * inter / (clan_sizes[i] + clus_sizes[C.indices[s:e]])
        best[i] = f1.max()
    return float(np.median(best)), float(best.mean())


def score_cuts(df, cut_cols):
    labelled = df[df["clan_acc"] != ""]
    true_labels = labelled["clan_acc"].to_numpy()
    logger.info(f"Scoring cuts on {len(labelled)} clan-labelled rows "
                f"({len(labelled)/len(df):.1%} of {len(df)}), "
                f"{pd.Series(true_labels).nunique()} clans")

    rows = []
    for col in cut_cols:
        pred = labelled[col].to_numpy()
        med_f1, mean_f1 = median_per_clan_f1(true_labels, pred)
        rows.append({
            "column": col,
            "n_clusters_all": int(df[col].nunique()),
            "n_clusters_labelled": int(pd.Series(pred).nunique()),
            "ARI": adjusted_rand_score(true_labels, pred),
            "AMI": adjusted_mutual_info_score(true_labels, pred),
            "median_per_clan_F1": med_f1,
            "mean_per_clan_F1": mean_f1,
        })
    return pd.DataFrame(rows).sort_values("n_clusters_all")


# ---------------------------------------------------------------------------
# 3. Sampled cophenetic correlation
# ---------------------------------------------------------------------------

def sampled_cophenetic(Z, M, n_sample=2500, seed=0):
    """
    Cophenetic distance for a random sample of observation pairs, without
    materialising the full condensed matrix.

    Walks the merges in order maintaining, for each live cluster, the sampled
    observations it contains.  When two clusters merge, every cross pair gets
    that merge's height as its cophenetic distance.
    """
    n = Z.shape[0] + 1
    rng = np.random.default_rng(seed)
    sample = np.sort(rng.choice(n, size=min(n_sample, n), replace=False))

    memb = {int(i): np.array([i], dtype=np.int64) for i in sample}
    pi, pj, ph = [], [], []

    for m in range(n - 1):
        a, b, h = int(Z[m, 0]), int(Z[m, 1]), float(Z[m, 2])
        A = memb.pop(a, None)
        B = memb.pop(b, None)
        if A is not None and B is not None:
            pi.append(np.repeat(A, len(B)))
            pj.append(np.tile(B, len(A)))
            ph.append(np.full(len(A) * len(B), h))
            memb[n + m] = np.concatenate([A, B])
        elif A is not None:
            memb[n + m] = A
        elif B is not None:
            memb[n + m] = B

    pi = np.concatenate(pi); pj = np.concatenate(pj); ph = np.concatenate(ph)
    logger.info(f"Cophenetic: {len(ph)} sampled pairs from {len(sample)} observations")

    # Input distance for the same pairs (1 - cosine on unit-normed directions)
    sims = np.einsum('ij,ij->i', M[pi], M[pj]).astype(np.float64)
    dist = np.clip(1.0 - sims, 0, 2)

    pearson  = float(np.corrcoef(ph, dist)[0, 1])
    spearman = float(np.corrcoef(pd.Series(ph).rank(), pd.Series(dist).rank())[0, 1])
    return {"cophenetic_vs_rawcosine_pearson": pearson,
            "cophenetic_vs_rawcosine_spearman": spearman,
            "n_pairs": int(len(ph))}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _concept_name(path: Path) -> str:
    parent = path.parent.name
    stem   = _VERSION_SUFFIX.sub('', path.stem)
    return parent if 'concept' in stem else stem


def load_directions_for(names, pattern):
    """Load unit-normed directions in the given name order."""
    by_name = {}
    for p in sorted(glob(pattern)):
        by_name[_concept_name(Path(p))] = p
    missing = [n for n in names if n not in by_name]
    if missing:
        raise KeyError(f"{len(missing)} names not found in --cav-pattern, e.g. {missing[:3]}")
    M = np.vstack([np.load(by_name[n]).astype(np.float32).ravel() for n in names])
    M /= np.linalg.norm(M, axis=1, keepdims=True)
    return M


def main():
    ap = argparse.ArgumentParser(
        description="Score a CAV dendrogram against Pfam clans.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--clusters", required=True,
                    help="Cluster CSV from cav_cluster.py.")
    ap.add_argument("--cav-pattern",
                    help="CAV glob, for cophenetic correlation. Needs the "
                         "<clusters>.linkage.npy / .names.txt sidecars.")
    ap.add_argument("--out", help="Write per-cut scores here (CSV).")
    ap.add_argument("--per-clan-out", help="Write per-clan contiguity here (CSV).")
    ap.add_argument("--n-sample", type=int, default=2500,
                    help="Observations sampled for cophenetic correlation.")
    ap.add_argument("--label", default=None, help="Name for this run in the output.")
    args = ap.parse_args()

    path = Path(args.clusters)
    label = args.label or path.stem
    df = pd.read_csv(path, dtype=str)
    for c in ("clan_acc", "clan_name"):
        if c in df.columns:
            df[c] = df[c].fillna("").astype(str).str.strip()
    cut_cols = [c for c in df.columns if c.startswith("cut_")]
    for c in cut_cols + ["dendrogram_order"]:
        df[c] = df[c].astype(int)
    logger.info(f"Loaded {len(df)} rows, {len(cut_cols)} cut columns from {path}")

    print(f"\n=== {label} ===")

    # --- 1. contiguity
    per_clan, summary = contiguity(df["dendrogram_order"].to_numpy(),
                                   df["clan_acc"].to_numpy())
    print("\n[1] Clan contiguity in leaf order  (norm_span 1.0 = perfectly adjacent)")
    for k, v in summary.items():
        print(f"    {k:28s} {v:.4g}")
    print("    -- 10 most contiguous clans --")
    print(per_clan.head(10).to_string())

    # --- 2. cut quality
    scores = score_cuts(df, cut_cols)
    print("\n[2] Cut quality vs Pfam clans (clan-labelled rows only)")
    print(scores.to_string(index=False,
                           float_format=lambda x: f"{x:.4f}"))
    best_ari = scores.loc[scores["ARI"].idxmax()]
    best_f1  = scores.loc[scores["median_per_clan_F1"].idxmax()]
    print(f"\n    best ARI  {best_ari['ARI']:.4f} at {best_ari['column']} "
          f"({best_ari['n_clusters_all']} clusters)")
    print(f"    best F1   {best_f1['median_per_clan_F1']:.4f} at {best_f1['column']} "
          f"({best_f1['n_clusters_all']} clusters)")

    # --- 3. cophenetic
    coph = None
    zpath = path.with_suffix(".linkage.npy")
    npath = path.with_suffix(".names.txt")
    if args.cav_pattern and zpath.exists() and npath.exists():
        Z = np.load(zpath)
        names = npath.read_text().split()
        M = load_directions_for(names, args.cav_pattern)
        coph = sampled_cophenetic(Z, M, n_sample=args.n_sample)
        print("\n[3] Cophenetic correlation vs RAW COSINE distance\n    (the tree's own input only for --distance cosine; for diffusion/snn\n     this measures how far the reweighting moved from raw cosine)")
        for k, v in coph.items():
            print(f"    {k:28s} {v:.4g}")
    else:
        print("\n[3] Cophenetic correlation skipped "
              "(need --cav-pattern plus .linkage.npy / .names.txt sidecars)")

    if args.out:
        out = scores.copy()
        out.insert(0, "run", label)
        for k, v in summary.items():
            out[k] = v
        if coph:
            for k, v in coph.items():
                out[k] = v
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out, index=False)
        logger.info(f"Wrote scores → {args.out}")
    if args.per_clan_out:
        per_clan.to_csv(args.per_clan_out)
        logger.info(f"Wrote per-clan contiguity → {args.per_clan_out}")


if __name__ == "__main__":
    main()
