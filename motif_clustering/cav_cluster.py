#!/usr/bin/env python3
"""
cav_cluster.py — Hierarchical clustering of CAV direction vectors.

Computes all-by-all cosine similarity (= dot product on unit-normalised vectors),
clusters with average linkage on (1 - similarity), orders rows by dendrogram leaf
order, then cuts the tree at thresholds 0.1 … 0.9 and writes one cluster-ID column
per threshold.

Output CSV columns:
  accession | clan_acc | clan_name | clan_description |
  pfam_short_name | pfam_long_name | pfam_long_description |
  cut_01 | cut_02 | … | cut_09

Usage
-----
python specific_scripts/cav_cluster.py \\
    --cav-pattern      "/path/to/cavs/*/L25_concept_v1.npy" \\
    --pfam-annotations pfamA.txt \\
    --clan-annotations Pfam-A.clans.tsv \\
    --out              results/pfam_clusters.csv
"""

import re
import argparse
import logging
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

_VERSION_SUFFIX = re.compile(r'_v\d+$')

# Sorted ascending; columns will be emitted in reverse (coarse → fine, left → right)
THRESHOLDS = sorted([round(t, 2) for t in np.arange(0.1, 1.0, 0.1)]
                    + [0.75, 0.85, 0.95])


# ---------------------------------------------------------------------------
# Data loading  (mirrors cav_viz_pattern.py)
# ---------------------------------------------------------------------------

def _concept_name(path: Path) -> str:
    parent = path.parent.name
    stem   = _VERSION_SUFFIX.sub('', path.stem)
    return parent if 'concept' in stem else stem


def load_directions(pattern: str) -> dict:
    paths = sorted(glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matched: {pattern}")
    result = {}
    for p in paths:
        path = Path(p)
        name = _concept_name(path)
        v    = np.load(path).astype(np.float64)
        norm = np.linalg.norm(v)
        if norm > 1e-10:
            result[name] = v / norm
        else:
            logger.warning(f"Skipping zero-norm vector: {path}")
    logger.info(f"Loaded {len(result)} CAV directions")
    return result


def load_pfam_annotations(pfam_txt: str) -> pd.DataFrame:
    """
    pfamA.txt — tab-separated, no header.
    col 0: accession, col 1: short_name, col 3: long_name, col 7: long_description.
    """
    rows = []
    with open(pfam_txt) as fh:
        for line in fh:
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            rows.append({
                "accession":        parts[0].strip(),
                "pfam_short_name":  parts[1].strip(),
                "pfam_long_name":   parts[3].strip() if len(parts) > 3 else "",
                "pfam_long_description": parts[7].strip() if len(parts) > 7 else "",
            })
    df = pd.DataFrame(rows).set_index("accession")
    logger.info(f"Loaded {len(df)} PFAM annotations")
    return df


def load_clan_annotations(clan_txt: str) -> pd.DataFrame:
    """
    Pfam-A.clans.tsv — tab-separated, no header.
    col 0: accession, col 1: clan_acc, col 2: clan_name,
    col 3: pfam_short_name (redundant), col 4: clan_description.
    """
    rows = []
    with open(clan_txt) as fh:
        for line in fh:
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            rows.append({
                "accession":       parts[0].strip(),
                "clan_acc":        parts[1].strip(),
                "clan_name":       parts[2].strip(),
                "clan_description": parts[4].strip() if len(parts) > 4 else "",
            })
    df = pd.DataFrame(rows).set_index("accession")
    logger.info(f"Loaded {len(df)} clan mappings ({df['clan_name'].nunique()} clans)")
    return df


# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------

def cluster_cavs(cav_dirs: dict) -> pd.DataFrame:
    """
    1. Build cosine-similarity matrix (dot product — vectors already unit-normed).
    2. Convert to distance: d = 1 - similarity, clip to [0, 2].
    3. Average-linkage hierarchical clustering.
    4. Order rows by dendrogram leaf order.
    5. Cut at each threshold in THRESHOLDS and add cluster-ID columns.

    Cluster IDs at each threshold are integers numbered in order of first
    appearance along the leaf-ordered rows (so cluster 1 always contains the
    first row, cluster 2 the next distinct group, etc.).

    Returns a DataFrame indexed by accession with columns cut_01 … cut_09.
    """
    names  = list(cav_dirs.keys())
    matrix = np.vstack([cav_dirs[n] for n in names])

    logger.info("Computing cosine similarity matrix …")
    sim  = matrix @ matrix.T
    dist = np.clip(1.0 - sim, 0, 2)
    np.fill_diagonal(dist, 0)

    logger.info("Running average-linkage hierarchical clustering …")
    condensed = squareform(dist, checks=False)
    Z         = linkage(condensed, method="average")

    # Leaf order from dendrogram
    order       = leaves_list(Z)
    names_ordered = [names[i] for i in order]
    logger.info(f"Leaf order computed for {len(names_ordered)} CAVs")

    rows = {name: {"dendrogram_order": i + 1} for i, name in enumerate(names_ordered)}

    rng = np.random.default_rng(42)

    # Build columns finest-first so we can reverse at the end for coarse→fine order
    cut_cols = {}
    for t in THRESHOLDS:
        col     = f"cut_{int(round(t * 100)):03d}"
        raw_ids = fcluster(Z, t=t, criterion="distance")   # 1-based scipy labels

        # Randomly shuffle the cluster ID labels for visual contrast
        unique_ids = np.unique(raw_ids)
        shuffled   = rng.permutation(len(unique_ids)) + 1   # still 1-based
        id_map     = dict(zip(unique_ids, shuffled))

        col_vals = {}
        for name in names_ordered:
            col_vals[name] = int(id_map[raw_ids[names.index(name)]])
        cut_cols[col] = col_vals

    # Reverse so coarsest threshold (fewest clusters) is leftmost column
    for col in reversed(list(cut_cols)):
        for name in names_ordered:
            rows[name][col] = cut_cols[col][name]

    result = pd.DataFrame.from_dict(rows, orient="index")
    result.index.name = "accession"
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Hierarchical clustering of CAV directions → annotated CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--cav-pattern", required=True,
                        help="Glob pattern matching CAV .npy files, "
                             "e.g. '/path/to/cavs/*/L25_concept_v1.npy'.")
    parser.add_argument("--pfam-annotations",
                        help="pfamA.txt — short name, long name, extended description.")
    parser.add_argument("--clan-annotations",
                        help="Pfam-A.clans.tsv — clan accession, clan name, description.")
    parser.add_argument("--out", required=True,
                        help="Output CSV path.")
    args = parser.parse_args()

    # Load CAVs and cluster
    cav_dirs   = load_directions(args.cav_pattern)
    cluster_df = cluster_cavs(cav_dirs)   # index = accession, ordered by dendrogram

    # Load optional annotation tables
    pfam  = load_pfam_annotations(args.pfam_annotations) if args.pfam_annotations else None
    clans = load_clan_annotations(args.clan_annotations) if args.clan_annotations else None

    # Build annotation block (left columns), preserving dendrogram row order
    annot_cols = {}
    for acc in cluster_df.index:
        row = {}
        if clans is not None and acc in clans.index:
            r = clans.loc[acc]
            row["clan_acc"]         = r["clan_acc"]
            row["clan_name"]        = r["clan_name"]
            row["clan_description"] = r["clan_description"]
        else:
            row["clan_acc"]         = ""
            row["clan_name"]        = ""
            row["clan_description"] = ""
        if pfam is not None and acc in pfam.index:
            r = pfam.loc[acc]
            row["pfam_short_name"]       = r["pfam_short_name"]
            row["pfam_long_name"]        = r["pfam_long_name"]
            row["pfam_long_description"] = r["pfam_long_description"]
        else:
            row["pfam_short_name"]       = ""
            row["pfam_long_name"]        = ""
            row["pfam_long_description"] = ""
        annot_cols[acc] = row

    annot_df = pd.DataFrame.from_dict(annot_cols, orient="index")
    annot_df.index.name = "accession"

    # Combine: annotation columns first, then cut columns
    out_df = pd.concat([annot_df, cluster_df], axis=1)
    out_df.index.name = "accession"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out)
    logger.info(f"Saved {len(out_df)} rows × {len(out_df.columns)} columns → {args.out}")

    # Quick summary: cluster counts at each threshold
    for col in cluster_df.columns:
        n = cluster_df[col].nunique()
        logger.info(f"  {col}: {n} clusters")


if __name__ == "__main__":
    main()
