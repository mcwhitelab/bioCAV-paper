#!/usr/bin/env python3
"""
cav_cluster.py — Hierarchical clustering of CAV direction vectors.

Builds a dendrogram over CAV directions and cuts it at multiple levels, writing
one cluster-ID column per level.  The linkage tree itself is the primary product
(leaf order drives heatmap row ordering), so three distance definitions are
available:

  --distance cosine     1 - cosine similarity on the raw CAV directions.
                        The original method.  Pair with --linkage average to
                        reproduce historical output, or --linkage ward.

  --distance diffusion  Walktrap's random-walk distance.  Builds a mutual-kNN
                        cosine graph, then computes the diffusion distance
                            r_ij = sqrt( sum_k (P^t_ik - P^t_jk)^2 / d(k) )
                        via the eigendecomposition of the normalised adjacency
                        (never materialising P^t).  Ward linkage on these
                        coordinates is exactly walktrap's agglomeration without
                        walktrap's adjacency restriction on merges — which buys
                        a tree with monotone heights that fcluster can cut by
                        distance.

  --distance snn        Shared-nearest-neighbour (Jaccard) reweighting of the
                        kNN graph.  Cheap denoising, but every non-neighbour
                        pair gets distance exactly 1.0, so leaf order among
                        unrelated CAVs is arbitrary.  Sanity check only.

Output CSV columns:
  accession | clan_acc | clan_name | clan_description |
  pfam_short_name | pfam_long_name | pfam_long_description |
  dendrogram_order | cut_* …

Side outputs (next to --out):
  <out>.linkage.npy   the linkage matrix Z, for cophenetic scoring
  <out>.names.txt     row names in linkage order (Z's leaf indices)
  <out>.cuts.csv      column -> cut height -> n_clusters

Usage
-----
python cav_cluster.py \\
    --cav-pattern      "/path/to/cavs/*/L25_concept_v1.npy" \\
    --pfam-annotations pfamA.txt \\
    --clan-annotations Pfam-A.clans.tsv \\
    --distance diffusion --linkage ward --cut-mode maxclust \\
    --out              results/pfam_clusters_diffusion.csv
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
from scipy.sparse import csr_matrix, diags, triu
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

_VERSION_SUFFIX = re.compile(r'_v\d+$')

# Original distance-criterion thresholds (calibrated for 1 - cosine, where the
# bulk of pairs sit near 0.97).  Meaningless on ward/diffusion height scales.
THRESHOLDS = sorted([round(t, 2) for t in np.arange(0.1, 1.0, 0.1)]
                    + [0.75, 0.85, 0.95])

# Cluster counts produced by the historical cosine+average+distance run on the
# 27k library.  Used by --cut-mode maxclust so every method is scored at
# matched resolutions.
DEFAULT_MAXCLUST = [314, 1326, 2769, 4629, 6756, 9172,
                    14214, 18984, 22788, 25278, 26640, 27272]

# Quantiles of the merge-height distribution, for --cut-mode quantile.
DEFAULT_QUANTILES = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40,
                     0.50, 0.65, 0.80, 0.90, 0.95, 0.99]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _concept_name(path: Path) -> str:
    parent = path.parent.name
    stem   = _VERSION_SUFFIX.sub('', path.stem)
    return parent if 'concept' in stem else stem


def load_directions(pattern: str):
    """Return (names, M) with M unit-normalised, float32, shape (n, d)."""
    paths = sorted(glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matched: {pattern}")
    names, vecs = [], []
    for p in paths:
        path = Path(p)
        v    = np.load(path).astype(np.float32).ravel()
        norm = np.linalg.norm(v)
        if norm <= 1e-10:
            logger.warning(f"Skipping zero-norm vector: {path}")
            continue
        names.append(_concept_name(path))
        vecs.append(v / norm)
    M = np.vstack(vecs)
    logger.info(f"Loaded {len(names)} CAV directions, dim {M.shape[1]}")
    return names, M


def load_pfam_annotations(pfam_txt: str) -> pd.DataFrame:
    """pfamA.txt — col 0 accession, 1 short_name, 3 long_name, 7 long_description."""
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
    """Pfam-A.clans.tsv — col 0 accession, 1 clan_acc, 2 clan_name, 4 clan_description."""
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
# kNN graph
# ---------------------------------------------------------------------------

def knn_graph(M, k, floor, mutual=True, block=2048):
    """
    Symmetric sparse cosine-similarity kNN graph.

    Returns (W, n_forced) where W is csr_matrix (n, n) with positive weights and
    n_forced is the number of nodes that would otherwise have been isolated and
    were rescued with their single nearest neighbour.
    """
    n = M.shape[0]
    rows, cols, vals = [], [], []
    top1 = np.empty((n, 2), dtype=np.float64)   # (neighbour index, similarity)

    for s in range(0, n, block):
        e = min(s + block, n)
        S = M[s:e] @ M.T                                   # (b, n) float32
        S[np.arange(e - s), np.arange(s, e)] = -np.inf     # mask self

        idx = np.argpartition(-S, kth=k, axis=1)[:, :k]
        v   = np.take_along_axis(S, idx, axis=1)

        best = np.argmax(v, axis=1)
        top1[s:e, 0] = idx[np.arange(e - s), best]
        top1[s:e, 1] = v[np.arange(e - s), best]

        keep = v >= floor
        r = np.repeat(np.arange(s, e), k).reshape(e - s, k)
        rows.append(r[keep]); cols.append(idx[keep]); vals.append(v[keep])

    rows = np.concatenate(rows); cols = np.concatenate(cols)
    vals = np.concatenate(vals).astype(np.float64)
    A = csr_matrix((vals, (rows, cols)), shape=(n, n))

    # Cosine is symmetric, so A[i,j] == A[j,i] whenever both directions were
    # kept.  With strictly positive weights, minimum() == "both present" and
    # maximum() == "either present".
    W = A.minimum(A.T) if mutual else A.maximum(A.T)
    W.eliminate_zeros()

    deg = np.asarray(W.sum(axis=1)).ravel()
    isolated = np.flatnonzero(deg == 0)
    if len(isolated):
        j = top1[isolated, 0].astype(np.int64)
        w = top1[isolated, 1]
        rescue = csr_matrix((np.concatenate([w, w]),
                             (np.concatenate([isolated, j]),
                              np.concatenate([j, isolated]))), shape=(n, n))
        W = W.maximum(rescue)
        W.eliminate_zeros()

    if W.data.min() <= 0:
        raise ValueError("kNN graph has non-positive weights; raise --knn-floor")

    nnz = W.nnz
    logger.info(f"kNN graph: k={k} floor={floor} mutual={mutual} — "
                f"{nnz} directed edges ({nnz / 2:.0f} undirected), "
                f"mean degree {nnz / n:.1f}, {len(isolated)} isolated nodes rescued")
    return W, len(isolated)


# ---------------------------------------------------------------------------
# Distance / embedding construction
# ---------------------------------------------------------------------------

def diffusion_coords(W, steps, n_eigs):
    """
    Coordinates whose Euclidean distance equals walktrap's random-walk distance
        r_ij = sqrt( sum_k (P^t_ik - P^t_jk)^2 / d(k) ),   P = D^-1 W

    P is similar to the symmetric S = D^-1/2 W D^-1/2, so the eigenvectors of S
    give P's spectral decomposition without ever forming P^t.
    """
    n = W.shape[0]
    d = np.asarray(W.sum(axis=1)).ravel()
    dinv_sqrt = 1.0 / np.sqrt(d)
    S = diags(dinv_sqrt) @ W @ diags(dinv_sqrt)

    n_comp, comp_labels = connected_components(W, directed=False)
    logger.info(f"Graph has {n_comp} connected component(s); "
                f"largest holds {np.bincount(comp_labels).max()} nodes")

    # Each connected component contributes an eigenvalue of exactly 1, and those
    # component-indicator eigenvectors sit at the top of the spectrum.  Asking
    # for only n_eigs total would spend the whole budget on them and leave
    # almost nothing for structure *within* the big component, so widen the
    # request by the component count.
    k = min(n_eigs + n_comp, n - 2)
    if n_comp > 1:
        logger.info(f"Widening eigenvector request to {k} so ~{n_eigs} non-trivial "
                    f"vectors survive the {n_comp} component indicators. "
                    f"To fragment the graph less, lower --knn-floor or use --plain-knn.")
    logger.info(f"Eigendecomposition: requesting {k} eigenpairs of a {n}x{n} sparse matrix …")
    evals, evecs = eigsh(S, k=k, which='LA')
    order = np.argsort(-evals)
    evals, evecs = evals[order], evecs[:, order]
    logger.info(f"Eigenvalues: max {evals[0]:.4f}, "
                f"#(>0.999) {(evals > 0.999).sum()} (component indicators), min {evals[-1]:.4f}")

    # Drop the global trivial eigenvector (psi_0 ∝ sqrt(d), which becomes a
    # constant coordinate after the 1/sqrt(d) scaling and contributes nothing
    # to distances).  Any remaining lambda≈1 vectors are component indicators
    # and are kept — they correctly place separate components far apart.
    evals, evecs = evals[1:], evecs[:, 1:]

    scale = np.sign(evals) * (np.abs(evals) ** steps)
    phi = (evecs * scale) * dinv_sqrt[:, None]

    kept = np.abs(scale) > 1e-12
    logger.info(f"Diffusion coords: t={steps}, {kept.sum()}/{len(scale)} "
                f"eigenvectors carry non-negligible weight after lambda^t")
    return np.ascontiguousarray(phi, dtype=np.float64)


def snn_distance(W):
    """Dense 1 - Jaccard(shared neighbours).  Non-neighbour pairs get exactly 1.0."""
    n = W.shape[0]
    A = (W > 0).astype(np.float64)
    inter = (A @ A.T).toarray()
    deg = np.asarray(A.sum(axis=1)).ravel()
    union = deg[:, None] + deg[None, :] - inter
    with np.errstate(divide='ignore', invalid='ignore'):
        jac = np.where(union > 0, inter / union, 0.0)
    dist = 1.0 - jac
    np.fill_diagonal(dist, 0.0)
    ties = np.isclose(dist, 1.0).sum() / (n * n)
    logger.info(f"SNN distance: {ties:.1%} of pairs are exactly 1.0 (tied)")
    return dist


def cosine_distance(M):
    logger.info("Computing dense cosine similarity matrix …")
    sim  = (M.astype(np.float64)) @ (M.astype(np.float64)).T
    dist = np.clip(1.0 - sim, 0, 2)
    np.fill_diagonal(dist, 0)
    return dist


# ---------------------------------------------------------------------------
# Cutting
# ---------------------------------------------------------------------------

def make_cuts(Z, n, cut_mode, cut_values):
    """
    Return (cuts, meta) where cuts maps column name -> 1-based label array in Z's
    row order, and meta is a list of (column, criterion, value, height, n_clusters).

    Columns are emitted coarse -> fine (fewest clusters leftmost).
    """
    heights = Z[:, 2]
    entries = []

    if cut_mode == "distance":
        for t in cut_values:
            labels = fcluster(Z, t=t, criterion="distance")
            entries.append((f"cut_{int(round(t * 100)):03d}", "distance", t, t, labels))

    elif cut_mode == "quantile":
        for q in cut_values:
            h = float(np.quantile(heights, q))
            labels = fcluster(Z, t=h, criterion="distance")
            entries.append((f"cut_q{int(round(q * 1000)):04d}", "quantile", q, h, labels))

    elif cut_mode == "maxclust":
        for kk in cut_values:
            kk = int(min(kk, n))
            labels = fcluster(Z, t=kk, criterion="maxclust")
            # Height at which this many clusters remain
            n_merges = n - len(np.unique(labels))
            h = float(heights[n_merges - 1]) if n_merges > 0 else 0.0
            entries.append((f"cut_k{kk:05d}", "maxclust", kk, h, labels))
    else:
        raise ValueError(f"Unknown cut mode: {cut_mode}")

    # Sort coarse -> fine
    entries.sort(key=lambda e: len(np.unique(e[4])))
    cuts, meta = {}, []
    for col, crit, val, h, labels in entries:
        cuts[col] = labels
        meta.append((col, crit, val, h, int(len(np.unique(labels)))))
    return cuts, meta


# ---------------------------------------------------------------------------
# Clustering driver
# ---------------------------------------------------------------------------

def cluster_cavs(names, M, args):
    n = len(names)

    if args.linkage_in:
        # Re-cut a tree that was already built.  The linkage matrix and its leaf
        # names fully determine every cut, so there is no need to reload the CAVs
        # or redo the eigendecomposition.
        Z = np.load(args.linkage_in)
        if Z.shape[0] != n - 1:
            raise ValueError(f"{args.linkage_in} has {Z.shape[0]} merges but "
                             f"{n} names were supplied")
        logger.info(f"Reusing linkage matrix from {args.linkage_in}")

    elif args.distance == "cosine":
        dist = cosine_distance(M)
        logger.info(f"Running {args.linkage}-linkage hierarchical clustering …")
        Z = linkage(squareform(dist, checks=False), method=args.linkage)
        del dist

    elif args.distance == "diffusion":
        W, _ = knn_graph(M, args.knn, args.knn_floor, mutual=not args.plain_knn)
        phi = diffusion_coords(W, args.walk_steps, args.n_eigs)
        logger.info(f"Running {args.linkage}-linkage on {phi.shape[1]}-d diffusion coords …")
        Z = linkage(phi, method=args.linkage)

    elif args.distance == "snn":
        W, _ = knn_graph(M, args.knn, args.knn_floor, mutual=not args.plain_knn)
        dist = snn_distance(W)
        logger.info(f"Running {args.linkage}-linkage on SNN distances …")
        Z = linkage(squareform(dist, checks=False), method=args.linkage)
        del dist
    else:
        raise ValueError(f"Unknown distance: {args.distance}")

    order = leaves_list(Z)
    names_ordered = [names[i] for i in order]
    logger.info(f"Leaf order computed for {n} CAVs; "
                f"merge heights span [{Z[:, 2].min():.4g}, {Z[:, 2].max():.4g}]")

    cut_mode = args.cut_mode
    if cut_mode == "auto":
        cut_mode = "distance" if (args.distance == "cosine"
                                  and args.linkage == "average") else "maxclust"
        logger.info(f"--cut-mode auto resolved to '{cut_mode}'")

    if args.cut_values:
        cut_values = [float(v) for v in args.cut_values]
        if cut_mode == "maxclust":
            cut_values = [int(v) for v in cut_values]
    elif cut_mode == "distance":
        cut_values = THRESHOLDS
    elif cut_mode == "quantile":
        cut_values = DEFAULT_QUANTILES
    else:
        cut_values = DEFAULT_MAXCLUST

    cuts, meta = make_cuts(Z, n, cut_mode, cut_values)

    rows = {name: {"dendrogram_order": i + 1} for i, name in enumerate(names_ordered)}
    rng = np.random.default_rng(42)
    pos = {name: i for i, name in enumerate(names)}

    for col in cuts:
        raw_ids = cuts[col]
        unique_ids = np.unique(raw_ids)
        if args.shuffle_ids:
            # Randomly permute labels for visual contrast in heatmaps.  Note this
            # destroys any correspondence between cluster IDs at different cut
            # levels — use --no-shuffle-ids when comparing across resolutions.
            mapped = rng.permutation(len(unique_ids)) + 1
        else:
            # Number in order of first appearance along the leaf order, so the
            # same nested group keeps a stable-ish ID across levels.
            first_seen = {}
            for name in names_ordered:
                cid = raw_ids[pos[name]]
                if cid not in first_seen:
                    first_seen[cid] = len(first_seen) + 1
            mapped = np.array([first_seen[c] for c in unique_ids])
        id_map = dict(zip(unique_ids, mapped))
        for name in names_ordered:
            rows[name][col] = int(id_map[raw_ids[pos[name]]])

    result = pd.DataFrame.from_dict(rows, orient="index")
    result.index.name = "accession"
    meta_df = pd.DataFrame(meta, columns=["column", "criterion", "value",
                                          "height", "n_clusters"])
    return result, Z, meta_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Hierarchical clustering of CAV directions → annotated CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--cav-pattern",
                        help="Glob pattern matching CAV .npy files, "
                             "e.g. '/path/to/cavs/*/L25_concept_v1.npy'. "
                             "Not needed with --linkage-in.")
    parser.add_argument("--linkage-in",
                        help="Reuse a saved linkage matrix (<out>.linkage.npy from "
                             "an earlier run) instead of rebuilding the tree. "
                             "Requires --names-in. Lets you re-cut an existing "
                             "dendrogram at new levels in seconds.")
    parser.add_argument("--names-in",
                        help="Leaf names for --linkage-in (<out>.names.txt).")
    parser.add_argument("--pfam-annotations",
                        help="pfamA.txt — short name, long name, extended description.")
    parser.add_argument("--clan-annotations",
                        help="Pfam-A.clans.tsv — clan accession, clan name, description.")
    parser.add_argument("--out", required=True, help="Output CSV path.")

    parser.add_argument("--distance", choices=["cosine", "diffusion", "snn"],
                        default="cosine",
                        help="Distance the dendrogram is built on (default: cosine).")
    parser.add_argument("--linkage", default="average",
                        choices=["average", "ward", "complete", "single"],
                        help="Linkage method (default: average). "
                             "ward is recommended for diffusion.")

    g = parser.add_argument_group("kNN graph (diffusion / snn)")
    g.add_argument("--knn", type=int, default=15,
                   help="Neighbours per node (default: 15).")
    g.add_argument("--knn-floor", type=float, default=0.15,
                   help="Drop kNN edges below this cosine similarity (default: 0.15).")
    g.add_argument("--plain-knn", action="store_true",
                   help="Symmetrise by union instead of mutual-kNN intersection.")
    g.add_argument("--walk-steps", type=int, default=4,
                   help="Random-walk length t for diffusion distance (default: 4).")
    g.add_argument("--n-eigs", type=int, default=200,
                   help="Eigenvectors for the diffusion embedding (default: 200).")

    c = parser.add_argument_group("cutting")
    c.add_argument("--cut-mode", default="auto",
                   choices=["auto", "distance", "quantile", "maxclust"],
                   help="How to cut the tree. 'distance' uses fixed heights "
                        "(only meaningful for cosine); 'quantile' uses quantiles "
                        "of the merge-height distribution; 'maxclust' targets "
                        "fixed cluster counts so methods are comparable. "
                        "'auto' picks distance for cosine+average, else maxclust.")
    c.add_argument("--cut-values", nargs="+",
                   help="Explicit cut values (heights, quantiles, or counts).")
    c.add_argument("--no-shuffle-ids", dest="shuffle_ids", action="store_false",
                   help="Number clusters by first appearance in leaf order instead "
                        "of randomly permuting them.")
    parser.set_defaults(shuffle_ids=True)

    args = parser.parse_args()

    if args.linkage_in:
        if not args.names_in:
            parser.error("--linkage-in requires --names-in")
        names = Path(args.names_in).read_text().split()
        M = None
        logger.info(f"Loaded {len(names)} leaf names from {args.names_in}")
    elif args.cav_pattern:
        names, M = load_directions(args.cav_pattern)
    else:
        parser.error("need --cav-pattern, or --linkage-in with --names-in")
    cluster_df, Z, meta_df = cluster_cavs(names, M, args)

    pfam  = load_pfam_annotations(args.pfam_annotations) if args.pfam_annotations else None
    clans = load_clan_annotations(args.clan_annotations) if args.clan_annotations else None

    annot_cols = {}
    for acc in cluster_df.index:
        row = {}
        if clans is not None and acc in clans.index:
            r = clans.loc[acc]
            row["clan_acc"]         = r["clan_acc"]
            row["clan_name"]        = r["clan_name"]
            row["clan_description"] = r["clan_description"]
        else:
            row["clan_acc"] = row["clan_name"] = row["clan_description"] = ""
        if pfam is not None and acc in pfam.index:
            r = pfam.loc[acc]
            row["pfam_short_name"]       = r["pfam_short_name"]
            row["pfam_long_name"]        = r["pfam_long_name"]
            row["pfam_long_description"] = r["pfam_long_description"]
        else:
            row["pfam_short_name"] = row["pfam_long_name"] = row["pfam_long_description"] = ""
        annot_cols[acc] = row

    annot_df = pd.DataFrame.from_dict(annot_cols, orient="index")
    annot_df.index.name = "accession"

    out_df = pd.concat([annot_df, cluster_df], axis=1)
    out_df.index.name = "accession"

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path)

    np.save(out_path.with_suffix(".linkage.npy"), Z)
    out_path.with_suffix(".names.txt").write_text("\n".join(names) + "\n")
    meta_df.to_csv(out_path.with_suffix(".cuts.csv"), index=False)

    logger.info(f"Saved {len(out_df)} rows × {len(out_df.columns)} columns → {out_path}")
    logger.info(f"Wrote linkage matrix, leaf names, and cut metadata alongside.")
    for _, r in meta_df.iterrows():
        logger.info(f"  {r['column']}: {r['n_clusters']} clusters "
                    f"(height {r['height']:.4g})")


if __name__ == "__main__":
    main()
