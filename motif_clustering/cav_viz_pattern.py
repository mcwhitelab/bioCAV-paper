#!/usr/bin/env python3
"""
cav_viz_pattern.py — Visualize CAV directions matched by a glob pattern.

Like cav_viz.py but instead of a library directory with the cavs/CONCEPT/concept_v1.npy
layout, this script accepts a shell glob pattern that directly matches .npy files,
e.g. ``cavs/PFAM/*/L25_concept_v1.npy``.

Each matched .npy file is one CAV; its concept name is taken from the parent
directory (e.g. cavs/PF01112/L25_concept_v1.npy → 'PF01112').

The CAV direction vectors themselves are embedded in 2D — no cell coordinates needed.

Usage
-----
# Static PNG:
python specific_scripts/cav_viz_pattern.py \\
    --cav-pattern "/path/to/cavs/*/L25_concept_v1.npy" \\
    --reducer     umap \\
    --out         results/figures/pfam_umap.png

# Interactive Plotly HTML (with PFAM hover annotations):
python specific_scripts/cav_viz_pattern.py \\
    --cav-pattern    "/path/to/cavs/*/L25_concept_v1.npy" \\
    --reducer        umap \\
    --pfam-annotations pfamA.txt \\
    --interactive \\
    --out            results/figures/pfam_umap.html
"""

import re
import argparse
import logging
from glob import glob
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

_VERSION_SUFFIX = re.compile(r'_v\d+$')


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _concept_name(path: Path) -> str:
    """
    Derive a concept name from a .npy file path.
    Uses the parent directory name when the file lives in a per-concept subdirectory
    (e.g. cavs/PF01112/L25_concept_v1.npy → 'PF01112'), otherwise falls back to
    the file stem with any trailing _v1-style suffix stripped.
    """
    parent = path.parent.name
    stem   = _VERSION_SUFFIX.sub('', path.stem)
    if 'concept' in stem:
        return parent
    return stem


def load_directions_from_pattern(pattern: str) -> Dict[str, np.ndarray]:
    """Load and unit-normalise CAV direction vectors from all files matching `pattern`."""
    paths = sorted(glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matched pattern: {pattern}")
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
    logger.info(f"Loaded {len(result)} CAV directions from pattern '{pattern}'")
    return result


def load_pfam_annotations(pfam_txt: str) -> pd.DataFrame:
    """
    Parse pfamA.txt (tab-separated, no header).
    Returns a DataFrame indexed by accession with columns:
      short_name, description, long_description.
    col 0: accession, col 1: short name, col 3: full English name,
    col 7: extended paragraph description.
    """
    rows = []
    with open(pfam_txt) as fh:
        for line in fh:
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            rows.append({
                "accession":        parts[0].strip(),
                "short_name":       parts[1].strip(),
                "long_name":        parts[3].strip() if len(parts) > 3 else "",
                "long_description": parts[7].strip() if len(parts) > 7 else "",
            })
    df = pd.DataFrame(rows).set_index("accession")
    logger.info(f"Loaded {len(df)} PFAM annotations from {pfam_txt}")
    return df


def load_clan_annotations(clan_txt: str) -> pd.DataFrame:
    """
    Parse a clan reference file (tab-separated, no header).
    col 0: PF accession, col 1: clan accession, col 2: clan name.
    Returns a DataFrame indexed by PF accession with columns: clan_acc, clan_name.
    """
    rows = []
    with open(clan_txt) as fh:
        for line in fh:
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            rows.append({
                "accession": parts[0].strip(),
                "clan_acc":  parts[1].strip(),
                "clan_name": parts[2].strip(),
            })
    df = pd.DataFrame(rows).set_index("accession")
    logger.info(f"Loaded {len(df)} clan mappings from {clan_txt} "
                f"({df['clan_name'].nunique()} unique clans)")
    return df


_CLAN_PALETTE = [
    "#808080", "#2f4f4f", "#556b2f", "#a0522d", "#228b22", "#7f0000",
    "#808000", "#483d8b", "#3cb371", "#008b8b", "#4682b4", "#000080",
    "#9acd32", "#daa520", "#8fbc8f", "#8b008b", "#b03060", "#ff0000",
    "#00ced1", "#ff8c00", "#ffd700", "#00ff00", "#8a2be2", "#00ff7f",
    "#dc143c", "#f4a460", "#9370db", "#0000ff", "#f08080", "#adff2f",
    "#da70d6", "#ff00ff", "#1e90ff", "#f0e68c", "#dda0dd", "#ff1493",
    "#afeeee", "#98fb98", "#87cefa", "#7fffd4", "#ffe4c4", "#ffb6c1",
]


def _clan_colors(clan_names: list) -> dict:
    """
    Assign a color from the fixed palette to each clan (sorted alphabetically).
    Cycles through the palette if there are more clans than colors.
    Unassigned / no-clan entries should use '#aaaaaa' (caller's responsibility).
    """
    unique = sorted(set(clan_names))
    return {name: _CLAN_PALETTE[i % len(_CLAN_PALETTE)]
            for i, name in enumerate(unique)}


# ---------------------------------------------------------------------------
# Shared: compute 2-D embedding
# ---------------------------------------------------------------------------

def _make_reducer(reducer: str, n_components: int = 2, metric: str = "euclidean",
                  n_neighbors: int = 30, min_dist: float = 0.1):
    if reducer == "umap":
        try:
            from umap import UMAP
        except ImportError:
            raise ImportError(
                "umap-learn is not installed. Use --reducer tsne or --reducer pca, "
                "or install with: pip install umap-learn"
            )
        return UMAP(n_components=n_components, random_state=42,
                    n_neighbors=n_neighbors, min_dist=min_dist, metric=metric)
    elif reducer == "tsne":
        from sklearn.manifold import TSNE
        return TSNE(n_components=n_components, random_state=42, metric=metric,
                    perplexity=min(n_neighbors, 50), max_iter=1000,
                    init="random", learning_rate="auto")
    elif reducer == "pca":
        return PCA(n_components=n_components, random_state=42)
    else:
        raise ValueError(f"Unknown reducer '{reducer}'. Choose: umap, tsne, pca")


def _cosine_dist_matrix(matrix: np.ndarray) -> np.ndarray:
    """Pairwise cosine distance for unit-normalised rows: 1 - (A @ A.T)."""
    sim  = matrix @ matrix.T
    dist = np.clip(1.0 - sim, 0, 2).astype(np.float64)
    np.fill_diagonal(dist, 0)
    return dist


def _cophenetic_dist_matrix(matrix: np.ndarray) -> np.ndarray:
    """
    Build the cophenetic distance matrix from average-linkage hierarchical
    clustering on cosine distances.  The cophenetic distance between two
    points is the height (distance threshold) at which they first merge in
    the dendrogram — it encodes the full tree structure and resolves clusters
    that look equidistant in raw cosine space.
    """
    from scipy.cluster.hierarchy import linkage, cophenet
    from scipy.spatial.distance import squareform

    logger.info("Building cosine distance matrix for cophenetic embedding …")
    dist_sq   = _cosine_dist_matrix(matrix)
    condensed = squareform(dist_sq, checks=False)

    logger.info("Running average-linkage clustering …")
    Z = linkage(condensed, method="average")

    logger.info("Extracting cophenetic distances …")
    _, coph_condensed = cophenet(Z, condensed)   # correlation coeff + condensed dists
    coph_square = squareform(coph_condensed)
    np.fill_diagonal(coph_square, 0)

    corr = float(_)
    logger.info(f"Cophenetic correlation with cosine distances: {corr:.4f}")
    return coph_square


def _embed(cav_dirs: Dict[str, np.ndarray], reducer: str, dims: int = 2,
           distance_mode: str = "approximate",
           n_neighbors: Optional[int] = None, min_dist: Optional[float] = None):
    """
    distance_mode:
      'approximate'  — UMAP/t-SNE on raw high-D vectors (fast, default)
      'precomputed'  — exact cosine distance matrix passed as precomputed
      'cophenetic'   — cophenetic distances from hierarchical clustering
                       passed as precomputed; best resolves dendrogram structure

    n_neighbors / min_dist override UMAP defaults. Cophenetic mode uses
    tighter defaults (n_neighbors=10, min_dist=0.0) unless overridden.
    """
    names  = list(cav_dirs.keys())
    matrix = np.vstack([cav_dirs[n] for n in names])

    # Sensible defaults per mode
    if distance_mode == "cophenetic":
        nn  = n_neighbors if n_neighbors is not None else 10
        md  = min_dist    if min_dist    is not None else 0.0
    else:
        nn  = n_neighbors if n_neighbors is not None else 30
        md  = min_dist    if min_dist    is not None else 0.1

    if distance_mode == "approximate":
        logger.info(f"Running {reducer.upper()} ({dims}D) on {len(names)} CAV vectors "
                    f"(n_neighbors={nn}, min_dist={md}) …")
        red    = _make_reducer(reducer, n_components=dims, n_neighbors=nn, min_dist=md)
        coords = red.fit_transform(matrix)

    elif distance_mode in ("precomputed", "cophenetic"):
        if reducer == "pca":
            raise ValueError(f"--distance-mode {distance_mode} is not supported with pca.")
        if distance_mode == "cophenetic":
            dist = _cophenetic_dist_matrix(matrix)
        else:
            logger.info(f"Building {len(names)}×{len(names)} cosine distance matrix …")
            dist = _cosine_dist_matrix(matrix)
        logger.info(f"Running {reducer.upper()} ({dims}D) on {distance_mode} distances "
                    f"(n_neighbors={nn}, min_dist={md}) …")
        red    = _make_reducer(reducer, n_components=dims, metric="precomputed",
                               n_neighbors=nn, min_dist=md)
        coords = red.fit_transform(dist)

    else:
        raise ValueError(f"Unknown distance_mode '{distance_mode}'. "
                         "Choose: approximate, precomputed, cophenetic")

    cols     = ["D1", "D2", "D3"][:dims]
    coord_df = pd.DataFrame(coords, columns=cols, index=names)
    return coord_df, red


def _axis_labels(reducer: str, red) -> Tuple[str, str]:
    xl, yl = {"umap": ("UMAP 1", "UMAP 2"),
               "tsne": ("t-SNE 1", "t-SNE 2"),
               "pca":  ("PC 1", "PC 2")}.get(reducer, ("Dim 1", "Dim 2"))
    if reducer == "pca" and hasattr(red, "explained_variance_ratio_"):
        var = red.explained_variance_ratio_ * 100
        xl, yl = f"PC 1 ({var[0]:.1f}%)", f"PC 2 ({var[1]:.1f}%)"
    return xl, yl


# ---------------------------------------------------------------------------
# Static plot (matplotlib)
# ---------------------------------------------------------------------------

def plot_direction_map(
    cav_pattern: str,
    out_path: str,
    reducer: str = "pca",
    figsize: Tuple[int, int] = (10, 9),
    label_points: bool = True,
):
    """Embed all matched CAV direction vectors in 2D and save a static PNG."""
    cav_dirs = load_directions_from_pattern(cav_pattern)
    coord_df, red = _embed(cav_dirs, reducer)
    xl, yl = _axis_labels(reducer, red)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(coord_df["D1"], coord_df["D2"], s=70, marker="o",
               color="steelblue", alpha=0.85, edgecolors="white",
               linewidths=0.5, zorder=3)

    if label_points:
        for name in coord_df.index:
            x, y = coord_df.loc[name, ["D1", "D2"]]
            ax.annotate(name, (x, y), fontsize=6, ha="left", va="bottom",
                        xytext=(3, 3), textcoords="offset points", color="dimgray")

    ax.set_xlabel(xl, fontsize=11)
    ax.set_ylabel(yl, fontsize=11)
    ax.set_title(f"CAV direction space ({reducer.upper()})\n{cav_pattern}", fontsize=12)
    if reducer == "pca":
        ax.axhline(0, color="lightgray", lw=0.5)
        ax.axvline(0, color="lightgray", lw=0.5)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved static plot: {out_path}")


# ---------------------------------------------------------------------------
# Interactive plot (Plotly)
# ---------------------------------------------------------------------------

def plot_direction_map_interactive(
    cav_pattern: str,
    out_path: str,
    reducer: str = "pca",
    dims: int = 2,
    pfam_annotations: Optional[str] = None,
    clan_annotations: Optional[str] = None,
    distance_mode: str = "approximate",
    n_neighbors: Optional[int] = None,
    min_dist: Optional[float] = None,
):
    """Embed CAV direction vectors in 2D or 3D and save an interactive Plotly HTML."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("plotly is not installed. Run: pip install plotly")

    if dims not in (2, 3):
        raise ValueError("--dims must be 2 or 3")

    cav_dirs = load_directions_from_pattern(cav_pattern)
    coord_df, red = _embed(cav_dirs, reducer, dims=dims, distance_mode=distance_mode,
                           n_neighbors=n_neighbors, min_dist=min_dist)
    xl, yl = _axis_labels(reducer, red)
    zl = yl.replace("2", "3")  # e.g. "UMAP 3", "PC 3"

    # Load annotations
    annot = load_pfam_annotations(pfam_annotations) if pfam_annotations else None
    clans = load_clan_annotations(clan_annotations) if clan_annotations else None

    # Assign per-point clan colors
    NO_CLAN = "no clan"
    clan_name_per_point = []
    for acc in coord_df.index:
        if clans is not None and acc in clans.index:
            clan_name_per_point.append(clans.loc[acc, "clan_name"])
        else:
            clan_name_per_point.append(NO_CLAN)

    color_map = _clan_colors([c for c in clan_name_per_point if c != NO_CLAN])
    color_map[NO_CLAN] = "#aaaaaa"
    point_colors = [color_map[c] for c in clan_name_per_point]

    # Build hover text and search strings
    hover_texts = []
    search_strings = []
    for acc, clan_name in zip(coord_df.index, clan_name_per_point):
        clan_line = ""
        if clan_name != NO_CLAN and clans is not None and acc in clans.index:
            clan_acc = clans.loc[acc, "clan_acc"]
            clan_line = f"Clan: {clan_name} ({clan_acc})<br>"

        if annot is not None and acc in annot.index:
            row = annot.loc[acc]
            long = row["long_description"]
            wrapped = "<br>".join(
                long[i:i+80] for i in range(0, min(len(long), 400), 80)
            ) + ("…" if len(long) > 400 else "")
            text = (f"<b>{acc}</b> · {row['short_name']}<br>"
                    f"<i>{row['long_name']}</i><br>"
                    f"{clan_line}<br>"
                    f"{wrapped}")
            s = f"{acc} {row['short_name']} {row['long_name']} {row['long_description']} {clan_name}"
        else:
            text = f"<b>{acc}</b><br>{clan_line}"
            s = f"{acc} {clan_name}"
        hover_texts.append(text)
        search_strings.append(s.lower())

    marker = dict(size=8 if dims == 2 else 5, color=point_colors,
                  opacity=0.85 if dims == 2 else 0.6,
                  line=dict(width=0.5 if dims == 2 else 0, color="white"))

    if dims == 3:
        trace = go.Scatter3d(
            x=coord_df["D1"], y=coord_df["D2"], z=coord_df["D3"],
            mode="markers",
            marker=marker,
            text=hover_texts,
            hovertemplate="%{text}<extra></extra>",
        )
        layout = go.Layout(
            title=dict(text=f"CAV direction space ({reducer.upper()}, 3D)", font_size=15),
            scene=dict(
                xaxis_title=xl, yaxis_title=yl, zaxis_title=zl,
                xaxis=dict(showgrid=True, zeroline=False),
                yaxis=dict(showgrid=True, zeroline=False),
                zaxis=dict(showgrid=True, zeroline=False),
            ),
            hoverlabel=dict(bgcolor="white", font_size=12),
            paper_bgcolor="white",
            width=950, height=800,
        )
    else:
        trace = go.Scatter(
            x=coord_df["D1"], y=coord_df["D2"],
            mode="markers",
            marker=marker,
            text=hover_texts,
            hovertemplate="%{text}<extra></extra>",
        )
        layout = go.Layout(
            title=dict(text=f"CAV direction space ({reducer.upper()})", font_size=15),
            xaxis=dict(title=xl, showgrid=False, zeroline=(reducer == "pca")),
            yaxis=dict(title=yl, showgrid=False, zeroline=(reducer == "pca")),
            hoverlabel=dict(bgcolor="white", font_size=12),
            plot_bgcolor="white",
            paper_bgcolor="white",
            width=950, height=800,
        )

    fig = go.Figure(data=[trace], layout=layout)

    # Render the plot div (no full HTML yet — we'll wrap it ourselves)
    plot_div = fig.to_html(
        include_plotlyjs="cdn",
        full_html=False,
        div_id="cavplot",
    )

    import json as _json
    search_js     = _json.dumps(search_strings)
    base_colors_js = _json.dumps(point_colors)   # original hex per point

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
  body {{ font-family: sans-serif; margin: 20px; }}
  #search-bar {{
    display: flex; align-items: center; gap: 10px; margin-bottom: 12px;
  }}
  #search-box {{
    padding: 7px 12px; font-size: 15px; border: 1px solid #ccc;
    border-radius: 6px; width: 320px;
  }}
  #match-count {{ color: #666; font-size: 13px; }}
  #clear-btn {{
    padding: 6px 12px; font-size: 13px; border: 1px solid #bbb;
    border-radius: 6px; cursor: pointer; background: #f5f5f5;
  }}
</style>
</head>
<body>
<div id="search-bar">
  <input id="search-box" type="text" placeholder="Search PFAM (e.g. zinc finger, kinase, PF00001)…" oninput="filterPlot(this.value)">
  <button id="clear-btn" onclick="clearSearch()">Clear</button>
  <span id="match-count"></span>
</div>
{plot_div}
<script>
const searchStrings = {search_js};
const baseColors   = {base_colors_js};
const DIM_COLOR    = 'rgba(200,200,200,0.15)';

function filterPlot(query) {{
  query = query.trim().toLowerCase();
  const n = searchStrings.length;
  const colors = new Array(n);
  let matches = 0;
  if (query === '') {{
    for (let i = 0; i < n; i++) colors[i] = baseColors[i];
    matches = n;
  }} else {{
    for (let i = 0; i < n; i++) {{
      const hit = searchStrings[i].includes(query);
      colors[i] = hit ? baseColors[i] : DIM_COLOR;
      if (hit) matches++;
    }}
  }}
  Plotly.restyle('cavplot', {{'marker.color': [colors]}});
  document.getElementById('match-count').textContent =
    query ? matches + ' match' + (matches !== 1 ? 'es' : '') : '';
}}

function clearSearch() {{
  const box = document.getElementById('search-box');
  box.value = '';
  filterPlot('');
}}
</script>
</body>
</html>"""

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text(html, encoding="utf-8")
    logger.info(f"Saved interactive {dims}D plot: {out_path}")


# ---------------------------------------------------------------------------
# Rotating 3-D GIF
# ---------------------------------------------------------------------------

def plot_direction_map_gif(
    cav_pattern: str,
    out_path: str,
    reducer: str = "umap",
    clan_annotations: Optional[str] = None,
    distance_mode: str = "approximate",
    n_neighbors: Optional[int] = None,
    min_dist: Optional[float] = None,
    n_frames: int = 72,       # 72 frames = 5° per step → full 360°
    fps: int = 20,
    elev: float = 25,         # fixed elevation angle
    point_size: int = 20,
    figsize: Tuple[int, int] = (7, 7),
):
    """
    Render a rotating 3-D scatter of CAV direction vectors and save as an
    animated GIF.  Axes, panes, grid lines, and tick labels are all hidden
    so only the coloured point cloud is visible — clean for presentations.

    Requires: matplotlib, Pillow  (pip install Pillow)
    """
    try:
        from matplotlib.animation import FuncAnimation, PillowWriter
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 – registers projection
    except ImportError as e:
        raise ImportError(f"Missing dependency: {e}. Run: pip install Pillow") from e

    cav_dirs = load_directions_from_pattern(cav_pattern)
    coord_df, _ = _embed(cav_dirs, reducer, dims=3, distance_mode=distance_mode,
                         n_neighbors=n_neighbors, min_dist=min_dist)

    # Clan colours (same logic as interactive plot)
    clans = load_clan_annotations(clan_annotations) if clan_annotations else None
    NO_CLAN = "no clan"
    clan_name_per_point = [
        clans.loc[acc, "clan_name"] if (clans is not None and acc in clans.index) else NO_CLAN
        for acc in coord_df.index
    ]
    color_map = _clan_colors([c for c in clan_name_per_point if c != NO_CLAN])
    color_map[NO_CLAN] = "#aaaaaa"
    point_colors = [color_map[c] for c in clan_name_per_point]

    x, y, z = coord_df["D1"].values, coord_df["D2"].values, coord_df["D3"].values

    fig = plt.figure(figsize=figsize, facecolor="white")
    ax  = fig.add_subplot(111, projection="3d", facecolor="white")

    ax.scatter(x, y, z, c=point_colors, s=point_size, alpha=0.6, linewidths=0)

    # Remove all axes decoration
    ax.set_axis_off()
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor("none")

    plt.tight_layout(pad=0)

    def _update(frame):
        ax.view_init(elev=elev, azim=frame * (360.0 / n_frames))
        return (fig,)

    logger.info(f"Rendering {n_frames}-frame GIF ({fps} fps) …")
    anim = FuncAnimation(fig, _update, frames=n_frames, interval=1000 // fps, blit=False)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    anim.save(out_path, writer=PillowWriter(fps=fps))
    plt.close()
    logger.info(f"Saved rotating GIF: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Visualize CAV directions matched by a glob pattern.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--cav-pattern", required=True,
                        help="Glob pattern matching CAV .npy files, "
                             "e.g. '/path/to/cavs/*/L25_concept_v1.npy'.")
    parser.add_argument("--out", required=True,
                        help="Output file path (.png for static, .html for interactive).")
    parser.add_argument("--reducer", default="pca",
                        choices=["umap", "tsne", "pca"],
                        help="Dimensionality reduction algorithm (default: pca).")
    parser.add_argument("--interactive", action="store_true",
                        help="Write an interactive Plotly HTML instead of a static PNG.")
    parser.add_argument("--dims", type=int, default=2, choices=[2, 3],
                        help="Embedding dimensions for interactive plot (default: 2).")
    parser.add_argument("--pfam-annotations",
                        help="Path to pfamA.txt for hover annotations "
                             "(accession, short name, description).")
    parser.add_argument("--clan-annotations",
                        help="Path to clan reference file for coloring points by clan.")
    parser.add_argument("--no-labels", action="store_true",
                        help="Suppress concept name labels on static plot.")
    parser.add_argument("--gif-out",
                        help="Also save a rotating 3-D GIF to this path (requires --reducer).")
    parser.add_argument("--n-neighbors", type=int, default=None,
                        help="UMAP n_neighbors (default: 30 for approximate/precomputed, "
                             "10 for cophenetic).")
    parser.add_argument("--min-dist", type=float, default=None,
                        help="UMAP min_dist (default: 0.1 for approximate/precomputed, "
                             "0.0 for cophenetic).")
    parser.add_argument("--distance-mode", default="approximate",
                        choices=["approximate", "precomputed", "cophenetic"],
                        help="How to compute distances for embedding. "
                             "'approximate': UMAP/t-SNE on raw vectors (fast, default). "
                             "'precomputed': exact all-by-all cosine distance matrix. "
                             "'cophenetic': distances from average-linkage dendrogram — "
                             "best resolves hierarchical cluster structure (not with pca).")
    parser.add_argument("--gif-frames", type=int, default=72,
                        help="Number of frames in the GIF (default: 72 = 5° per step).")
    parser.add_argument("--gif-fps", type=int, default=20,
                        help="Frames per second of the GIF (default: 20).")
    args = parser.parse_args()

    if args.gif_out:
        plot_direction_map_gif(
            cav_pattern=args.cav_pattern,
            out_path=args.gif_out,
            reducer=args.reducer,
            clan_annotations=args.clan_annotations,
            distance_mode=args.distance_mode,
            n_neighbors=args.n_neighbors,
            min_dist=args.min_dist,
            n_frames=args.gif_frames,
            fps=args.gif_fps,
        )

    if args.interactive:
        plot_direction_map_interactive(
            cav_pattern=args.cav_pattern,
            out_path=args.out,
            reducer=args.reducer,
            dims=args.dims,
            pfam_annotations=args.pfam_annotations,
            clan_annotations=args.clan_annotations,
            distance_mode=args.distance_mode,
            n_neighbors=args.n_neighbors,
            min_dist=args.min_dist,
        )
    else:
        plot_direction_map(
            cav_pattern=args.cav_pattern,
            out_path=args.out,
            reducer=args.reducer,
            label_points=not args.no_labels,
        )


if __name__ == "__main__":
    main()
