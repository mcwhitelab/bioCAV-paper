#!/usr/bin/env python3
"""
cav_continuum_viz.py — Visualise the normal→cancer transcriptional continuum
for a set of CAV pairs.

For each pair the figure shows:
  TOP STRIP  — every cell placed on the L2 score axis, coloured blue (normal)
               or red (cancer).  Lets you see how well the two populations
               separate and whether intermediate cells exist.
  CURVES     — smoothed expression of the top-N genes along the same axis.
               Positive-r genes slope upward left→right (more expressed in
               cancer); negative-r genes slope downward.

Pairs are specified via --pairs (e.g. the four skin/melanoma cell types).
One panel column is drawn per pair so they share the same x-axis framing and
can be compared directly.

Usage
-----
python specific_scripts/cav_continuum_viz.py \\
    --coords      results/hierarchy/cell_coordinates.tsv \\
    --h5ad        data/cells.h5ad \\
    --gene-corr-dir  results/gene_correlation/ \\
    --group-col   cell_type \\
    --context-col tissue \\
    --condition-col disease \\
    --level       L2 \\
    --pairs  endothelial_cell__skin_epidermis__normal_vs_melanoma \\
             epithelial_cell__skin_epidermis__normal_vs_melanoma \\
             fibroblast__skin_epidermis__normal_vs_melanoma \\
             T_cell__skin_epidermis__normal_vs_melanoma \\
    --out  results/figures/cav_continuum_melanoma.png
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Condition colours
NORMAL_COLOR    = "#3a6fad"   # blue
CANCER_COLOR    = "#c0392b"   # red
STRIP_ALPHA     = 0.6
CURVE_POS_CMAP  = "YlOrRd"   # warm  — positive-r genes
CURVE_NEG_CMAP  = "YlGnBu"   # cool  — negative-r genes


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def norm(v: str) -> str:
    return v.lower().replace("_", " ").strip()


def safe_labels(series) -> np.ndarray:
    return series.astype(str).replace("nan", "unknown").values


def cells_matching(obs: pd.DataFrame,
                   group_col: str, context_col: str, condition_col: str,
                   group_val: str, context_val: str, cond_val: str) -> np.ndarray:
    def match(col, val):
        if not col or col not in obs.columns:
            return np.ones(len(obs), bool)
        labels = safe_labels(obs[col])
        nval = norm(val)
        return np.array([norm(l) == nval for l in labels])
    return (match(group_col, group_val) &
            match(context_col, context_val) &
            match(condition_col, cond_val))


def gaussian_smooth(x: np.ndarray, y: np.ndarray,
                    x_grid: np.ndarray, bw: float) -> np.ndarray:
    """Nadaraya-Watson Gaussian kernel smoother."""
    out = np.zeros(len(x_grid))
    for i, xi in enumerate(x_grid):
        w = np.exp(-0.5 * ((x - xi) / bw) ** 2)
        if w.sum() > 1e-12:
            out[i] = (w * y).sum() / w.sum()
    return out


def minmax(arr: np.ndarray, clip_pct: float = 99.0) -> np.ndarray:
    """Min-max scale with percentile clipping to avoid outlier squash."""
    lo  = np.percentile(arr, 100 - clip_pct)
    hi  = np.percentile(arr, clip_pct)
    arr = np.clip(arr, lo, hi)
    return (arr - lo) / (hi - lo + 1e-9)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_coords(coords_path: str, level: str) -> pd.DataFrame:
    logger.info(f"Loading coordinates: {coords_path}")
    df = pd.read_csv(coords_path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    level_cols = [c for c in df.columns if c.startswith(level + "__")]
    logger.info(f"  {len(df)} cells, {len(level_cols)} {level} axes")
    return df


def load_expression(h5ad_path: str) -> Tuple[np.ndarray, List[str],
                                               List[str], Dict[str, str]]:
    import anndata as ad
    logger.info(f"Loading h5ad: {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)

    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)

    gene_ids   = list(adata.var_names)
    cell_ids   = list(adata.obs_names)
    logger.info(f"  {X.shape[0]} cells × {X.shape[1]} genes")

    gene_name_map: Dict[str, str] = {}
    for candidate in ["feature_name", "gene_name", "gene_symbol",
                      "symbol", "name", "Gene", "gene"]:
        if candidate in adata.var.columns:
            gene_name_map = dict(zip(gene_ids,
                                     adata.var[candidate].astype(str).values))
            logger.info(f"  Gene names from column '{candidate}'")
            break

    return X, gene_ids, cell_ids, gene_name_map


def load_obs(h5ad_path: str, cols: List[str]) -> pd.DataFrame:
    import anndata as ad
    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs[cols].copy()
    adata.file.close()
    obs.index = obs.index.astype(str)
    return obs


def load_gene_corr(tsv_path: Path, top_n: int,
                   gene_name_map: Dict[str, str],
                   top_n_each: int = None) -> pd.DataFrame:
    """Load a gene correlation TSV, return selected rows with display names.

    If top_n_each is set, returns top-N by positive r + top-N by negative r.
    Otherwise returns top-N by |r| (already sorted that way by cav_gene_correlation).
    """
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as e:
        logger.warning(f"  Could not read {tsv_path.name}: {e}")
        return pd.DataFrame()
    if df.empty:
        logger.warning(f"  Gene corr file is empty: {tsv_path.name}")
        return pd.DataFrame()

    if top_n_each is not None:
        pos = df[df["r"] >= 0].nlargest(top_n_each, "r")
        neg = df[df["r"] <  0].nsmallest(top_n_each, "r")
        df  = pd.concat([pos, neg]).reset_index(drop=True)
    else:
        df = df.head(top_n).copy()
    if "gene_name" not in df.columns:
        df["gene_name"] = df["gene"].map(gene_name_map).fillna(df["gene"])
    df["display"] = df["gene_name"].where(
        df["gene_name"].notna() & (df["gene_name"] != ""),
        df["gene"]
    )
    return df


# ---------------------------------------------------------------------------
# Per-pair data assembly
# ---------------------------------------------------------------------------

def get_pair_data(pair_label: str,
                  coords: pd.DataFrame,
                  X_expr: np.ndarray,
                  gene_ids: List[str],
                  gene_name_map: Dict[str, str],
                  obs: pd.DataFrame,
                  group_col: str, context_col: str, condition_col: str,
                  level: str,
                  gene_corr_dir: Path,
                  top_n: int,
                  max_cells: int,
                  baseline_value: str,
                  top_n_each: int = None) -> Optional[dict]:
    """
    Returns a dict with everything needed to draw one panel, or None on failure.
    """
    # Parse pair_label: group__context__baseline_vs_condition
    try:
        left, right = pair_label.rsplit("_vs_", 1)
        parts       = left.split("__")
        cond_val    = right
        b_val       = parts[-1]
        ctx_val     = parts[-2] if len(parts) >= 3 else ""
        g_val       = parts[0]
    except Exception:
        logger.error(f"Cannot parse pair label '{pair_label}' — "
                     f"expected format group__context__baseline_vs_condition")
        return None

    # Cell masks
    base_mask = cells_matching(obs, group_col, context_col, condition_col,
                               g_val, ctx_val, b_val)
    cond_mask = cells_matching(obs, group_col, context_col, condition_col,
                               g_val, ctx_val, cond_val)
    pair_mask = base_mask | cond_mask
    n_pair    = pair_mask.sum()

    if n_pair < 10:
        logger.warning(f"  {pair_label}: only {n_pair} cells — skipping")
        return None

    logger.info(f"  {pair_label}: {base_mask.sum()} normal + "
                f"{cond_mask.sum()} cancer = {n_pair} cells")

    # Align coords & obs on shared string IDs
    shared_ids = obs.index.intersection(coords.index)
    obs_a      = obs.reindex(shared_ids)
    coords_a   = coords.reindex(shared_ids)

    # Re-derive masks on aligned subset
    base_mask_a = cells_matching(obs_a, group_col, context_col, condition_col,
                                 g_val, ctx_val, b_val)
    cond_mask_a = cells_matching(obs_a, group_col, context_col, condition_col,
                                 g_val, ctx_val, cond_val)
    pair_mask_a = base_mask_a | cond_mask_a

    # Score column
    score_col = f"{level}__{g_val}__{ctx_val}__{cond_val}"
    if score_col not in coords_a.columns:
        score_col = f"{level}__{g_val}__{ctx_val}__{b_val}"
    if score_col not in coords_a.columns:
        available = [c for c in coords_a.columns if level in c][:4]
        logger.warning(f"  Score column not found. Sample cols: {available}")
        return None

    scores_all = coords_a[score_col].values[pair_mask_a].astype(np.float32)
    labels_all = np.where(
        cond_mask_a[pair_mask_a], cond_val, b_val
    )

    # Expression: align h5ad cell_ids with obs/coords
    expr_cell_ids = list(obs.index)   # obs was loaded from h5ad, same order
    id_to_row     = {cid: i for i, cid in enumerate(expr_cell_ids)}

    pair_ids  = shared_ids[pair_mask_a]
    expr_rows = np.array([id_to_row[cid] for cid in pair_ids
                          if cid in id_to_row])
    valid_ids = [cid for cid in pair_ids if cid in id_to_row]

    if len(expr_rows) < 10:
        logger.warning(f"  {pair_label}: too few expression rows")
        return None

    # Re-align scores/labels to cells that also have expression
    keep_mask = np.array([cid in id_to_row for cid in pair_ids])
    scores_all = scores_all[keep_mask]
    labels_all = labels_all[keep_mask]

    # Subsample if needed
    if len(scores_all) > max_cells:
        rng  = np.random.default_rng(42)
        idx  = rng.choice(len(scores_all), max_cells, replace=False)
        scores_all = scores_all[idx]
        labels_all = labels_all[idx]
        expr_rows  = expr_rows[idx]

    X_pair = X_expr[expr_rows]   # (n_cells, n_genes)

    # Gene correlation file
    tsv_name = f"{pair_label}.tsv"
    tsv_path = gene_corr_dir / tsv_name
    if not tsv_path.exists():
        logger.warning(f"  Gene corr file not found: {tsv_path}")
        return None

    gene_df = load_gene_corr(tsv_path, top_n, gene_name_map,
                              top_n_each=top_n_each)
    if gene_df.empty:
        return None

    # Pull expression for top genes
    gene_id_to_col = {g: i for i, g in enumerate(gene_ids)}
    gene_expr = {}
    for _, row in gene_df.iterrows():
        gid = row["gene"]
        if gid in gene_id_to_col:
            gene_expr[gid] = X_pair[:, gene_id_to_col[gid]]

    return {
        "pair_label": pair_label,
        "group":      g_val,
        "context":    ctx_val,
        "baseline":   b_val,
        "condition":  cond_val,
        "scores":     scores_all,
        "labels":     labels_all,
        "gene_df":    gene_df,
        "gene_expr":  gene_expr,
        "n_normal":   (labels_all == b_val).sum(),
        "n_cancer":   (labels_all == cond_val).sum(),
    }


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def draw_strip(ax, scores: np.ndarray, labels: np.ndarray,
               baseline: str, condition: str, x_lim: tuple = None):
    """Top strip: cells as dots on L2 axis, coloured by condition."""
    rng    = np.random.default_rng(0)
    jitter = rng.uniform(-0.4, 0.4, len(scores))

    is_base = labels == baseline
    ax.scatter(scores[is_base],  jitter[is_base],
               c=NORMAL_COLOR, s=5, alpha=STRIP_ALPHA, linewidths=0,
               label=baseline)
    ax.scatter(scores[~is_base], jitter[~is_base],
               c=CANCER_COLOR, s=5, alpha=STRIP_ALPHA, linewidths=0,
               label=condition)

    ax.set_xlim(x_lim if x_lim else (scores.min() - 0.05, scores.max() + 0.05))
    ax.set_ylim(-1, 1)
    ax.axis("off")


def draw_curves(ax, scores: np.ndarray, gene_df: pd.DataFrame,
                gene_expr: dict, bw_frac: float = 0.08,
                clip_pct: float = 99.0,
                x_lim: tuple = None):
    """Smoothed expression curves for each gene along the L2 axis."""
    x_min, x_max = scores.min(), scores.max()
    bw     = (x_max - x_min) * bw_frac
    x_grid = np.linspace(x_min, x_max, 120)

    pos_genes = gene_df[gene_df["r"] >= 0]
    neg_genes = gene_df[gene_df["r"] <  0]

    pos_cmap = plt.colormaps.get_cmap(CURVE_POS_CMAP)
    neg_cmap = plt.colormaps.get_cmap(CURVE_NEG_CMAP)

    legend_handles = []

    def _draw_group(gdf, cmap, start=0.35, end=0.95):
        n = len(gdf)
        for i, (_, row) in enumerate(gdf.iterrows()):
            gid = row["gene"]
            if gid not in gene_expr:
                continue
            expr_raw = gene_expr[gid].astype(np.float32)
            expr_sc  = minmax(expr_raw, clip_pct=clip_pct)

            color = cmap(start + (end - start) * (i / max(n - 1, 1)))
            y_sm  = gaussian_smooth(scores, expr_sc, x_grid, bw)
            line, = ax.plot(x_grid, y_sm, color=color, lw=1.6, alpha=0.85)

            label = f"{row['display']}  (r={row['r']:.2f})"
            legend_handles.append(Line2D([0], [0], color=color, lw=2,
                                         label=label))

    _draw_group(pos_genes, pos_cmap)
    _draw_group(neg_genes, neg_cmap)

    # Vertical line at score=0 if it's in range
    if x_min < 0 < x_max:
        ax.axvline(0, color="#aaaaaa", lw=0.8, ls="--", zorder=0)

    plot_xmin = x_lim[0] if x_lim else x_min - 0.05
    plot_xmax = x_lim[1] if x_lim else x_max + 0.05
    ax.set_xlim(plot_xmin, plot_xmax)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("L2 score  (← normal    cancer →)", fontsize=8)
    ax.set_ylabel("Scaled expression", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)

    return legend_handles


# ---------------------------------------------------------------------------
# Gene strips style
# ---------------------------------------------------------------------------

def draw_gene_strips(ax, scores: np.ndarray, gene_df: pd.DataFrame,
                     gene_expr: dict, clip_pct: float = 99.0,
                     x_lim: tuple = None):
    """
    Stacked expression strips — one per gene, cells coloured by expression.
    No smoothing, no extrapolation: every dot is a real cell at its real L2 score.
    Positive-r genes use a warm colormap; negative-r use a cool colormap.
    Gene labels sit on the y-axis.
    """
    rng      = np.random.default_rng(1)
    n_genes  = len(gene_df)
    strip_h  = 0.8   # half-height of each strip in data coords

    pos_cmap = plt.colormaps.get_cmap("YlOrRd")
    neg_cmap = plt.colormaps.get_cmap("YlGnBu")

    yticks, ylabels = [], []

    for i, (_, row) in enumerate(gene_df.iterrows()):
        gid   = row["gene"]
        y_ctr = n_genes - i          # top gene at top
        if gid not in gene_expr:
            continue

        expr_raw = gene_expr[gid].astype(np.float32)
        expr_sc  = minmax(expr_raw, clip_pct=clip_pct)

        cmap   = pos_cmap if row["r"] >= 0 else neg_cmap
        colors = cmap(0.2 + 0.75 * expr_sc)   # avoid very pale end

        jitter = rng.uniform(-strip_h * 0.45, strip_h * 0.45, len(scores))
        ax.scatter(scores, y_ctr + jitter, c=colors,
                   s=3, alpha=0.6, linewidths=0, rasterized=True)

        yticks.append(y_ctr)
        ylabels.append(f"{row['display']}  (r={row['r']:.2f})")

    # Dividing lines between strips
    for i in range(n_genes + 1):
        y = n_genes - i + 0.5
        ax.axhline(y, color="#e0e0e0", lw=0.5, zorder=0)

    # Vertical line at score=0
    x_min = scores.min()
    x_max = scores.max()
    if x_min < 0 < x_max:
        ax.axvline(0, color="#aaaaaa", lw=0.8, ls="--", zorder=1)

    plot_xmin = x_lim[0] if x_lim else x_min - 0.05
    plot_xmax = x_lim[1] if x_lim else x_max + 0.05
    ax.set_xlim(plot_xmin, plot_xmax)
    ax.set_ylim(0.5, n_genes + 0.5)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("L2 score  (← normal    cancer →)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right", "left"]].set_visible(False)

    # Small colour bar legend for expression
    sm_pos = plt.cm.ScalarMappable(cmap=pos_cmap,
                                    norm=plt.Normalize(vmin=0, vmax=1))
    sm_neg = plt.cm.ScalarMappable(cmap=neg_cmap,
                                    norm=plt.Normalize(vmin=0, vmax=1))
    return sm_pos, sm_neg


# ---------------------------------------------------------------------------
# Main figure
# ---------------------------------------------------------------------------

def make_figure(pair_data_list: List[dict],
                out_path: str,
                panel_w: float = 4.5,
                strip_h: float = 0.6,
                curves_h: float = 4.0,
                title: str = "",
                clip_pct: float = 99.0,
                shared_xaxis: bool = False,
                style: str = "strips"):

    n_pairs = len(pair_data_list)
    fig_w   = panel_w * n_pairs
    fig_h   = strip_h + curves_h + 1.2   # extra for legends

    # Shared x-axis range across all panels
    if shared_xaxis:
        all_scores = np.concatenate([pd_["scores"] for pd_ in pair_data_list])
        x_lim = (all_scores.min() - 0.2, all_scores.max() + 0.2)
    else:
        x_lim = None

    fig = plt.figure(figsize=(fig_w, fig_h))
    top_margin = 0.88 if not title else 0.85
    if title:
        fig.suptitle(title, fontsize=12, y=0.98)

    outer_gs = gridspec.GridSpec(
        1, n_pairs, figure=fig,
        left=0.06, right=0.97, top=top_margin, bottom=0.22,
        wspace=0.35,
    )

    for col_idx, pd_ in enumerate(pair_data_list):
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1,
            subplot_spec=outer_gs[col_idx],
            height_ratios=[strip_h, curves_h],
            hspace=0.05,
        )

        ax_strip  = fig.add_subplot(inner[0])
        ax_curves = fig.add_subplot(inner[1])

        # Title for panel
        title_str = (f"{pd_['group'].replace('_', ' ')}\n"
                     f"{pd_['context'].replace('_', ' ')}\n"
                     f"n={pd_['n_normal']} normal, "
                     f"n={pd_['n_cancer']} {pd_['condition'].replace('_', ' ')}")
        ax_strip.set_title(title_str, fontsize=8, pad=4)

        draw_strip(ax_strip, pd_["scores"], pd_["labels"],
                   pd_["baseline"], pd_["condition"], x_lim=x_lim)

        if style == "strips":
            draw_gene_strips(ax_curves, pd_["scores"],
                             pd_["gene_df"], pd_["gene_expr"],
                             clip_pct=clip_pct, x_lim=x_lim)
        else:
            legend_handles = draw_curves(ax_curves, pd_["scores"],
                                         pd_["gene_df"], pd_["gene_expr"],
                                         clip_pct=clip_pct, x_lim=x_lim)
            ax_curves.legend(
                handles=legend_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.28),
                fontsize=6.5,
                frameon=False,
                ncol=1,
            )

    # Shared strip legend (normal / cancer dots) — bottom right
    shared_legend = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=NORMAL_COLOR,
               markersize=6, label="normal"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=CANCER_COLOR,
               markersize=6, label="cancer"),
    ]
    fig.legend(handles=shared_legend, loc="lower right",
               bbox_to_anchor=(0.99, 0.01), fontsize=8, frameon=False,
               title="Cell strip", title_fontsize=8)

    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved figure: {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="CAV continuum visualisation — cells on L2 axis + "
                    "smoothed gene expression curves.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--coords",        required=True,
                        help="cell_coordinates.tsv from cav_hierarchy.py")
    parser.add_argument("--h5ad",          required=True,
                        help="AnnData h5ad with expression in .X")
    parser.add_argument("--gene-corr-dir", required=True,
                        help="Directory of per-pair TSVs from cav_gene_correlation.py")
    parser.add_argument("--pairs",         required=True, nargs="+",
                        help="Pair labels to plot, e.g. "
                             "T_cell__skin_epidermis__normal_vs_melanoma")
    parser.add_argument("--group-col",     default="cell_type")
    parser.add_argument("--context-col",   default="tissue")
    parser.add_argument("--condition-col", default="disease")
    parser.add_argument("--baseline-value", default="normal")
    parser.add_argument("--level",         default="L2",
                        choices=["L0", "L1", "L2", "delta"])
    parser.add_argument("--top-n",         type=int, default=8,
                        help="Top genes by |r| to plot (default: 8). "
                             "Overridden by --top-n-each.")
    parser.add_argument("--top-n-each",    type=int, default=None,
                        help="Plot top-N positive AND top-N negative genes "
                             "separately (e.g. --top-n-each 4 gives 4 up + 4 down). "
                             "Overrides --top-n.")
    parser.add_argument("--max-cells",     type=int, default=1000,
                        help="Max cells per pair to plot (default: 1000)")
    parser.add_argument("--bw-frac",       type=float, default=0.08,
                        help="Kernel bandwidth as fraction of score range "
                             "(default: 0.08 — increase to smooth more)")
    parser.add_argument("--out",           required=True,
                        help="Output figure path (.png)")
    parser.add_argument("--title",         default="",
                        help="Figure title")
    parser.add_argument("--panel-width",   type=float, default=4.5)
    parser.add_argument("--curves-height", type=float, default=4.0)
    parser.add_argument("--style",
                        choices=["strips", "curves"], default="strips",
                        help="'strips' (default): one expression dot-strip per gene, "
                             "no extrapolation. 'curves': smoothed expression curves.")
    parser.add_argument("--clip-pct",      type=float, default=99.0,
                        help="Percentile for expression clipping before scaling "
                             "(default: 99 — raise to show more outlier signal, "
                             "lower to compress more)")
    parser.add_argument("--shared-xaxis",  action="store_true",
                        help="Use the same x-axis range across all panels.")
    args = parser.parse_args()

    # Load shared data
    coords = load_coords(args.coords, args.level)

    obs_cols = list({args.group_col, args.context_col, args.condition_col})
    obs      = load_obs(args.h5ad, obs_cols)

    X_expr, gene_ids, expr_cell_ids, gene_name_map = load_expression(args.h5ad)
    # obs and X_expr share the same row order (both from h5ad)

    gene_corr_dir = Path(args.gene_corr_dir)

    # Build per-pair data
    pair_data_list = []
    for pair_label in args.pairs:
        logger.info(f"Processing pair: {pair_label}")
        pd_ = get_pair_data(
            pair_label        = pair_label,
            coords            = coords,
            X_expr            = X_expr,
            gene_ids          = gene_ids,
            gene_name_map     = gene_name_map,
            obs               = obs,
            group_col         = args.group_col,
            context_col       = args.context_col,
            condition_col     = args.condition_col,
            level             = args.level,
            gene_corr_dir     = gene_corr_dir,
            top_n             = args.top_n,
            top_n_each        = args.top_n_each,
            max_cells         = args.max_cells,
            baseline_value    = args.baseline_value,
        )
        if pd_ is not None:
            pair_data_list.append(pd_)

    if not pair_data_list:
        logger.error("No pairs loaded — check pair labels and file paths.")
        sys.exit(1)

    logger.info(f"Drawing {len(pair_data_list)} panels → {args.out}")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    make_figure(
        pair_data_list,
        out_path      = args.out,
        panel_w       = args.panel_width,
        curves_h      = args.curves_height,
        title         = args.title,
        clip_pct      = args.clip_pct,
        shared_xaxis  = args.shared_xaxis,
        style         = args.style,
    )


if __name__ == "__main__":
    main()
