#!/usr/bin/env python3
"""
Render Figure 4 Panel A: a single-angle (azimuth 330 deg) 3-D view of the
CAV direction-space UMAP, colored by Pfam clan, no title/legend (legend
removed per user -- manual labels added afterward).

Coordinates are decoded directly from the embedded Plotly data in the
precomputed motif_clustering/results/figures/pfam_umap_3d.html (the
default/approximate-distance-mode variant -- the "_sim" cophenetic
variant packs everything into a near-uniform disk with no visible
structure) rather than recomputed. Run from the figures/ directory:

    python figure_data/pfam_clan_tsne/render_pfam_clan_umap_3d.py
"""
import base64
import json
import re
import sys
from pathlib import Path

from PIL import Image, ImageChops

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Resolved relative to this file (<repo>/figures/figure_data/pfam_clan_tsne/)
# so the script runs from any checkout rather than one absolute location.
REPO_ROOT = Path(__file__).resolve().parents[2]
UMAP_3D_HTML = REPO_ROOT / "motif_clustering/results/figures/pfam_umap_3d.html"
CLAN_TSV = REPO_ROOT / "motif_clustering/results/figures/Pfam-A.clans.tsv"
OUT_PNG = REPO_ROOT / "figures/figures/pfam_clan_umap_3d_mosaic.png"

AZIMUTH = 330
ELEVATION = 20
GREY_ALPHA = 0.40
COLOR_ALPHA = 0.85


def _decode_plotly_array(field):
    return np.frombuffer(base64.b64decode(field["bdata"]), dtype=field["dtype"])


def load_umap_3d_coords(html_path: Path):
    """Extract x/y/z/color/text from a Plotly HTML's single scatter3d trace."""
    html = html_path.read_text()
    start = html.find("Plotly.newPlot")
    arr_start = html.find("[", start)
    depth, i, in_str, esc = 0, arr_start, False, False
    while i < len(html):
        c = html[i]
        if in_str:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == '"':
                in_str = False
        else:
            if c == '"':
                in_str = True
            elif c == "[":
                depth += 1
            elif c == "]":
                depth -= 1
                if depth == 0:
                    break
        i += 1
    trace = json.loads(html[arr_start : i + 1])[0]
    x = _decode_plotly_array(trace["x"])
    y = _decode_plotly_array(trace["y"])
    z = _decode_plotly_array(trace["z"])
    colors = trace["marker"]["color"]
    text = trace["text"]
    return x, y, z, colors, text


def main():
    clans = pd.read_csv(
        CLAN_TSV, sep="\t", header=None,
        names=["accession", "clan_acc", "clan_name", "clan_desc", "pfam_name"],
    ).set_index("accession")

    x, y, z, colors, text = load_umap_3d_coords(UMAP_3D_HTML)
    accessions = [re.search(r"<b>(PF\d+)</b>", t).group(1) for t in text]
    clan_names = [
        clans.loc[a, "clan_name"] if a in clans.index and pd.notna(clans.loc[a, "clan_name"]) else "no clan"
        for a in accessions
    ]
    df = pd.DataFrame({"x": x, "y": y, "z": z, "color": colors, "clan": clan_names})
    df["is_grey"] = df["clan"] == "no clan"

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    grey, color = df[df.is_grey], df[~df.is_grey]
    ax.scatter(grey["x"], grey["y"], grey["z"], s=10, c="#dddddd", alpha=GREY_ALPHA,
              edgecolors="none", depthshade=False)
    ax.scatter(color["x"], color["y"], color["z"], s=22, c=color["color"], alpha=COLOR_ALPHA,
              edgecolors="none", depthshade=False)
    ax.view_init(elev=ELEVATION, azim=AZIMUTH)
    ax.set_xlabel("UMAP 1", fontsize=9, labelpad=-8)
    ax.set_ylabel("UMAP 2", fontsize=9, labelpad=-8)
    ax.set_zlabel("UMAP 3", fontsize=9, labelpad=-8)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor((0.8, 0.8, 0.8, 0.5))

    plt.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    plt.close()

    # matplotlib's 3D axes leave asymmetric top/bottom margin that
    # bbox_inches="tight" doesn't fully remove (content used ~91% of
    # width but only ~74% of height) -- auto-crop remaining whitespace so
    # the panel fills the space it's given when embedded in the figure.
    _trim_whitespace(OUT_PNG)
    print(f"Saved {OUT_PNG}", file=sys.stderr)


def _trim_whitespace(png_path: Path, pad: int = 15):
    im = Image.open(png_path).convert("RGB")
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    bbox = ImageChops.difference(im, bg).getbbox()
    if bbox:
        l, t, r, b = bbox
        l, t = max(0, l - pad), max(0, t - pad)
        r, b = min(im.width, r + pad), min(im.height, b + pad)
        im.crop((l, t, r, b)).save(png_path)


if __name__ == "__main__":
    main()
