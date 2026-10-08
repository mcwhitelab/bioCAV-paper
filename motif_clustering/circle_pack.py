#!/usr/bin/env python3
"""
circle_pack.py — Circle-packed view of a CAV motif dendrogram, cut at several
heights.  Python port of pfam_circlepack.R + the interactive HTML builder, so
the whole thing runs on the cluster against any clustering cav_cluster.py
produces.

The cut_* columns all come from one dendrogram and are strictly nested, so they
form a hierarchy: origin -> cut levels (coarse to fine) -> individual motifs.
That hierarchy is circle-packed once and cached; the cache is then rendered
into a self-contained interactive HTML.

Mirrors the R in the details that matter:
  * cluster node ids are  f"{clusternum}{cutcolumn}"  (e.g. "7cut_k00314"),
    unique across cut levels
  * --hide-top-n suppresses the OUTLINE of the n largest top-cut clusters and
    of their children at the next cut down.  The packing is unchanged and every
    motif inside is still drawn — it only stops the giant rings swallowing the
    figure.  This is a draw-time filter on the payload's cluster array.
  * clan colours come from a fixed clan->colour table so a clan keeps its
    colour across every figure and cut set.  Pass the same --clan-colours CSV
    the R figures used to keep the two in sync; otherwise one is generated and
    written out for reuse.

The HTML is built by substituting into a template that must contain the
placeholder `/*__DATA__*/null`, plus `__TITLE__` and `__SUBTITLE__`.

Usage
-----
python circle_pack.py \\
    --clusters results/figures/circlepack_input_diffusion.csv \\
    --cuts cut_k00314 cut_k02769 cut_k04629 cut_k09172 cut_k18984 \\
    --template interactive_template.html \\
    --hide-top-n 3 \\
    --out results/figures/circlepack_diffusion.html
"""

import re
import sys
import json
import argparse
import colorsys
import logging
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "vendor"))
import circlify  # noqa: E402  (vendored, see vendor/README.md)

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

LABEL_FIELDS = ["accession", "pfam_short_name", "pfam_long_name",
                "clan_acc", "clan_name"]


def cut_display_name(col):
    """'cut_095' -> 'cut 0.95';  'cut_k00314' -> 'k = 314'."""
    m = re.fullmatch(r"cut_k0*(\d+)", col)
    if m:
        return f"k = {int(m.group(1)):,}"
    m = re.fullmatch(r"cut_(\d+)", col)
    if m:
        return f"cut {int(m.group(1)) / 100:.2f}"
    return col


# ---------------------------------------------------------------------------
# Hierarchy + packing
# ---------------------------------------------------------------------------

def check_nesting(df, cuts):
    for coarse, fine in zip(cuts, cuts[1:]):
        bad = int((df.groupby(fine)[coarse].nunique() > 1).sum())
        if bad:
            raise SystemExit(
                f"{bad} clusters in {fine} span more than one {coarse} cluster. "
                f"The cut columns are not nested, so they do not describe a "
                f"hierarchy. Are they all from the same dendrogram?")
    logger.info(f"Cut columns nest cleanly: {' -> '.join(cuts)}")


def build_hierarchy(df, cuts):
    """Nested list for circlify. Internal datum = number of motifs beneath."""
    def recurse(sub, depth):
        if depth == len(cuts):
            return [{"id": acc, "datum": 1} for acc in sub["accession"]]
        out = []
        for cid, grp in sub.groupby(cuts[depth], sort=True):
            out.append({"id": f"{cid}{cuts[depth]}",
                        "datum": int(len(grp)),
                        "children": recurse(grp, depth + 1)})
        return out
    return recurse(df, 0)


def pack(df, cuts, cache):
    """Circle-pack the hierarchy; cache because it is the slow step."""
    if cache and Path(cache).exists():
        lay = pd.read_csv(cache)
        logger.info(f"Loaded cached layout from {cache} "
                    f"({len(lay)} nodes, max depth {lay['depth'].max()})")
        return lay

    data = build_hierarchy(df, cuts)
    logger.info(f"Packing {len(df):,} motifs across {len(cuts)} cut levels "
                f"(this is the slow step) …")
    circles = circlify.circlify(data, show_enclosure=False)
    lay = pd.DataFrame({
        "name":  [c.ex["id"] for c in circles],
        "x":     [c.x for c in circles],
        "y":     [c.y for c in circles],
        "r":     [c.r for c in circles],
        "depth": [c.level for c in circles],
    })
    logger.info(f"Packed {len(lay)} nodes, max depth {lay['depth'].max()}")
    if cache:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        lay.to_csv(cache, index=False)
        logger.info(f"Cached layout → {cache}")
    return lay


# ---------------------------------------------------------------------------
# Colours
# ---------------------------------------------------------------------------

def clan_palette(clan_ids, path, seed=42):
    """
    Stable clan -> colour table.  Reuses `path` if it exists so a clan keeps its
    colour across figures; otherwise generates one and writes it there.

    812 clans is far past any palette that stays distinguishable, so hues are
    sampled rather than assigned: the point is seeing which branches share a
    clan, not looking any one clan up.
    """
    if path and Path(path).exists():
        tab = pd.read_csv(path)
        missing = set(clan_ids) - set(tab["clan_acc"])
        if missing:
            logger.warning(f"{len(missing)} clans absent from {path}; "
                           f"they will render grey")
        logger.info(f"Loaded {len(tab)} clan colours from {path}")
        return dict(zip(tab["clan_acc"], tab["colour"]))

    rng = np.random.default_rng(seed)
    hues = rng.permutation(np.linspace(0, 1, len(clan_ids), endpoint=False))
    sat  = rng.uniform(0.45, 0.75, len(clan_ids))
    lit  = rng.uniform(0.35, 0.65, len(clan_ids))
    cols = ["#%02x%02x%02x" % tuple(int(round(255 * v))
                                    for v in colorsys.hls_to_rgb(h, l, s))
            for h, s, l in zip(hues, sat, lit)]
    table = dict(zip(clan_ids, cols))
    if path:
        pd.DataFrame({"clan_acc": clan_ids, "colour": cols}).to_csv(path, index=False)
        logger.info(f"Wrote {len(table)} clan colours → {path} (reuse this to "
                    f"keep colours stable across figures)")
    return table


# ---------------------------------------------------------------------------
# Payload
# ---------------------------------------------------------------------------

def hidden_cluster_nodes(df, cuts, hide_top_n):
    """The n largest top-cut clusters plus their children at the next cut."""
    if hide_top_n <= 0:
        return set(), 0
    biggest = df[cuts[0]].value_counts().head(hide_top_n)
    inside = df[df[cuts[0]].isin(biggest.index)]
    kids = sorted(inside[cuts[1]].unique())
    nodes = {f"{c}{cuts[0]}" for c in biggest.index} | {f"{k}{cuts[1]}" for k in kids}
    logger.info(f"Hiding outlines of the {hide_top_n} largest {cuts[0]} "
                f"cluster(s) [{', '.join(map(str, biggest.index))}] plus their "
                f"{len(kids)} {cuts[1]} children; {biggest.sum():,} of "
                f"{len(df):,} motifs ({100 * biggest.sum() / len(df):.1f}%) sit "
                f"inside and are still drawn")
    return nodes, int(biggest.sum())


def build_payload(df, lay, cuts, clan_ids, clan_colours, hidden):
    n_cut = len(cuts)
    drawn = lay[(lay["depth"] >= 1) & (~lay["name"].isin(hidden))]
    clusters = drawn[drawn["depth"] <= n_cut]

    tips = (lay[lay["depth"] == n_cut + 1][["name", "x", "y", "r"]]
            .merge(df, left_on="name", right_on="accession", how="left"))
    if tips["accession"].isna().any():
        raise SystemExit("Some packed tips are missing from the cluster table.")

    for f in LABEL_FIELDS:
        if tips[f].astype(str).str.contains("|", regex=False).any():
            raise SystemExit(f"Column {f} contains '|', which the label format "
                             f"uses as its separator.")

    idx = {c: i for i, c in enumerate(clan_ids)}
    colour_idx = [idx.get(c, -1) if c else -1 for c in tips["clan_acc"]]

    labels = (tips[LABEL_FIELDS].astype(str)
              .apply(lambda r: "|".join(r), axis=1).tolist())

    # one row per tip, its cluster id at each cut, flattened row-major
    cid = tips[cuts].astype(int).to_numpy().reshape(-1).tolist()

    logger.info(f"Payload: {len(tips):,} tips, {len(clusters):,} drawn clusters "
                f"({len(drawn) - len(clusters):,} tip circles, "
                f"{len(lay) - len(drawn):,} nodes hidden)")

    return {
        "meta": {"cutNames": [cut_display_name(c) for c in cuts],
                 "nTips": int(len(tips)), "nClusters": int(len(clusters))},
        "palette": [clan_colours.get(c, "#c9c8c0") for c in clan_ids],
        "cl": {"x": clusters["x"].round(2).tolist(),
               "y": clusters["y"].round(2).tolist(),
               "r": clusters["r"].round(3).tolist(),
               "d": clusters["depth"].astype(int).tolist()},
        "tp": {"x": tips["x"].round(2).tolist(),
               "y": tips["y"].round(2).tolist(),
               "r": tips["r"].round(3).tolist(),
               "c": colour_idx},
        "lab": labels,
        "cid": cid,
    }


def render(template_path, payload, title, subtitle, out):
    template = Path(template_path).read_text()
    pieces = template.split("/*__DATA__*/null")
    if len(pieces) != 2:
        raise SystemExit(f"{template_path} must contain exactly one "
                         f"'/*__DATA__*/null' placeholder (found "
                         f"{len(pieces) - 1}).")
    # split rather than re.sub: the payload is megabytes and must not be read as
    # a regex replacement, where backslashes would be reinterpreted
    html = (pieces[0] + json.dumps(payload, separators=(",", ":")) + pieces[1]
            ).replace("__TITLE__", title).replace("__SUBTITLE__", subtitle)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(html)
    logger.info(f"{Path(out).name}: {payload['meta']['nTips']} tips, "
                f"{payload['meta']['nClusters']} clusters, "
                f"{Path(out).stat().st_size / 1024**2:.1f} MB")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Circle-packed view of a CAV motif dendrogram.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--clusters", required=True, help="Cluster CSV from cav_cluster.py.")
    ap.add_argument("--cuts", nargs="+", help="Cut columns, coarse to fine "
                                              "(default: the 5 coarsest).")
    ap.add_argument("--template", help="HTML template containing /*__DATA__*/null, "
                                       "__TITLE__ and __SUBTITLE__.")
    ap.add_argument("--out", help="Output HTML path (needs --template).")
    ap.add_argument("--payload-out", help="Also write the raw payload JSON here.")
    ap.add_argument("--layout-cache", help="Cache the packing here "
                                           "(default: alongside --clusters).")
    ap.add_argument("--clan-colours", help="clan_acc,colour CSV. Reused if present, "
                                           "written if not, so colours stay stable.")
    ap.add_argument("--hide-top-n", type=int, default=0,
                    help="Draw the n largest top-cut clusters, and their children "
                         "at the next cut, without outlines (default: 0 = all).")
    ap.add_argument("--title", default="Motif circle pack")
    args = ap.parse_args()

    if args.out and not args.template:
        ap.error("--out needs --template")
    if not args.out and not args.payload_out:
        ap.error("nothing to write: give --out (with --template) or --payload-out")

    df = pd.read_csv(args.clusters, dtype=str).fillna("")
    cut_cols = [c for c in df.columns if c.startswith("cut_")]
    cuts = args.cuts or cut_cols[:5]
    missing = [c for c in cuts if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing columns: {missing}. Available: {cut_cols}")
    for f in LABEL_FIELDS:
        if f not in df.columns:
            raise SystemExit(f"--clusters is missing the {f} column.")

    check_nesting(df, cuts)

    cache = args.layout_cache or str(
        Path(args.clusters).with_suffix("").as_posix()
        + "_layout_" + cuts[0].replace("cut_", "") + ".csv")
    lay = pack(df, cuts, cache)

    clan_ids = sorted(c for c in df["clan_acc"].unique() if c)
    colours = clan_palette(clan_ids, args.clan_colours)
    hidden, n_hidden_motifs = hidden_cluster_nodes(df, cuts, args.hide_top_n)

    payload = build_payload(df, lay, cuts, clan_ids, colours, hidden)

    if args.payload_out:
        Path(args.payload_out).write_text(json.dumps(payload, separators=(",", ":")))
        logger.info(f"Wrote payload → {args.payload_out}")

    if args.out:
        names = payload["meta"]["cutNames"]
        subtitle = (f"{len(df):,} motifs &middot; {len(clan_ids)} clans &middot; "
                    f"nested at {', '.join(names)}")
        if args.hide_top_n > 0:
            subtitle += (f" &middot; {args.hide_top_n} largest {names[0]} "
                         f"cluster{'s' if args.hide_top_n > 1 else ''} "
                         f"({n_hidden_motifs:,} motifs) drawn without outlines")
        render(args.template, payload, args.title, subtitle, args.out)

    print(df[list(cuts)].nunique().to_string())


if __name__ == "__main__":
    main()
