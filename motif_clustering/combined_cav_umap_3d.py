#!/usr/bin/env python3
"""
combined_cav_umap_3d.py — One combined 2D or 3D UMAP of all EC + GO
(MF/BP/CC) CAV concept vectors in a single embedding, colored by category,
with hover text showing each concept's plain-English name (not just its
numeric ID).

Unlike cav_viz_pattern.py (one category per plot), this loads all four
categories into one shared embedding space so their relative arrangement
is directly comparable, then colors points by category rather than by
clan/Pfam annotation.

Usage
-----
python combined_cav_umap_3d.py \\
    --ec-pattern "$PROFAB_EC_SETS/ec_dataset_part*/ecNo_*/random_positive_train_max1000_cav/concept_v1.npy" \\
    --go-pattern "$PROFAB_GO_SETS/go_dataset_part*/GO_*/random_positive_train_max1000_cav/concept_v1.npy" \\
    --mf-list go_eval/results/go_namespace_lists/mf_ids.txt \\
    --bp-list go_eval/results/go_namespace_lists/bp_ids.txt \\
    --cc-list go_eval/results/go_namespace_lists/cc_ids.txt \\
    --ec-enzyme-dat  ec_eval/static/ec_annotations/enzyme.dat \\
    --ec-enzclass    ec_eval/static/ec_annotations/enzclass.txt \\
    --go-obo         go_eval/static/go_annotations/go.obo \\
    --dims 3 \\
    --out    go_eval/results/figures/combined_umap_3d.html
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cav_viz_pattern import load_directions_from_pattern, _embed, _axis_labels

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ec_eval" / "scripts"))
from ec_names import load_ec_names

CATEGORY_COLORS = {
    "EC": "#D55E00",
    "MF": "#0072B2",
    "BP": "#009E73",
    "CC": "#CC79A7",
}


def load_go_names(go_obo: str) -> dict:
    names = {}
    cur_id = None
    with open(go_obo) as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur_id = None
            elif line.startswith("id: "):
                cur_id = line[4:].strip()
            elif line.startswith("name: ") and cur_id:
                names[cur_id] = line[6:].strip()
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ec-pattern", required=True)
    parser.add_argument("--go-pattern", required=True)
    parser.add_argument("--mf-list", required=True)
    parser.add_argument("--bp-list", required=True)
    parser.add_argument("--cc-list", required=True)
    parser.add_argument("--ec-enzyme-dat", required=True)
    parser.add_argument("--ec-enzclass", required=True)
    parser.add_argument("--go-obo", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--coords-csv-out", default=None,
                        help="Optional: also write category, concept, display_name, "
                             "D1, D2[, D3] to this CSV path (e.g. for a static "
                             "ggplot render of the same embedding).")
    parser.add_argument("--dims", type=int, default=3, choices=[2, 3])
    parser.add_argument("--reducer", default="umap", choices=["umap", "tsne", "pca"])
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.1)
    args = parser.parse_args()

    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("plotly is not installed. Run: pip install plotly")

    ec_dirs = load_directions_from_pattern(args.ec_pattern)
    go_dirs_all = load_directions_from_pattern(args.go_pattern)

    def restrict(d, list_path):
        keep = {line.strip() for line in open(list_path) if line.strip()}
        return {k: v for k, v in d.items() if k in keep}

    mf_dirs = restrict(go_dirs_all, args.mf_list)
    bp_dirs = restrict(go_dirs_all, args.bp_list)
    cc_dirs = restrict(go_dirs_all, args.cc_list)

    print(f"EC: {len(ec_dirs)}  MF: {len(mf_dirs)}  BP: {len(bp_dirs)}  CC: {len(cc_dirs)}")

    ec_names  = load_ec_names(args.ec_enzyme_dat, args.ec_enzclass)
    go_names  = load_go_names(args.go_obo)
    print(f"Name lookups: {len(ec_names)} EC, {len(go_names)} GO")

    all_dirs   = {}
    categories = {}
    for cat, d in [("EC", ec_dirs), ("MF", mf_dirs), ("BP", bp_dirs), ("CC", cc_dirs)]:
        for name, vec in d.items():
            key = f"{cat}:{name}"   # disambiguate in case of name collisions across categories
            all_dirs[key] = vec
            categories[key] = cat

    print(f"Total combined: {len(all_dirs)} CAV vectors")

    coord_df, red = _embed(all_dirs, args.reducer, dims=args.dims,
                           n_neighbors=args.n_neighbors, min_dist=args.min_dist)
    xl, yl = _axis_labels(args.reducer, red)

    def normalize_concept(cat, raw_name):
        """load_directions_from_pattern returns raw directory names
        (e.g. 'ecNo_1-1-1-1', 'GO_0004022') -- convert to the clean
        lookup-key format ('1.1.1.1', 'GO:0004022')."""
        if cat == "EC":
            return raw_name.replace("ecNo_", "", 1).replace("-", ".")
        return raw_name.replace("GO_", "GO:", 1)

    cats = [categories[k] for k in coord_df.index]
    raw_names     = [k.split(":", 1)[1] for k in coord_df.index]
    concept_names = [normalize_concept(cat, raw) for cat, raw in zip(cats, raw_names)]

    def display_name(cat, concept):
        if cat == "EC":
            return ec_names.get(concept, "")
        return go_names.get(concept, "")

    hover_texts = []
    for concept, cat in zip(concept_names, cats):
        label = "EC " + concept if cat == "EC" else concept
        desc = display_name(cat, concept)
        text = f"<b>{label}</b><br>{desc}<br>Category: {cat}" if desc else \
               f"<b>{label}</b><br>Category: {cat}"
        hover_texts.append(text)

    fig = go.Figure()
    for cat in ["EC", "MF", "BP", "CC"]:
        mask = [c == cat for c in cats]
        idx = coord_df.index[mask]
        trace_kwargs = dict(
            mode="markers",
            name=f"{cat} (n={sum(mask)})",
            marker=dict(size=3 if args.dims == 3 else 5,
                       color=CATEGORY_COLORS[cat], opacity=0.6, line=dict(width=0)),
            text=[hover_texts[i] for i, m in enumerate(mask) if m],
            hovertemplate="%{text}<extra></extra>",
        )
        if args.dims == 3:
            fig.add_trace(go.Scatter3d(
                x=coord_df.loc[idx, "D1"], y=coord_df.loc[idx, "D2"], z=coord_df.loc[idx, "D3"],
                **trace_kwargs,
            ))
        else:
            fig.add_trace(go.Scatter(
                x=coord_df.loc[idx, "D1"], y=coord_df.loc[idx, "D2"],
                **trace_kwargs,
            ))

    title = f"Combined CAV space ({args.reducer.upper()}, {args.dims}D) — EC + GO (MF/BP/CC), n={len(all_dirs)}"
    if args.dims == 3:
        zl = yl.replace("2", "3")
        layout = dict(
            scene=dict(
                xaxis_title=xl, yaxis_title=yl, zaxis_title=zl,
                xaxis=dict(showgrid=True, zeroline=False),
                yaxis=dict(showgrid=True, zeroline=False),
                zaxis=dict(showgrid=True, zeroline=False),
            ),
        )
    else:
        layout = dict(
            xaxis=dict(title=xl, showgrid=False, zeroline=False),
            yaxis=dict(title=yl, showgrid=False, zeroline=False),
            plot_bgcolor="white",
        )

    fig.update_layout(
        title=dict(text=title, font_size=15),
        legend=dict(itemsizing="constant"),
        hoverlabel=dict(bgcolor="white", font_size=12),
        paper_bgcolor="white",
        width=1000, height=800,
        **layout,
    )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(args.out, include_plotlyjs="cdn")
    print(f"Saved {args.out}")

    if args.coords_csv_out:
        import pandas as pd
        out_df = coord_df.copy()
        out_df["category"] = cats
        out_df["concept"] = concept_names
        out_df["display_name"] = [display_name(c, n) for c, n in zip(cats, concept_names)]
        Path(args.coords_csv_out).parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(args.coords_csv_out, index=False)
        print(f"Saved {args.coords_csv_out}")


if __name__ == "__main__":
    main()
