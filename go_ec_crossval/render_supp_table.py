#!/usr/bin/env python3
"""
render_supp_table.py — Render supp_table_rank_vs_training_size.csv as a
clean table image for direct inclusion in the supplement.

Usage
-----
python render_supp_table.py --csv results/supp_table_rank_vs_training_size.csv
"""

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


COL_LABELS = {
    "Rank bucket": "Rank of annotated\nEC partner",
    "n_pairs": "N pairs",
    "min": "Min",
    "q1": "Q1",
    "median": "Median",
    "q3": "Q3",
    "max": "Max",
    "mean": "Mean",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--out", default=None,
                        help="Output PNG path (default: <csv dir>/supp_table_rank_vs_training_size.png)")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    out_path = args.out or str(csv_path.with_suffix(".png"))

    df = pd.read_csv(csv_path)
    for col in ["min", "median", "max"]:
        df[col] = df[col].map(lambda v: f"{v:g}")

    col_keys   = list(df.columns)
    col_labels = [COL_LABELS.get(c, c) for c in col_keys]

    fig, ax = plt.subplots(figsize=(7.5, 0.5 + 0.4 * len(df)))
    ax.axis("off")

    n_cols = len(col_keys)
    first_col_w = 0.22
    other_w = (1 - first_col_w) / (n_cols - 1)
    col_widths = [first_col_w] + [other_w] * (n_cols - 1)

    tbl = ax.table(
        cellText=df.values,
        colLabels=col_labels,
        cellLoc="center",
        loc="center",
        colWidths=col_widths,
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.8)

    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if row == 0:
            cell.set_facecolor("#e8e8e8")
            cell.set_text_props(weight="bold")
        else:
            cell.set_facecolor("white")

    plt.title(
        "Supplemental Table. Annotated EC partner's training-set size\n"
        "(n_positive), by rank in the GO↔EC cosine validation",
        fontsize=10, pad=14,
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
