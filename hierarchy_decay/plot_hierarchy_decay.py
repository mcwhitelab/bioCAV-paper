#!/usr/bin/env python3
"""
plot_hierarchy_decay.py — Supplemental figure: CAV projection score for a
handful of validation-positive proteins, scored against progressively
broader ancestor CAVs in the EC hierarchy. One line per protein; x-axis is
hierarchy level, labeled by EC notation depth (X, X.X, X.X.X, X.X.X.X).

Usage
-----
python plot_hierarchy_decay.py \\
    --ec-csv results/ec_hierarchy_decay.csv \\
    --out    results/hierarchy_decay.png
"""

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ec-csv", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ec_df = pd.read_csv(args.ec_csv)

    fig, ax = plt.subplots(figsize=(4.5, 4))

    for i, key in enumerate(ec_df["true_ec"].unique()):
        sub = ec_df[ec_df["true_ec"] == key].sort_values("level")
        legend_label = f"{sub['protein_id'].iloc[0]} ({key})"
        ax.plot(sub["level"], sub["score"], marker="o", markersize=4,
                color=PALETTE[i % len(PALETTE)], linewidth=1.4,
                alpha=0.85, label=legend_label)

    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(["X", "X.X", "X.X.X", "X.X.X.X"], fontsize=8)
    ax.set_xlabel("EC hierarchy level", fontsize=9)
    ax.set_ylabel("CAV projection score", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")

    plt.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
