#!/usr/bin/env python3
"""
plot_go_ec_cosine_validation.py — Figure for go_ec_cosine_validation.py's
output: matched vs. background cosine-similarity distribution, and the
rank of each GO term's annotated EC partner among all trained EC CAVs.

Usage
-----
python plot_go_ec_cosine_validation.py --results-dir results
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True,
                        help="Dir containing go_ec_cosine_validation.csv and "
                             "go_ec_background_sims.npy (from go_ec_cosine_validation.py).")
    parser.add_argument("--out", default=None,
                        help="Output PNG path (default: <results-dir>/go_ec_cosine_validation.png).")
    parser.add_argument("--max-rank-show", type=int, default=20,
                        help="Ranks beyond this are grouped into a '>N' bin (default: 20).")
    parser.add_argument("--bg-sample", type=int, default=20000,
                        help="Max background points to sample for the histogram (default: 20000).")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    out_path = args.out or str(results_dir / "go_ec_cosine_validation.png")

    df = pd.read_csv(results_dir / "go_ec_cosine_validation.csv")
    bg = np.load(results_dir / "go_ec_background_sims.npy")

    rng = np.random.default_rng(0)
    bg_sample = rng.choice(bg, size=min(args.bg_sample, len(bg)), replace=False)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    ax = axes[0]
    bins = np.linspace(-0.3, 1.0, 60)
    ax.hist(bg_sample, bins=bins, alpha=0.6, density=True,
            label=f"Background (non-matching)\nn={len(bg):,} pairs", color="#999999")
    ax.hist(df["matched_sim"], bins=bins, alpha=0.7, density=True,
            label=f"Annotated GO↔EC match\nn={len(df)} pairs", color="#D55E00")
    ax.axvline(df["matched_sim"].median(), color="#D55E00", ls="--", lw=1)
    ax.axvline(np.median(bg), color="#666666", ls="--", lw=1)
    ax.set_xlabel("Cosine similarity between GO-CAV and EC-CAV")
    ax.set_ylabel("Density")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Concept-vector similarity:\nannotated match vs. background", fontsize=10)

    ax = axes[1]
    max_show = args.max_rank_show
    ranks = df["rank"].clip(upper=max_show + 1)
    counts = ranks.value_counts().sort_index()
    labels = [str(i) for i in range(1, max_show + 1)] + [f">{max_show}"]
    vals = [counts.get(i, 0) for i in range(1, max_show + 1)] + [counts.get(max_show + 1, 0)]
    ax.bar(range(len(vals)), vals, color="#0072B2")
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=6, rotation=90)
    n_ec = df["n_ec_cavs"].iloc[0]
    ax.set_xlabel(f"Rank of annotated EC partner\n(out of {n_ec} trained EC CAVs)")
    ax.set_ylabel("Number of GO terms")
    ax.spines[["top", "right"]].set_visible(False)
    pct1 = (df["rank"] == 1).mean() * 100
    ax.set_title(f"Rank 1: {pct1:.1f}% of {len(df)} GO terms", fontsize=10)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
