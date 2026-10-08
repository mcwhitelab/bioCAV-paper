#!/usr/bin/env python3
"""
scramble_cluster_ids.py — Randomly permute cluster-ID labels in a cluster CSV.

Cluster IDs numbered sequentially along the leaf order make adjacent clusters
adjacent numbers, so an Excel colour scale renders them as a smooth gradient
and cluster breaks are invisible.  Permuting the labels gives neighbouring
clusters unrelated numbers, so conditional formatting shows a hard colour break
at every boundary.

This only relabels — cluster membership, row order and every other column are
untouched.  Note that permuted IDs are not comparable across cut levels (the
same group gets an unrelated number in each column); use the unscrambled file
for cross-resolution analysis.

Usage
-----
python scramble_cluster_ids.py --in clusters.csv --out clusters_scrambled.csv
python scramble_cluster_ids.py --in clusters.csv          # in place (--out omitted)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def scramble(df, seed=42, columns=None):
    cols = columns or [c for c in df.columns if c.startswith("cut_")]
    rng = np.random.default_rng(seed)
    out = df.copy()
    for col in cols:
        vals = out[col].astype(int).to_numpy()
        uniq = np.unique(vals)
        # 1-based labels, permuted
        id_map = dict(zip(uniq, rng.permutation(len(uniq)) + 1))
        out[col] = [id_map[v] for v in vals]
    return out, cols


def main():
    ap = argparse.ArgumentParser(
        description="Permute cluster-ID labels for visual contrast in Excel.")
    ap.add_argument("--in", dest="inp", required=True, help="Cluster CSV.")
    ap.add_argument("--out", help="Output CSV (default: overwrite input).")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--columns", nargs="+",
                    help="Columns to scramble (default: every cut_* column).")
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    out, cols = scramble(df, seed=args.seed, columns=args.columns)
    dest = Path(args.out or args.inp)
    out.to_csv(dest, index=False)
    print(f"Scrambled {len(cols)} column(s) in {len(out)} rows → {dest}")
    for c in cols:
        print(f"  {c}: {out[c].nunique()} clusters")


if __name__ == "__main__":
    main()
