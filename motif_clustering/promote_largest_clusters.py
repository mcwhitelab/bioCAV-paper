#!/usr/bin/env python3
"""
Rewrite a circle-pack input CSV so the N largest top-level clusters no longer
exist: their families are re-assigned one cut level up before any layout is
computed.

For a family inside one of the N largest clusters at the first cut column, every
cut column takes the value of the next one down (cut[0] <- cut[1], cut[1] <-
cut[2], ...); the last column keeps its own value, so the family still has one
entry per level.  Shifted values are offset by 1,000,000 so they cannot collide
with a real cluster number at the column they move into.  Families outside the
N largest clusters are untouched.

The result is an ordinary clustering table: feed it to the standard circle-pack
script with no drop/hide (hide_top_n <- 0) and the large clusters are never
part of the hierarchy the layout sees.

Usage: promote_largest_clusters.py IN.csv OUT.csv [--top-n 3]
       [--cuts cut_k00314 cut_k02769 cut_k04629 cut_k09172 cut_k18984]
"""
import argparse
import pandas as pd

OFFSET = 1_000_000

ap = argparse.ArgumentParser()
ap.add_argument("infile"); ap.add_argument("outfile")
ap.add_argument("--top-n", type=int, default=3)
ap.add_argument("--cuts", nargs="+",
                default=["cut_k00314", "cut_k02769", "cut_k04629", "cut_k09172", "cut_k18984"])
a = ap.parse_args()

df = pd.read_csv(a.infile)
top = a.cuts[0]
sizes = df[top].value_counts()
biggest = sizes.index[:a.top_n]
m = df[top].isin(biggest)
orig = df.loc[m, a.cuts].copy()
for i in range(len(a.cuts) - 1):
    df.loc[m, a.cuts[i]] = orig[a.cuts[i + 1]].values + OFFSET
df.loc[m, a.cuts[-1]] = orig[a.cuts[-1]].values + OFFSET

new_sizes = df[top].value_counts()
assert not df[top].isin(biggest).any(), "a largest cluster survived"
assert len(new_sizes) == len(sizes) - a.top_n + orig[a.cuts[1]].nunique()
df.to_csv(a.outfile, index=False)
print(f"{a.infile}: promoted {int(m.sum())} families out of {a.top_n} largest {top} clusters "
      f"(sizes {list(sizes.values[:a.top_n])}); top level now {len(new_sizes)} clusters, "
      f"largest {int(new_sizes.iloc[0])} families -> {a.outfile}")
