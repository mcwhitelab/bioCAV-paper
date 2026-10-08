#!/usr/bin/env python
"""
Screen of sequence-family cutoffs for preview_C_by_families (exploratory, 2026-09-29).

Pairs: known-positive-excluded GO validation pairs, strict (_pf) candidate filter,
gated ranks as in panel C; BLAST = tie-averaged DIAMOND nearest-neighbour rank,
"no prediction" = no hit to the true term's training positives.
n_families = sequence families among the term's CAV training positives
(../diversity/term_diversity.py).

Writes
  family_screen_bins.csv    top-10 / Rank-1 per non-cumulative family-count bin
  family_screen_cutoff.csv  top-10 / Rank-1 below vs at-or-above each cutoff K
95% intervals: bootstrap over GO terms (pairs within a term are not independent),
1,000 resamples, within ontology x group.
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1] / "figures" / "figure_data_staging" / "noknown"
DIV = HERE.parent / "diversity"
ONTS = ["mf", "bp", "cc"]
BINS = [1, 5, 10, 15, 20, 30, 40, 60, 100, 1000]
BIN_LABELS = ["1-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40-59", "60-99", ">=100"]
CUTOFFS = [5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100]
NBOOT = 1000


def load():
    frames = []
    for o in ONTS:
        r = pd.read_csv(STAGE / f"go_specificity_ranks_allterms_{o}.csv") \
            .merge(pd.read_csv(HERE / f"blast_ranks_{o}.csv"), on=["protein_id", "go_term"]) \
            .merge(pd.read_csv(DIV / f"term_diversity_{o}.csv"), on="go_term")
        r["ontology"] = o.upper()
        frames.append(r)
    r = pd.concat(frames, ignore_index=True)
    cp, tp, bp = r.llr > 0, r.tool_predicted.astype(bool), r.blast_predicted.astype(bool)
    out = pd.DataFrame({"ontology": r.ontology, "go_term": r.go_term, "n_families": r.n_families})
    for m, rank, pred in [("CAV", r.cav_rank_ont_filt_pf, cp),
                          ("DeepGoSE", r.tool_rank_ont_filt_pf, tp),
                          ("BLAST", r.blast_rank_mid, bp)]:
        out[f"{m}_top10"] = ((rank <= 10) & pred).astype(float)
        out[f"{m}_rank1"] = ((rank < 1.5) & pred).astype(float)
    return out


def summarise(d, rng):
    """Pair-level rates with a term-level bootstrap interval."""
    cols = [c for c in d.columns if c.endswith(("_top10", "_rank1"))]
    per_term = d.groupby("go_term")[cols].agg(["sum", "count"])
    sums = np.stack([per_term[(c, "sum")].values for c in cols], axis=1)
    cnts = per_term[(cols[0], "count")].values
    est = sums.sum(0) / cnts.sum()
    nt = len(cnts)
    idx = rng.integers(0, nt, size=(NBOOT, nt))
    boot = sums[idx].sum(1) / cnts[idx].sum(1)[:, None]
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    rows = []
    for i, c in enumerate(cols):
        m, metric = c.split("_")
        rows.append({"method": m, "metric": metric, "rate": est[i], "lo": lo[i], "hi": hi[i],
                     "n_pairs": int(cnts.sum()), "n_terms": nt})
    return rows


def main():
    rng = np.random.default_rng(42)
    d = load()
    d = pd.concat([d, d.assign(ontology="All")], ignore_index=True)
    d["bin"] = pd.cut(d.n_families, BINS, right=False, labels=BIN_LABELS)

    bins = []
    for (o, b), g in d.groupby(["ontology", "bin"], observed=True):
        for row in summarise(g, rng):
            bins.append({"ontology": o, "bin": b, **row})
    pd.DataFrame(bins).to_csv(HERE / "family_screen_bins.csv", index=False, float_format="%.4f")

    cut = []
    for o, g in d.groupby("ontology"):
        for k in CUTOFFS:
            for side, sub in [("below", g[g.n_families < k]), ("at_or_above", g[g.n_families >= k])]:
                if sub.go_term.nunique() < 3:
                    continue
                for row in summarise(sub, rng):
                    cut.append({"ontology": o, "cutoff": k, "side": side, **row})
    pd.DataFrame(cut).to_csv(HERE / "family_screen_cutoff.csv", index=False, float_format="%.4f")

    b = pd.DataFrame(bins)
    t = b[(b.metric == "top10")].pivot_table(index=["ontology", "bin", "n_terms", "n_pairs"],
                                             columns="method", values="rate", observed=True)
    pd.set_option("display.width", 200)
    print((100 * t).round(1).to_string())


if __name__ == "__main__":
    main()
