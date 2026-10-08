#!/usr/bin/env python
"""Summarise the GO/EC C-scan: per arm, median per-term AUROC/AUPR (GO by ontology, EC by level),
median rank of the true term (+ % rank 1, % top 10), and median cos to the shipped CAV.
Writes summary_cscan.tsv and prints it."""
import pandas as pd
from pathlib import Path
H = Path(__file__).resolve().parent
t = pd.read_csv(H / "per_term_cscan.tsv", sep="\t"); p = pd.read_csv(H / "per_pair_cscan.tsv", sep="\t")
f = pd.read_csv(H / "fits_cscan.tsv", sep="\t")
arms = list(dict.fromkeys(t.arm))
rows = []
for (lib, ont), g in t.groupby(["library", "ontology"]):
    for arm in arms:
        a = g[g.arm == arm]; q = p[(p.library == lib) & (p.ontology == ont) & (p.arm == arm)]
        cos = f[(f.library == lib) & (f.ontology == ont)].get(f"cos_{arm}_vs_shipped")
        rows.append({"library": lib, "set": ont, "arm": arm, "n_terms": a.term.nunique(),
                     "auroc_med": a.auroc.median(), "aupr_med": a.aupr.median(),
                     "rank_med": q["rank"].median(), "rank1_pct": 100 * (q["rank"] == 1).mean() if q["rank"].notna().any() else float("nan"),
                     "top10_pct": 100 * (q["rank"] <= 10).mean() if q["rank"].notna().any() else float("nan"),
                     "cos_to_shipped_med": cos.median() if cos is not None else (1.0 if arm == "shipped" else float("nan"))})
s = pd.DataFrame(rows); s.to_csv(H / "summary_cscan.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(s.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
