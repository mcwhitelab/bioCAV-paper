#!/usr/bin/env python
"""Rank 1 / top 10 per CC term class and per secreted / cell-envelope term: CAV and BLAST
(Fig 2B tiers; ungated = miss) vs the three signal-peptide baselines (always a call; pairs
without UniProt features = miss). Extracellular class is also reported without apoplast
(GO:0048046), which is kept out of callouts (separate project)."""
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent
d = pd.read_csv(HERE / "sp_baseline_ranks_cc.csv")
V = ["sp", "sp_king", "sp_tm_king"]
def summ(g):
    s = {"pairs": len(g), "terms": g.go_term.nunique(), "frac_sp": round(g.sp.mean() * 100, 1)}
    hits = {"CAV": (g.cav_rank_ont_filt_pf, g.llr > 0), "BLAST": (g.blast_rank_mid, g.blast_predicted == True)}
    hits.update({v: (g[f"{v}_rank"].fillna(1e9), g.has_features) for v in V})
    for k, (r, ok) in hits.items():
        s[f"{k}_r1"] = round(((r < 1.5) & ok).mean() * 100, 1)
        s[f"{k}_t10"] = round(((r <= 10) & ok).mean() * 100, 1)
    return pd.Series(s)
pd.set_option("display.width", 250)
d["class2"] = d.term_class.where(d.go_term != "GO:0048046", "extracellular (apoplast only)")
out = []
for name, key in [("term_class", "term_class"), ("class_ex_apoplast", "class2")]:
    a = d.groupby(key).apply(summ); a.insert(0, "grouping", name); out.append(a.reset_index().rename(columns={key: "class"}))
res = pd.concat(out); res.to_csv(HERE / "sp_baseline_summary.csv", index=False); print(res.to_string(index=False))
ex = d[d.term_class.str.startswith("extra") | d.go_term.isin(["GO:0005778", "GO:0020002"])]
t = ex.groupby(["go_term", "go_name"]).apply(summ).reset_index(); t.to_csv(HERE / "sp_baseline_secreted_terms.csv", index=False)
print(t.to_string(index=False))
# within secreted plant proteins: is CAV beating the profile baseline by separating terms the profile can't?
pl = d[(d.kingdom == "Viridiplantae") & d.sp]
print("\nplant SP+ queries:", len(pl)); print(pl.groupby("term_class").apply(summ).to_string())
