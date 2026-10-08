#!/usr/bin/env python
"""Split the Fig 2B CC pairs by what kind of CC term the true term is, and compare CAV vs BLAST
vs DeepGo-SE per class (same tiers as figures.R: rank < 1.5 = Rank 1; ungated -> no prediction).
Classes (GO is_a/part_of closure, checked in this order):
  complex        is_a GO:0032991 protein-containing complex
  extracellular  is_a/part_of GO:0005576 extracellular region, or is_a GO:0005618 cell wall
  membrane       is_a GO:0016020 membrane
  other compartment  everything else (organelles, lumens, cytoskeleton, ...)"""
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent; PAPER = HERE.parents[2]
obo = PAPER / "go_eval/static/go_annotations/go.obo"
par, name, cur = {}, {}, None
for line in open(obo):
    line = line.strip()
    if line == "[Term]": cur = None
    elif line.startswith("id: GO:"): cur = line[4:]; par[cur] = set()
    elif cur and line.startswith("name: "): name[cur] = line[6:]
    elif cur and line.startswith("is_a: "): par[cur].add(("is_a", line[6:16]))
    elif cur and line.startswith("relationship: part_of "): par[cur].add(("part_of", line[22:32]))
    elif line.startswith("[Typedef]"): cur = None
def anc(t, rels):
    seen, st = set(), [t]
    while st:
        x = st.pop()
        for r, p in par.get(x, ()):
            if r in rels and p not in seen: seen.add(p); st.append(p)
    return seen
def cls(t):
    if "GO:0032991" in anc(t, {"is_a"}) | {t}: return "complex"
    a = anc(t, {"is_a", "part_of"}) | {t}
    if "GO:0005576" in a or "GO:0005618" in anc(t, {"is_a"}) | {t}: return "extracellular / cell wall"
    if "GO:0016020" in anc(t, {"is_a"}) | {t}: return "membrane"
    return "other compartment"
d = pd.read_csv(PAPER / "figures/figure_data/go_rank_families_cc.csv")
d["term_class"] = d.go_term.map(cls); d["go_name"] = d.go_term.map(name)
d["fam"] = pd.cut(d.n_families, [-1, 19, 59, 1e9], labels=["<20", "20-59", ">=60"])
def r1(rank, ok): return ((rank < 1.5) & ok).mean() * 100
def t10(rank, ok): return ((rank <= 10) & ok).mean() * 100
def summ(g):
    return pd.Series({"pairs": len(g), "terms": g.go_term.nunique(),
        "CAV_r1": r1(g.cav_rank_ont_filt_pf, g.llr > 0), "BLAST_r1": r1(g.blast_rank_mid, g.blast_predicted == True),
        "DGS_r1": r1(g.tool_rank_ont_filt_pf, g.tool_predicted == True),
        "CAV_t10": t10(g.cav_rank_ont_filt_pf, g.llr > 0), "BLAST_t10": t10(g.blast_rank_mid, g.blast_predicted == True),
        "DGS_t10": t10(g.tool_rank_ont_filt_pf, g.tool_predicted == True),
        "CAVwin_BLASTlow": (((g.cav_rank_ont_filt_pf < 1.5) & (g.llr > 0)) & ~((g.blast_rank_mid <= 10) & (g.blast_predicted == True))).sum(),
        "BLASTwin_CAVlow": (((g.blast_rank_mid < 1.5) & (g.blast_predicted == True)) & ~((g.cav_rank_ont_filt_pf <= 10) & (g.llr > 0))).sum()})
pd.set_option("display.width", 200)
a = d.groupby("term_class").apply(summ).round(1); print(a.to_string()); a.to_csv(HERE / "cc_class_summary.csv")
b = d.groupby(["term_class", "fam"], observed=True).apply(summ).round(1); print(b.to_string()); b.to_csv(HERE / "cc_class_by_family.csv")
t = d.groupby(["term_class", "go_term", "go_name"]).apply(summ).round(1).reset_index().sort_values("pairs", ascending=False)
t.to_csv(HERE / "cc_term_summary.csv", index=False)
d.to_csv(HERE / "cc_pairs_with_class.csv", index=False)
