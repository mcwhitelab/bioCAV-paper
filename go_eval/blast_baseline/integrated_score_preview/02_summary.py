#!/usr/bin/env python
"""Rank-tier composition per ontology x family bin for CAV, DeepGo-SE, BLAST and the
integrated scores (01_integrated_ranks.py). Tier rules as figures.R tier_of()."""
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
TIERS = ["Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction"]
def tier(rank, pred):
    return np.select([~pred, rank < 1.5, rank <= 3, rank <= 10],
                     ["No prediction", "Rank 1", "Rank 2-3", "Rank 4-10"], "Rank >10")
long = []
for ont in ["mf", "bp", "cc"]:
    d = pd.read_csv(HERE / f"integrated_ranks_{ont}.csv")
    d["fam"] = pd.cut(d.n_families, [-np.inf, 20, 60, np.inf], right=False,
                      labels=["<20 families", "20-59 families", ">=60 families"])
    either = (d.llr > 0) | d.blast_predicted
    m = {"CAV": (d.cav_rank, d.llr > 0),
         "DeepGoSE": (d.tool_rank_ont_filt_pf, d.tool_predicted.astype(bool)),
         "BLAST": (d.blast_rank_mid, d.blast_predicted.astype(bool)),
         "CAV x BLAST (rank product)": (d.rankprod_rank, either),
         "CAV x BLAST (soft, eps 0.1)": (d["soft_eps0.1_rank"], either),
         "soft eps 0.01": (d["soft_eps0.01_rank"], either),
         "soft eps 1": (d["soft_eps1_rank"], either),
         "Best of CAV or BLAST, per pair": (d.oracle_min_rank, either)}
    for k, (r, p) in m.items():
        long.append(pd.DataFrame({"ontology": ont.upper(), "fam": d.fam.astype(str),
                                  "method": k, "tier": tier(r.values, p.values)}))
L = pd.concat(long)
comp = (L.groupby(["ontology", "fam", "method", "tier"]).size().rename("n").reset_index())
comp["prop"] = comp.n / comp.groupby(["ontology", "fam", "method"]).n.transform("sum")
comp.to_csv(HERE / "integrated_tier_composition.csv", index=False)
L["r1"] = L.tier == "Rank 1"; L["t10"] = L.tier.isin(TIERS[:3])
s = L.groupby(["ontology", "method"], sort=False).agg(n=("r1", "size"), rank1_pct=("r1", "mean"),
                                                       top10_pct=("t10", "mean")).reset_index()
s[["rank1_pct", "top10_pct"]] = (100 * s[["rank1_pct", "top10_pct"]]).round(1)
s.to_csv(HERE / "integrated_summary.csv", index=False)
print(s.to_string(index=False))
sf = L.groupby(["ontology", "fam", "method"], sort=False).agg(r1=("r1", "mean"), t10=("t10", "mean")).reset_index()
sf[["r1","t10"]] = (100*sf[["r1","t10"]]).round(1)
sf.to_csv(HERE / "integrated_summary_by_family.csv", index=False)
print(sf[sf.method.isin(["CAV","BLAST","CAV x BLAST (rank product)","CAV x BLAST (soft, eps 0.1)"])]
      .pivot_table(index=["ontology","fam"], columns="method", values="t10").to_string())
