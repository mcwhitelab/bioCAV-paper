#!/usr/bin/env python
"""
GO callout candidates on the C=1e-4 library (2026-10-01), with the rules recovered from the
old-library callouts (both reproduce the old lists exactly on the data they were picked from):
  strict   CAV rank 1 (`cav_rank_ont_filt_pf`) and LLR > 0, DeepGo-SE emitted nothing for the
           term, and no DIAMOND hit to the term's training positives. Gives the old
           callouts/mf_callouts.csv (11) from the old live data.
  cav_only CAV rank 1 and LLR > 0, DeepGo-SE not in its top 10 (or not predicted), BLAST mid
           rank not in its top 10 (or no hit). Gives the old bp_cav_only_rank1_cases.csv (62)
           from the obsolete-included data it was picked on.
Inputs: figures/figure_data/go_rank_families_{ont}.csv (C=1e-4, live).
Writes strict_candidates.csv, cav_only_candidates.csv and jack_queries.txt (the new strict
candidates plus the 60 old strict candidates, for calibrating the jackhmmer recipe).
"""
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
OLD_JACK = Path("/xdisk/clairemcwhite/clairemcwhite/go_train_positives_embed/jackhmmer")
names, cur = {}, None
for l in open(PAPER / "go_eval/static/go_annotations/go.obo"):
    if l.startswith("id: GO:"): cur = l[4:].strip()
    elif l.startswith("name: ") and cur: names[cur] = l[6:].strip(); cur = None

strict, cav_only = [], []
for o in ["mf", "bp", "cc"]:
    d = pd.read_csv(PAPER / f"figures/figure_data/go_rank_families_{o}.csv").drop_duplicates(["protein_id", "go_term"])
    d["ont"], d["name"] = o, d.go_term.map(names)
    c1 = (d.cav_rank_ont_filt_pf == 1) & (d.llr > 0)
    strict.append(d[c1 & (d.tool_predicted == False) & (~d.blast_predicted)])
    cav_only.append(d[c1 & ~((d.tool_rank_ont_filt_pf <= 10) & d.tool_predicted)
                         & ~((d.blast_rank_mid <= 10) & d.blast_predicted)])
strict, cav_only = pd.concat(strict), pd.concat(cav_only)
strict.to_csv(HERE / "strict_candidates.csv", index=False)
cav_only.to_csv(HERE / "cav_only_candidates.csv", index=False)
old = pd.read_csv(OLD_JACK / "strict_candidates.csv")
q = sorted(set(strict.protein_id) | set(old.protein_id))
(HERE / "jack_queries.txt").write_text("\n".join(q) + "\n")
print("strict", strict.ont.value_counts().to_dict(), "| cav_only", cav_only.ont.value_counts().to_dict(),
      "| jackhmmer queries", len(q), f"({len(set(old.protein_id))} old)")
