#!/usr/bin/env python
"""
Export Fig 2C data (2026-09-29): rank tiers of CAV, DeepGo-SE and a DIAMOND
nearest-neighbour baseline, split by the number of sequence families among the
true term's CAV training positives.

Inputs (known-positive-excluded pairs, strict _pf candidate filter, obsolete GO
terms excluded):
  noobsolete/go_specificity_ranks_allterms_{ont}.csv  CAV / DeepGo-SE ranks
  noobsolete/blast_ranks_{ont}.csv                    BLAST ranks (blast_ranks.py)
  ../diversity/term_diversity_{ont}.csv               n_families (term_diversity.py)
Output: figures/figure_data/go_rank_families_{ont}.csv, one row per pair.
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1] / "figures" / "figure_data"

for ont in ["mf", "bp", "cc"]:
    r = pd.read_csv(HERE / "noobsolete" / f"go_specificity_ranks_allterms_{ont}.csv")
    b = pd.read_csv(HERE / "noobsolete" / f"blast_ranks_{ont}.csv")
    d = pd.read_csv(HERE.parent / "diversity" / f"term_diversity_{ont}.csv")
    # one exact duplicate pair per ontology upstream (two annotation IDs that
    # normalise to the same accession, e.g. P0DTC2 and a -PRO_ chain); drop it
    r = r.drop_duplicates(["protein_id", "go_term"])
    b = b.drop_duplicates(["protein_id", "go_term"])
    m = r.merge(b, on=["protein_id", "go_term"], validate="one_to_one") \
         .merge(d[["go_term", "n_train", "n_families"]], on="go_term", validate="many_to_one")
    assert len(m) == len(r) == len(b)
    out = m[["protein_id", "go_term", "n_train", "n_families", "n_terms_ont",
             "llr", "cav_rank_ont_filt_pf",
             "tool_predicted", "tool_rank_ont_filt_pf",
             "blast_predicted", "blast_rank_mid", "blast_true_pident"]]
    out.to_csv(OUT / f"go_rank_families_{ont}.csv", index=False)
    print(f"{ont}: {len(out)} pairs, {out.go_term.nunique()} terms -> go_rank_families_{ont}.csv")
