#!/usr/bin/env python
"""
Export Fig 3A data (2026-09-29): leak-free EC level-3/4 pairs with CAV and
DIAMOND nearest-neighbour ranks and CLEAN-SupCon correctness.

Input: figures/figure_data_staging/noknown/ec_ranks_noknown.csv (ec_ranks_noknown.py).
Pairs whose protein is already a proFAB positive (train or test span) of the true
EC are dropped. Ranks use the strict candidate filter (same-level EC CAVs, minus
the protein's other EC-Bench ECs and its other proFAB-positive ECs at that level).
Output: figures/figure_data/ec_ranks_samelevel_noknown.csv
"""
from pathlib import Path

import pandas as pd

PAPER = Path(__file__).resolve().parents[2]
src = PAPER / "figures" / "figure_data_staging" / "noknown" / "ec_ranks_noknown.csv"
dst = PAPER / "figures" / "figure_data" / "ec_ranks_samelevel_noknown.csv"

d = pd.read_csv(src, dtype={"ec_number": str})
n_all = len(d)
assert not d.duplicated(["protein_id", "ec_number"]).any()
d = d[~d.known_positive]
out = d[["protein_id", "ec_number", "level", "n_cand", "llr", "cav_rank_level_filt_pf",
         "blast_predicted", "blast_rank_mid", "blast_true_pident",
         "clean_supcon_correct", "ecbench_blastp_correct"]]
out.to_csv(dst, index=False)
print(f"{n_all} pairs -> {len(out)} leak-free ({out.level.value_counts().sort_index().to_dict()}) -> {dst.name}")
