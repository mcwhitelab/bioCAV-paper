#!/usr/bin/env python
"""
Fig 3A ranks with the C=1e-4 converged EC library (2026-10-01, preview only).

Runs ec_eval/blast_baseline/ec_ranks_noknown.py unchanged except that
  - the EC tree (CAV vectors, spans, known-positive sets) is --lib's: old = proFAB
    (952 shipped ECs), new = the C=1e-4 overlay (953: + EC 3.4.22.28, added 2026-10-01);
  - BLAST hits/members for new are extended to EC 3.4.22.28 (go_eval/rank_C1e-4/
    extend_blast.sh), so BLAST ranks over the same candidates; CLEAN columns unchanged;
  - the llr gate is recomputed with that library (Gaussian fits to the CAV scores of the
    EC's cached proFAB test positive/negative embeddings, shared scaler_v1.pkl), not
    carried from the shipped file;
  - everything is written to --out (figure_data/ and the noknown staging dir untouched).
  --lib old  regression check: must reproduce figure_data/ec_ranks_samelevel_noknown.csv
"""
import argparse, sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy import stats

PAPER = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper")
PROFAB = Path("/xdisk/clairemcwhite/proFAB/EC_sets")
LIBS = {"old": PROFAB,
        "new": Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay/EC_sets")}
BLAST = {"old": PAPER / "ec_eval" / "blast_baseline",
         "new": PAPER / "ec_eval" / "blast_baseline" / "C1e-4_extra"}
SUB = "random_positive_train_max1000_cav"

ap = argparse.ArgumentParser()
ap.add_argument("--lib", choices=LIBS, required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(PAPER / "ec_eval" / "blast_baseline"))
import ec_ranks_noknown as E  # noqa: E402


E.EC_SETS = LIBS[a.lib]
E.HERE = BLAST[a.lib]  # test_vs_ectrain.tsv, ec_term_members.tsv
E.STAGE = OUT
E.FIG_DATA = OUT

# llr gate with this library
scaler = joblib.load(E.SCALER)
pairs = pd.read_csv(PAPER / "figures" / "figure_data" / "ec_specificity_ranks.csv", dtype={"ec_number": str})
embs, sids = E.load_sequence_embeddings(E.BENCH / "test_ec.fasta.pkl")
idx = {}
for i, s in enumerate(sids):
    idx.setdefault(E.norm(s), i)
Xval = scaler.transform(embs)
ec_dir = {d.name[5:].replace("-", "."): d for d in LIBS[a.lib].glob("ec_dataset_part*/ecNo_*")}
llr = {}
for ec, grp in pairs.groupby("ec_number"):
    d = ec_dir[ec]
    w = np.load(d / SUB / "concept_v1.npy").ravel().astype(np.float64)
    sc = {}
    for k, f in [("pos", "random_positive_test_max1000.span.fasta.pkl"),
                 ("neg", "random_negative_test_max1000.span.fasta.pkl")]:
        if not (d / f).exists():
            sc = None; break
        e, _ = E.load_sequence_embeddings(d / f)
        sc[k] = scaler.transform(e) @ w
    for p in grp.protein_id:
        if sc is None:
            llr[(p, ec)] = np.nan; continue
        v = Xval[idx[p]] @ w
        llr[(p, ec)] = (stats.norm.logpdf(v, sc["pos"].mean(), sc["pos"].std()) -
                        stats.norm.logpdf(v, sc["neg"].mean(), sc["neg"].std()))
pairs["llr_live"] = pairs["llr"]
pairs["llr"] = [llr[(p, e)] for p, e in zip(pairs.protein_id, pairs.ec_number)]
ok = pairs.llr.notna()
print(f"llr computed for {ok.sum()}/{len(pairs)} pairs; max|diff| vs live {np.abs(pairs.llr - pairs.llr_live)[ok].max():.2e}; "
      f"gate agreement {((pairs.llr > 0) == (pairs.llr_live > 0))[ok].mean():.4f}")
pairs.to_csv(OUT / "ec_specificity_ranks.csv", index=False)

E.main()

# as export_fig3A.py
d = pd.read_csv(OUT / "ec_ranks_noknown.csv", dtype={"ec_number": str})
d = d[~d.known_positive]
# 2026-10-06: drop UniProt fragment entries (partial sequences, e.g. C0HLV0 = 27 aa of a ~10x
# larger endoglucanase); list from ec_eval/static/ec_annotations/EC-Bench/fetch_uniprot_fragments.py
frag = pd.read_csv(PAPER / "ec_eval/static/ec_annotations/EC-Bench/uniprot_fragments.tsv", sep="\t")
frag = set(frag.accession[frag.fragment])
n0 = len(d); d = d[~d.protein_id.isin(frag)]
print(f"dropped {n0 - len(d)} UniProt-fragment pairs ({n0} -> {len(d)})")
d[["protein_id", "ec_number", "level", "n_cand", "llr", "cav_rank_level_filt_pf",
   "blast_predicted", "blast_rank_mid", "blast_true_pident",
   "clean_supcon_correct", "ecbench_blastp_correct"]].to_csv(OUT / "ec_ranks_samelevel_noknown.csv", index=False)
print(f"{len(d)} leak-free pairs -> {OUT / 'ec_ranks_samelevel_noknown.csv'}")
