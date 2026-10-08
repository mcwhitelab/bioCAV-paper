#!/usr/bin/env python
"""
Fig 2C ranks with the C=1e-4 converged GO library (2026-10-01, preview only).

Same pairs, candidates and filters as the live Fig 2C (known-positive-excluded pairs,
strict _pf filter, obsolete terms excluded); only the CAV vectors change:
  1. LLR gate recomputed with the given library: Gaussian fits to the CAV scores of the
     term's cached proFAB test positive/negative embeddings, scored at the validation
     protein (eval_go_temporal.py's definition, shared scaler_v1.pkl).
  2. compute_ranks_all_terms.main() with every CAV of the given library as candidates
     (old: the 6,677 shipped; new: the overlay's 6,678, i.e. + GO:0032436, added
     2026-10-01). The proFAB-positive (_pf) filter also reads the library's tree.
  3. BLAST ranks recomputed over the same candidates (blast_baseline/blast_ranks.py;
     new: hits/members extended to GO:0032436 by extend_blast.sh), then merged with the
     family counts as export_fig2C_families.py.

  --lib old  proFAB; regression check: must reproduce the live llr, allterms and BLAST ranks
  --lib new  the C=1e-4 overlay of proFAB (distribution/retrain_go_ec_C1e-4/04_build_overlay.py)
Writes only into --out (pairs file copied there first); figure_data/ is not touched.
"""
import argparse, glob, os, shutil, sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy import stats

PAPER = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper")
GE = PAPER / "go_eval"
PROFAB = Path("/xdisk/clairemcwhite/proFAB/GO_sets")
LIBS = {"old": PROFAB,
        "new": Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay/GO_sets")}
BLAST = {"old": GE / "blast_baseline", "new": GE / "blast_baseline" / "C1e-4_extra"}
PAIRS_SRC = PAPER / "figures" / "figure_data_staging" / "noknown"
SUB = "random_positive_train_max1000_cav"
ONTS = ["mf", "bp", "cc"]

ap = argparse.ArgumentParser()
ap.add_argument("--lib", choices=LIBS, required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
OUT = Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
os.environ["FIG_DATA_DIR"] = str(OUT)
sys.path.insert(0, str(GE / "rank_all_terms"))
import compute_ranks_all_terms as C  # noqa: E402  (reads FIG_DATA_DIR at import)
assert C.FIG_DATA == OUT
BASE = LIBS[a.lib]
C.GO_SETS = BASE  # load_profab_positives: spans of every CAV term of this library

TERMS = sorted(Path(d).parent.name for d in glob.glob(str(BASE / "go_dataset_part*" / "GO_*" / SUB)))
term_dir = {Path(d).name: Path(d) for d in glob.glob(str(BASE / "go_dataset_part*" / "GO_*"))}
cav_cache = {}


def cav(term_dirname):
    if term_dirname not in cav_cache:
        cav_cache[term_dirname] = np.load(term_dir[term_dirname] / SUB / "concept_v1.npy").ravel().astype(np.float64)
    return cav_cache[term_dirname]


def load_cavs():
    ids = [t.replace("GO_", "GO:") for t in TERMS]
    return ids, np.stack([cav(t) for t in TERMS]), ["shared"] * len(ids)


C.load_cavs = load_cavs
scaler = joblib.load(C.SCALER)

# 1. LLR gate
for ont in ONTS:
    pairs = pd.read_csv(PAIRS_SRC / f"go_specificity_ranks_{ont}.csv")
    embs, sids = C.load_sequence_embeddings(
        GE / "results" / f"temporal_eval_{ont}" / "eval_embeddings" / "val_proteins.span.fasta.pkl")
    idx = {}
    for i, s in enumerate(sids):
        idx[s] = i
        for p in s.split("|"):
            if p and p not in C.SKIP:
                idx.setdefault(p, i)
    Xval = scaler.transform(embs)
    llr = {}
    for g, grp in pairs.groupby("go_term"):
        t = g.replace("GO:", "GO_")
        w = cav(t)
        sc = {}
        for k, f in [("pos", "random_positive_test_max1000.span.fasta.pkl"),
                     ("neg", "random_negative_test_max1000.span.fasta.pkl")]:
            e, _ = C.load_sequence_embeddings(term_dir[t] / f)
            sc[k] = scaler.transform(e) @ w
        mp, sp, mn, sn = sc["pos"].mean(), sc["pos"].std(), sc["neg"].mean(), sc["neg"].std()
        for p in grp.protein_id:
            v = Xval[idx[p]] @ w
            llr[(p, g)] = stats.norm.logpdf(v, mp, sp) - stats.norm.logpdf(v, mn, sn)
    new = np.array([llr[(p, g)] for p, g in zip(pairs.protein_id, pairs.go_term)])
    d = np.abs(new - pairs.llr.values)
    print(f"[{ont}] llr vs live: max|diff| {d.max():.2e}; gate (llr>0) agreement "
          f"{((new > 0) == (pairs.llr.values > 0)).mean():.4f}; "
          f"pass rate live {(pairs.llr > 0).mean():.3f} -> {a.lib} {(new > 0).mean():.3f}", flush=True)
    pairs["llr_live"] = pairs["llr"]
    pairs["llr"] = new
    pairs.to_csv(OUT / f"go_specificity_ranks_{ont}.csv", index=False)

# 2. ranks over all live terms of the ontology
C.main()

# 3. BLAST ranks over the same candidates, then the Fig 2C table (as export_fig2C_families.py)
sys.path.insert(0, str(GE / "blast_baseline"))
import blast_ranks as Bl  # noqa: E402
(OUT / "blast").mkdir(exist_ok=True)
shutil.copy(BLAST[a.lib] / "term_members.tsv", OUT / "blast" / "term_members.tsv")
os.environ["HITS_FILE"] = str(BLAST[a.lib] / "val_vs_train.tsv")
Bl.HERE, Bl.STAGE = OUT / "blast", OUT
Bl.main()
for ont in ONTS:
    r = pd.read_csv(OUT / f"go_specificity_ranks_allterms_{ont}.csv").drop_duplicates(["protein_id", "go_term"])
    b = pd.read_csv(OUT / "blast" / f"blast_ranks_{ont}.csv").drop_duplicates(["protein_id", "go_term"])
    dv = pd.read_csv(GE / "diversity" / f"term_diversity_{ont}.csv")
    m = r.merge(b, on=["protein_id", "go_term"], validate="one_to_one") \
         .merge(dv[["go_term", "n_train", "n_families"]], on="go_term", validate="many_to_one")
    assert len(m) == len(r) == len(b)
    m[["protein_id", "go_term", "n_train", "n_families", "n_terms_ont", "llr", "cav_rank_ont_filt_pf",
       "tool_predicted", "tool_rank_ont_filt_pf", "blast_predicted", "blast_rank_mid",
       "blast_true_pident"]].to_csv(OUT / f"go_rank_families_{ont}.csv", index=False)
print("done ->", OUT)
