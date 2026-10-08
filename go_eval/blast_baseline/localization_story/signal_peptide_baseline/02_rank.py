#!/usr/bin/env python
"""Sequence-feature baseline for the Fig 2B CC ranking (preview, 2026-10-06).
Question: do the CAV wins on secreted / cell-envelope terms (plant-type cell wall, fungal-type
cell wall, cell outer membrane) need more than 'has a signal peptide, is from kingdom K'?

Per CC term t, from its CAV training positives (C1e-4_extra/term_members.tsv, same set BLAST
searches against): P(f | t) over a discrete feature profile f, add-0.5 smoothed.
A query is scored for every candidate term by P(f_query | t) (= the term's likelihood ratio
vs background, up to a per-query constant) and the true term is mid-ranked (ties averaged)
among exactly the Fig 2B candidates (same live terms, ancestor / GPAD / proFAB exclusions as
integrated_score_preview/01_integrated_ranks.py; n_cand asserted equal).
Feature profiles (variants):
  sp        signal peptide yes/no
  sp_king   signal peptide x kingdom (8 classes)
  sp_tm_king signal peptide x TM helices (0 / 1 / >=2) x kingdom
The baseline always makes a call (no 'no prediction' tier). Nothing is fitted to the eval pairs."""
import glob, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata
HERE = Path(__file__).resolve().parent; PAPER = HERE.parents[3]; GE = PAPER / "go_eval"; BB = GE / "blast_baseline"
STAGE = PAPER / "figures/figure_data_staging/C1e-4"; LIVE = PAPER / "figures/figure_data"
BASE = Path("/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay/GO_sets")
SUB = "random_positive_train_max1000_cav"
sys.path.insert(0, str(GE / "rank_all_terms"))
import os; os.environ["FIG_DATA_DIR"] = str(STAGE)
import compute_ranks_all_terms as C  # noqa: E402
C.GO_SETS = BASE
ns, parents = C.parse_obo(C.OBO)
dirs = sorted(glob.glob(str(BASE / "go_dataset_part*" / "GO_*" / SUB)))
go_ids = [Path(d).parent.name.replace("GO_", "GO:") for d in dirs]
col = {g: i for i, g in enumerate(go_ids)}
cav_ns = np.array([ns.get(g) for g in go_ids], dtype=object)
cc_mask = cav_ns == "cc"

gpad = pd.concat([pd.read_csv(GE / "static/go_annotations" / f"{o}.gpad", sep="\t", comment="!", header=None,
                              usecols=[0, 1, 3], names=["db", "pid", "go"]) for o in C.ONTS])
gpad = gpad[gpad.db == "UniProtKB"]; gpad["pid"] = gpad.pid.astype(str).str.split("-").str[0]
known_true = {}
for pid, gs in gpad.groupby("pid").go:
    s = set(gs)
    for t in list(s): s |= C.ancestors(t, parents)
    known_true[pid] = s
prs = pd.read_csv(STAGE / "go_specificity_ranks_cc.csv").drop_duplicates(["protein_id", "go_term"])
profab_pos = C.load_profab_positives(set(prs.protein_id))

feat = pd.read_csv(HERE / "features.csv").set_index("protein_id")
feat["tm"] = np.minimum(feat.n_tm, 2)
PROF = {"sp": ["sp"], "sp_king": ["sp", "kingdom"], "sp_tm_king": ["sp", "tm", "kingdom"]}
mem = pd.read_csv(BB / "C1e-4_extra/term_members.tsv", sep="\t", names=["go_term", "s"])
mem = mem[mem.go_term.isin(col) & (mem.go_term.map(ns) == "cc")]
mem = mem[mem.s.isin(feat.index)]
L = {}  # variant -> (profile index, matrix [n_terms x n_profiles] of P(f|t))
for v, cols in PROF.items():
    key = feat[cols].astype(str).agg("|".join, axis=1)
    cats = sorted(key.unique()); ci = {c: i for i, c in enumerate(cats)}
    M = np.full((len(go_ids), len(cats)), 0.5)
    k = mem.s.map(key).map(ci).values; t = mem.go_term.map(col).values
    np.add.at(M, (t, k), 1)
    M /= M.sum(1, keepdims=True)
    L[v] = (key.map(ci), np.log(M))

fam = pd.read_csv(LIVE / "go_rank_families_cc.csv")
integ = pd.read_csv(HERE.parents[1] / "integrated_score_preview/integrated_ranks_cc.csv")[["protein_id", "go_term", "n_cand"]]
rows = []
for r in prs.itertuples(index=False):
    p, g = r.protein_id, r.go_term; j = col[g]
    drop = (C.ancestors(g, parents) | known_true.get(p, set())) - {g}
    for t in profab_pos.get(p, ()): drop |= {t} | C.ancestors(t, parents)
    drop -= {g}
    keep = cc_mask.copy(); keep[[col[t] for t in drop if t in col]] = False
    cand = np.where(keep)[0]; jj = int(np.searchsorted(cand, j)); assert cand[jj] == j
    row = {"protein_id": p, "go_term": g, "n_cand_sp": len(cand), "has_features": p in feat.index}
    if p in feat.index:
        row.update(sp=bool(feat.at[p, "sp"]), n_tm=int(feat.at[p, "n_tm"]), kingdom=feat.at[p, "kingdom"])
        for v, (qi, LM) in L.items():
            s = LM[cand, qi[p]]
            row[f"{v}_rank"] = float(rankdata(-s, method="average")[jj])
    rows.append(row)
out = pd.DataFrame(rows).merge(integ, on=["protein_id", "go_term"], validate="one_to_one")
assert (out.n_cand_sp == out.n_cand).all(), "candidate sets differ from integrated_ranks_cc"
out = out.drop(columns="n_cand_sp").merge(fam, on=["protein_id", "go_term"], validate="one_to_one")
cls = pd.read_csv(HERE.parent / "cc_pairs_with_class.csv")[["protein_id", "go_term", "term_class", "go_name"]]
out = out.merge(cls, on=["protein_id", "go_term"], validate="one_to_one")
out.to_csv(HERE / "sp_baseline_ranks_cc.csv", index=False)
print(len(out), "pairs;", (~out.has_features).sum(), "without UniProt features (counted as misses)")
