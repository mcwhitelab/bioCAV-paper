#!/usr/bin/env python
"""Sample the GO/EC terms for the lbfgs-convergence pilot (seed 42).

GO: 20 terms per ontology from the Fig 2 leak-free (noknown) evaluated set,
stratified by training n_positive: <100 / 100-499 / >=500, 7/7/6, topped up at random to 20.
No leak-free Fig 2 term has more than 297 training positives, so in practice the GO
sample is two strata (26 terms <100, 34 terms 100-297).
EC: every EC term in the Fig 3A pair list (ec_specificity_ranks_samelevel.csv, 49 terms).
Writes terms.tsv: library, ontology, term, span, cav_dir, n_positive, stratum.
"""
import json, glob
from pathlib import Path
import numpy as np, pandas as pd

P = Path("/groups/clairemcwhite/cav_workspace/bioCAV-paper")
GO = Path("/xdisk/clairemcwhite/proFAB/GO_sets"); EC = Path("/xdisk/clairemcwhite/proFAB/EC_sets")
rng = np.random.default_rng(42)

def npos(cav_dir):
    return json.load(open(cav_dir / "report_v1.json"))["concept_metrics"]["n_positive"]

def stratum(n):
    return "lt100" if n < 100 else ("100-499" if n < 500 else "ge500")

rows = []
for ont in ["mf", "bp", "cc"]:
    terms = pd.read_csv(P / f"go_eval/results/temporal_eval_{ont}_noknown/eval_temporal_per_term_summary.tsv", sep="\t").go_term
    recs = []
    for t in terms:
        d = glob.glob(str(GO / "go_dataset_part*" / t.replace(":", "_")))
        assert len(d) == 1, t
        cav = Path(d[0]) / "random_positive_train_max1000_cav"
        n = npos(cav); recs.append((t, Path(d[0]), cav, n, stratum(n)))
    df = pd.DataFrame(recs, columns=["term", "dir", "cav_dir", "n_positive", "stratum"])
    print(ont, len(df), df.stratum.value_counts().to_dict())
    per = {"lt100": 7, "100-499": 7, "ge500": 6}
    pick = []
    for s, k in per.items():
        sub = df[df.stratum == s]
        pick.append(sub.iloc[rng.choice(len(sub), min(k, len(sub)), replace=False)])
    pick = pd.concat(pick)
    if len(pick) < 20:   # top up from the other strata
        rest = df[~df.term.isin(pick.term)]
        pick = pd.concat([pick, rest.iloc[rng.choice(len(rest), 20 - len(pick), replace=False)]])
    for r in pick.itertuples():
        rows.append(("GO", ont, r.term, str(r.dir / "random_positive_train_max1000.span"), str(r.cav_dir), r.n_positive, r.stratum))

ec = pd.read_csv(P / "figures/figure_data/ec_specificity_ranks_samelevel.csv")
for e in sorted(ec.ec_number.unique()):
    d = glob.glob(str(EC / "ec_dataset_part*" / ("ecNo_" + e.replace(".", "-"))))
    assert len(d) == 1, e
    cav = Path(d[0]) / "random_positive_train_max1000_cav"
    n = npos(cav)
    rows.append(("EC", f"level{e.count('.') + 1}", e, str(Path(d[0]) / "random_positive_train_max1000.span"), str(cav), n, stratum(n)))

out = pd.DataFrame(rows, columns=["library", "ontology", "term", "span", "cav_dir", "n_positive", "stratum"])
out.to_csv("terms.tsv", sep="\t", index=False)
print(out.groupby(["library", "ontology"]).size().to_string())
print("total positives to embed:", out.n_positive.sum())
