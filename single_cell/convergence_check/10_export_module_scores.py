#!/usr/bin/env python3
"""PREVIEW data for a compact Fig 6D replacement: per-cell module scores vs the L2 CAV
score, fibroblast/colorectum, paired-donor cells, current (stored) CAVs.
Modules = the named blocks of figures/previews/fib_coexpr_top20_stored.tsv
(ECM/quiescent fibroblast: the 9 annotated ECM genes of its 20-gene block; ADAMDEC1+ mucosal fibroblast, 6 genes; NSG1 excluded by
the membership rule). Score = mean over module genes of z-scored log1p(CP10K) across these cells.
CAVEAT (state in Methods): module genes were selected by |CAV r|, so score-vs-L2 correlation is
partly by construction; the non-circular results are the module co-expression structure and
the within-condition correlations.
Writes figures/figure_data_staging/singlecell_modules/{sc_fib_module_scores.csv, sc_fib_module_stats.tsv}"""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from pathlib import Path
from scipy.stats import pearsonr
LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
OUT = Path("../../figures/figure_data_staging/singlecell_modules"); OUT.mkdir(parents=True, exist_ok=True)
COL = "L2__fibroblast__colorectum__colorectal_cancer"
mods = pd.read_csv("../../figures/previews/fib_coexpr_top20_stored.tsv", sep="\t")
# ECM program = the 9 annotated ECM/quiescent genes (the set highlighted in the Fig 6 scatter
# preview and marked in the first heatmap), all of which sit in the heatmap's ECM block.
ECM9 = ["MGP", "DCN", "OGN", "C3", "CCDC80", "MFAP5", "SFRP2", "DPT", "LUM"]
assert set(ECM9) <= set(mods.gene[mods.module == "ECM / quiescent fibroblast"])
MODULES = {"ecm": ECM9,
           "adamdec1": mods.gene[mods.module == "ADAMDEC1+ mucosal fibroblast"].tolist()}
norm = lambda s: s.astype(str).str.lower().str.replace("_", " ")
a = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r")
obs = a.obs[["cell_type", "tissue", "disease", "donor_id"]].copy(); obs.index = obs.index.astype(str)
m = (norm(obs.cell_type) == "fibroblast") & (norm(obs.tissue) == "colorectum") & norm(obs.disease).isin(["normal", "colorectal cancer"])
sub = obs[m]; dc = sub.groupby("donor_id", observed=True).disease.nunique()
cells = sub.index[sub.donor_id.isin(dc[dc == 2].index)]
name = pd.Series(a.var["feature_name"].astype(str).values, index=a.var_names.astype(str))
X = a[cells, :].X; X = X.toarray() if sp.issparse(X) else np.asarray(X); a.file.close()
X = np.log1p(X / np.maximum(X.sum(1, keepdims=True), 1) * 1e4)
var = {g: i for i, g in enumerate(name.index)}
coords = pd.read_csv(LIB / "results" / "hierarchy" / "cell_coordinates.tsv", sep="\t", index_col=0, usecols=["cell_id", COL])
coords.index = coords.index.astype(str)
df = pd.DataFrame({"cell_id": cells, "donor_id": obs.loc[cells, "donor_id"].astype(str).values,
                   "condition": np.where(norm(obs.loc[cells, "disease"]) == "normal", "normal", "cancer"),
                   "l2_score": coords.loc[cells, COL].values})
rows = []
for k, genes in MODULES.items():
    idx = [var[name[name == g].index[0]] for g in genes]
    Z = X[:, idx]; Z = (Z - Z.mean(0)) / Z.std(0)
    df[f"{k}_score"] = Z.mean(1)
    for grp, sel in [("all", np.ones(len(df), bool)), ("normal", df.condition == "normal"), ("cancer", df.condition == "cancer")]:
        r, p = pearsonr(df.l2_score[sel], df[f"{k}_score"][sel])
        rows.append({"module": k, "n_genes": len(genes), "genes": ",".join(genes), "cells": grp, "n": int(sel.sum()), "pearson_r": r, "p": p})
    # per-donor within-condition-pooled r (paired donors), to show it is not one donor
    for d, g in df.groupby("donor_id"):
        r, p = pearsonr(g.l2_score, g[f"{k}_score"])
        rows.append({"module": k, "n_genes": len(genes), "genes": "", "cells": f"donor:{d}", "n": len(g), "pearson_r": r, "p": p})
    rows.append({"module": k, "n_genes": len(genes), "genes": "", "cells": "mean score normal-cancer",
                 "n": len(df), "pearson_r": df.loc[df.condition == "normal", f"{k}_score"].mean() - df.loc[df.condition == "cancer", f"{k}_score"].mean(), "p": np.nan})
df.to_csv(OUT / "sc_fib_module_scores.csv", index=False)
S = pd.DataFrame(rows); S.to_csv(OUT / "sc_fib_module_stats.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(S.drop(columns="genes").round(4).to_string(index=False))
