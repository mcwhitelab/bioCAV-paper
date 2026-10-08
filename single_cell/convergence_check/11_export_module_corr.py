#!/usr/bin/env python3
"""PREVIEW data for a compact 'module summary' of the fibroblast top/bottom-20 heatmap
(current CAVs): long-format gene x gene Pearson co-expression (log1p CP10K, 2000 paired-donor
fibroblasts) for the named-module genes, with module and CAV r per gene.
ECM module restricted to the 9 annotated ECM genes (as in fig6D_fib_modules).
Writes figures/figure_data_staging/singlecell_modules/sc_fib_module_corr.csv, sc_fib_module_genes.csv"""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from pathlib import Path
LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
OUT = Path("../../figures/figure_data_staging/singlecell_modules")
PAIR = "fibroblast__colorectum__normal_vs_colorectal_cancer"
ECM9 = ["MGP", "DCN", "OGN", "C3", "CCDC80", "MFAP5", "SFRP2", "DPT", "LUM"]
t = pd.read_csv("../../figures/previews/fib_coexpr_top20_stored.tsv", sep="\t")
t = t[t.module.notna() & (t.module != "")]
t = t[(t.module != "ECM / quiescent fibroblast") | t.gene.isin(ECM9)]
cavr = pd.read_csv(f"gene_corr_stored/{PAIR}.tsv", sep="\t").set_index("gene_name").r
genes = pd.DataFrame({"gene": t.gene.values, "module": t.module.values, "order": range(len(t)),
                      "cav_r": cavr.reindex(t.gene).values})
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
idx = [var[name[name == g].index[0]] for g in genes.gene]
C = pd.DataFrame(X[:, idx], columns=genes.gene.values).corr()
long = C.stack().rename("r").reset_index(); long.columns = ["gene1", "gene2", "r"]
long.to_csv(OUT / "sc_fib_module_corr.csv", index=False); genes.to_csv(OUT / "sc_fib_module_genes.csv", index=False)
print(genes.groupby("module", sort=False).agg(n=("gene", "size"), mean_cav_r=("cav_r", "mean")).round(3))
