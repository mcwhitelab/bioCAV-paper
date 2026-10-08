#!/usr/bin/env python3
"""
export_fig6_fib_modules.py -- data for the lower half of Figure 6D (fibroblast/colorectum),
which replaced the continuum strips on 2026-10-05 (METHODS §6.4 items 22b-22d).

Fibroblast/colorectum is Fig 6's contrast case: its L2 axis barely separates normal from
cancer cells (AUROC 0.56 held-out) and pseudobulk DE finds ~nothing, but the genes most
correlated with the axis form coherent fibroblast programs. Panel D plots two of them against
the L2 score:
  ECM / quiescent fibroblast: the 9 annotated ECM genes of the ECM block of the clustered
    top/bottom-20 co-expression heatmap (convergence_check/08_fib_module_heatmap.py)
  ADAMDEC1+ mucosal fibroblast: that heatmap's ADAMDEC1 block (NSG1 excluded by the
    membership rule)
Cells: paired-donor fibroblast/colorectum cells (same selection as export_fig5_case_studies.py).
Score: mean over the module genes of z-scored log1p(CP10K) across these cells.
L2: --coords (default results/hierarchy/cell_coordinates.tsv; Fig 6 uses convergence_check/coords_conv1.tsv).
Caveat: module genes were chosen by |CAV r|, so the score-vs-L2 trend is partly by construction.

Writes <out_dir>/sc_fib_module_scores.csv (per cell) and sc_fib_module_stats.tsv (Pearson r of
score vs L2: all / per condition / per donor; normal - cancer mean difference).
Usage (from single_cell/): python scripts/export_fig6_fib_modules.py [out_dir=../figures/figure_data] [--coords=<tsv>]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.stats import pearsonr

LIB = Path("cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
_args = [a for a in sys.argv[1:] if not a.startswith("--coords=")]
OUT = Path(_args[0] if _args else "../figures/figure_data")
# --coords=<tsv>: hierarchy coordinates to read L2 from (default: the library file). Fig 6 uses
# convergence_check/coords_conv1.tsv (C=1 PCA-128 refit to lbfgs convergence) since 2026-10-07.
COORDS = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--coords=")),
              str(LIB / "results" / "hierarchy" / "cell_coordinates.tsv"))
COL = "L2__fibroblast__colorectum__colorectal_cancer"
MODULES = {
    "ecm": ["MGP", "DCN", "OGN", "C3", "CCDC80", "MFAP5", "SFRP2", "DPT", "LUM"],
    "adamdec1": ["BMP5", "ADAM28", "CXCL14", "CCL13", "ADAMDEC1", "STMN2", "HAPLN1"],  # HAPLN1 enters with the converged C=1 fit
}


def norm(s):
    return s.astype(str).str.lower().str.replace("_", " ")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r")
    obs = a.obs[["cell_type", "tissue", "disease", "donor_id"]].copy(); obs.index = obs.index.astype(str)
    m = ((norm(obs.cell_type) == "fibroblast") & (norm(obs.tissue) == "colorectum") &
         norm(obs.disease).isin(["normal", "colorectal cancer"]))
    sub = obs[m]; dc = sub.groupby("donor_id", observed=True).disease.nunique()
    cells = sub.index[sub.donor_id.isin(dc[dc == 2].index)]
    name = pd.Series(a.var["feature_name"].astype(str).values, index=a.var_names.astype(str))
    X = a[cells, :].X; X = X.toarray() if sp.issparse(X) else np.asarray(X); a.file.close()
    X = np.log1p(X / np.maximum(X.sum(1, keepdims=True), 1) * 1e4)
    var = {g: i for i, g in enumerate(name.index)}
    coords = pd.read_csv(COORDS, sep="\t",
                         index_col=0, usecols=["cell_id", COL])
    coords.index = coords.index.astype(str)
    df = pd.DataFrame({"cell_id": cells, "donor_id": obs.loc[cells, "donor_id"].astype(str).values,
                       "condition": np.where(norm(obs.loc[cells, "disease"]) == "normal", "normal", "cancer"),
                       "l2_score": coords.loc[cells, COL].values})
    rows = []
    for k, genes in MODULES.items():
        idx = [var[name[name == g].index[0]] for g in genes]
        Z = X[:, idx]; Z = (Z - Z.mean(0)) / Z.std(0)
        df[f"{k}_score"] = Z.mean(1)
        groups = [("all", np.ones(len(df), bool)), ("normal", (df.condition == "normal").values),
                  ("cancer", (df.condition == "cancer").values)]
        groups += [(f"donor:{d}", (df.donor_id == d).values) for d in sorted(df.donor_id.unique())]
        for grp, sel in groups:
            r, p = pearsonr(df.l2_score[sel], df[f"{k}_score"][sel])
            rows.append({"module": k, "n_genes": len(genes), "genes": ",".join(genes) if grp == "all" else "",
                         "cells": grp, "n": int(sel.sum()), "pearson_r": r, "p": p})
        rows.append({"module": k, "n_genes": len(genes), "genes": "", "cells": "mean score normal-cancer", "n": len(df),
                     "pearson_r": df.loc[df.condition == "normal", f"{k}_score"].mean()
                                  - df.loc[df.condition == "cancer", f"{k}_score"].mean(), "p": np.nan})
    df.to_csv(OUT / "sc_fib_module_scores.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "sc_fib_module_stats.tsv", sep="\t", index=False)
    print(f"Wrote sc_fib_module_scores.csv ({len(df)} cells), sc_fib_module_stats.tsv to {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
