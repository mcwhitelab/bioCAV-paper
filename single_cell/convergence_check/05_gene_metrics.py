#!/usr/bin/env python3
"""Gene-level comparison of arms for the 3 Fig 6 pairs: top-5/bottom-5 overlap
with the stored library, Spearman(CAV r, pseudobulk log2fc_mle) on the DE gene set,
Spearman(CAV r, label r), Spearman(CAV r, stored CAV r), and the named callouts."""
import pandas as pd, numpy as np
from pathlib import Path
from scipy.stats import spearmanr
R = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d/results")
PAIRS = ["neutrophil__breast__normal_vs_breast_cancer", "epithelial_cell__lung10x__normal_vs_lung_cancer",
         "fibroblast__colorectum__normal_vs_colorectal_cancer"]
ARMS = ["stored", "orig", "conv1", "conv4", "np_orig", "np_conv1", "np_conv4"]
CALL = ["ADAMDEC1", "MGP", "DCN", "OGN", "C3", "CCDC80", "MFAP5", "SFRP2", "DPT", "LUM"]
rows, callrows, tops = [], [], []
# genes the live run (gene_correlation_paired_donors) had; the stored rerun must match it
for p in PAIRS:
    live = pd.read_csv(R / "gene_correlation_paired_donors" / f"{p}.tsv", sep="\t").set_index("gene")
    de = pd.read_csv(R / "de_pseudobulk_paired_donors/min10" / f"{p}.tsv", sep="\t").set_index("gene")
    lab = pd.read_csv(R / "label_correlation_paired_donors" / f"{p}.tsv", sep="\t").set_index("gene")
    g = {a: pd.read_csv(f"gene_corr_{a}/{p}.tsv", sep="\t").set_index("gene") for a in ARMS}
    st = g["stored"]
    assert np.allclose(st.r.reindex(live.index), live.r, atol=1e-4, equal_nan=True), f"stored != live for {p}"
    top = lambda d: (set(d.r.nlargest(5).index), set(d.r.nsmallest(5).index))
    t0 = top(st)
    for a in ARMS:
        d = g[a]; tt = top(d)
        i = d.index.intersection(de.index); j = d.index.intersection(lab.index)
        rows.append({"pair": p.split("__normal_vs")[0], "arm": a,
                     "rho_vs_stored_r": spearmanr(d.r, st.r.reindex(d.index), nan_policy="omit")[0],
                     "top5_overlap": len(tt[0] & t0[0]), "bot5_overlap": len(tt[1] & t0[1]),
                     "rho_DE_lfc_mle": spearmanr(d.r.loc[i], de.log2fc_mle.loc[i], nan_policy="omit")[0],
                     "rho_label_r": spearmanr(d.r.loc[j], lab.label_r.loc[j], nan_policy="omit")[0],
                     "n_DE_genes": len(i), "max_abs_r": d.r.abs().max()})
        nm = d.set_index("gene_name").r
        tops.append({"pair": p.split("__normal_vs")[0], "arm": a,
                     "top5": ",".join(d.loc[d.r.nlargest(5).index, "gene_name"]),
                     "bottom5": ",".join(d.loc[d.r.nsmallest(5).index, "gene_name"])})
        if p.startswith("fibroblast"):
            labn = lab.set_index("gene_name").label_r
            for gname in CALL:
                callrows.append({"gene": gname, "arm": a, "cav_r": nm.get(gname, np.nan), "label_r": labn.get(gname, np.nan),
                                 "rank_of_" + str(len(nm)): int((nm.rank(ascending=False)).get(gname, -1)) if gname in nm else None})
M = pd.DataFrame(rows); T = pd.DataFrame(tops); C = pd.DataFrame(callrows)
M.to_csv("gene_metrics.tsv", sep="\t", index=False); T.to_csv("top_genes.tsv", sep="\t", index=False); C.to_csv("fibroblast_callouts.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(M.round(3).to_string()); print(); print(T.to_string()); print()
print(C.pivot(index="gene", columns="arm", values="cav_r")[ARMS].join(C.groupby("gene").label_r.first()).round(3).to_string())
