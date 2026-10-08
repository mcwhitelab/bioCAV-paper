#!/usr/bin/env python3
"""PREVIEW: clustered gene x gene co-expression heatmap of the top 20 + bottom 20
fibroblast/colorectum genes by paired-donor CAV r, per arm. Same statistic as the
fibroblast-subpop supplement panel B: Pearson r of log1p(CP10K) across the
paired-donor fibroblast cells. Average linkage on 1 - r.
Writes figures/previews/fib_coexpr_top20_<arm>.png and fib_coexpr_top20_<arm>.tsv"""
import sys
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
PAIR = "fibroblast__colorectum__normal_vs_colorectal_cancer"
OUT = Path("../../figures/previews")
ECM = {"MGP", "DCN", "OGN", "C3", "CCDC80", "MFAP5", "SFRP2", "DPT", "LUM"}
# arms: CAV arm names (rank by CAV r), or "de_stat" / "de_lfc" (rank by pseudobulk DESeq2
# min10: signed Wald statistic log2fc_mle/lfcSE, or raw log2fc_mle)
ARMS = sys.argv[1:] or ["np_conv4", "stored"]
DE = LIB / "results" / "de_pseudobulk_paired_donors" / "min10" / f"{PAIR}.tsv"
norm = lambda s: s.astype(str).str.lower().str.replace("_", " ").str.strip()

a = ad.read_h5ad(LIB / "data" / "cells.h5ad", backed="r")
obs = a.obs[["cell_type", "tissue", "disease", "donor_id"]].copy(); obs.index = obs.index.astype(str)
m = (norm(obs.cell_type) == "fibroblast") & (norm(obs.tissue) == "colorectum") & norm(obs.disease).isin(["normal", "colorectal cancer"])
sub = obs[m]; dc = sub.groupby("donor_id", observed=True).disease.nunique()
cells = sub.index[sub.donor_id.isin(dc[dc == 2].index)]
X = a[cells, :].X; X = X.toarray() if sp.issparse(X) else np.asarray(X)
X = np.log1p(X / np.maximum(X.sum(1, keepdims=True), 1) * 1e4)
var = pd.Index(a.var_names.astype(str)); a.file.close()
print(f"{len(cells)} paired-donor fibroblast cells", file=sys.stderr)

# co-expression: purple (anti) - white - green (co), kept apart from the blue/orange direction bar
cmap = plt.get_cmap("PRGn")
# Module bar: each hierarchical cluster (fcluster) is named by the marker set it hits most
# (>= 2 markers); clusters that hit none stay unlabelled. Colours avoid blue/orange/green/purple.
MODULES = [  # name, markers, colour
    ("ECM / quiescent fibroblast", {"DCN", "LUM", "OGN", "MGP", "DPT", "CCDC80", "SFRP2", "MFAP5", "C3", "FBN1", "PRELP", "IGFBP6", "EFEMP1", "C1R", "GSN"}, "#8C564B"),
    ("ADAMDEC1+ mucosal fibroblast", {"ADAMDEC1", "CXCL14", "CCL13", "BMP5", "ADAM28", "NSG1"}, "#E6C229"),
    ("SFRP1+ fibroblast", {"SFRP1", "FBLN1", "TNXB", "CCN5", "TMEM100", "SMIM14"}, "#D9B38C"),
    ("CAF-like (COMP/MMP11)", {"COMP", "MMP11", "CST1"}, "#7F0000"),
    ("T cell", {"CD3D", "CD3E", "CD2", "CD27", "TRBC2", "CD69", "LTB", "CYTIP", "ACAP1", "STK17A"}, "#D62728"),
    ("B / plasma", {"IGKC", "IGHA1", "IGHA2", "JCHAIN", "MZB1", "CD79A", "MS4A1", "IGLC2", "IGLC3"}, "#E377C2"),
    ("NK / cytotoxic", {"GZMA", "CTSW", "KLRB1", "PSTPIP1"}, "#FF9896"),
    ("myeloid", {"FCER1G", "TYROBP", "AIF1", "CYBB", "C1QA", "C1QB", "MSR1", "FPR3", "S100A9", "SPP1"}, "#252525"),
    ("epithelial", {"SFN", "SPINT2", "DSG2", "LGALS2", "MAL2", "PDZK1IP1"}, "#8C8C8C"),
    ("endothelial", {"ICAM2", "CAV1", "CLEC14A", "RAMP3"}, "#17BECF"),
    ("generic non-lymphocyte", {"S100A6", "S100A10", "S100A11", "ANXA2", "LGALS3", "CD63", "LAPTM4A", "IFITM3"}, "#D0D0D0"),
]
MOD_COLOR = {n: c for n, _, c in MODULES}
def name_cluster(genes):
    hits = {n: len(set(genes) & m) for n, m, _ in MODULES}
    if hits["T cell"] >= 2 and len(set(genes) & {"CD79A", "MS4A1"}) >= 1 and hits["T cell"] >= hits["B / plasma"]:
        return "lymphocyte (T + B)"
    best = max(hits, key=hits.get)
    return best if hits[best] >= 2 else None
MOD_COLOR["lymphocyte (T + B)"] = "#D62728"
for arm in ARMS:
    if arm.startswith("de_"):
        g = pd.read_csv(DE, sep="\t").dropna(subset=["log2fc_mle", "lfcSE"])
        g["r"] = g.log2fc_mle / g.lfcSE if arm == "de_stat" else g.log2fc_mle
        key = "Wald stat" if arm == "de_stat" else "log2FC"
        cavr = pd.read_csv(f"gene_corr_stored/{PAIR}.tsv", sep="\t").set_index("gene").r
        sig = (g.padj < 0.05) & (g.log2fc_mle.abs() >= 1)
    else:
        g = pd.read_csv(f"gene_corr_{arm}/{PAIR}.tsv", sep="\t"); key = "CAV r"
    sel = pd.concat([g.nlargest(20, "r"), g.nsmallest(20, "r")])
    E = pd.DataFrame(X[:, var.get_indexer(sel.gene)], columns=sel.gene_name.values)
    C = E.corr()
    Z = linkage(squareform(1 - C.values, checks=False), "average")
    cl = fcluster(Z, t=0.85, criterion="distance")   # 1-r < 0.85 ~ r > 0.15 within a module
    if arm.startswith("de_"):
        lab = [f"{n} ({r:+.1f}; CAV {cavr.get(gid, np.nan):+.2f}{'; DE' if s_ else ''})"
               for n, r, gid, s_ in zip(sel.gene_name, sel.r, sel.gene, sig.loc[sel.index])]
    else:
        lab = [f"{n} ({r:+.2f})" for n, r in zip(sel.gene_name, sel.r)]
    C.index = C.columns = lab
    genes_by_cl = {k: [n for n, c in zip(sel.gene_name, cl) if c == k] for k in set(cl)}
    cl_name = {k: name_cluster(v) for k, v in genes_by_cl.items()}
    mod = [cl_name[c] for c in cl]
    side = pd.DataFrame({key: ["#0072B2" if r > 0 else "#D55E00" for r in sel.r],
                         "module": [MOD_COLOR[m] if m else "white" for m in mod]}, index=lab)
    cg = sns.clustermap(C, row_linkage=Z, col_linkage=Z, cmap=cmap, vmin=-1, vmax=1, center=0,
                        row_colors=side, col_colors=side, figsize=(11, 11), linewidths=0.2,
                        dendrogram_ratio=0.1, cbar_pos=(0.02, 0.83, 0.02, 0.12),
                        xticklabels=True, yticklabels=True, cbar_kws={"label": "co-expression r"})
    cg.ax_heatmap.tick_params(labelsize=7)
    cg.fig.suptitle(f"Fibroblast — colorectum: top/bottom 20 genes by {key} ({arm}); "
                    f"Pearson co-expression, {len(cells)} paired-donor cells", y=1.0, fontsize=10)
    present = [m for m in dict.fromkeys(mod) if m]
    cg.ax_heatmap.legend(handles=[Patch(color="#0072B2", label=f"{key} > 0 (up in cancer)"),
                                  Patch(color="#D55E00", label=f"{key} < 0 (down in cancer)")]
                                 + [Patch(color=MOD_COLOR[m], label=m) for m in present],
                         loc="upper left", bbox_to_anchor=(1.20, 1.14), fontsize=7, frameon=False)
    cg.savefig(OUT / f"fib_coexpr_top20_{arm}.png", dpi=200, bbox_inches="tight"); plt.close("all")
    order = [lab[i] for i in cg.dendrogram_row.reordered_ind]
    pd.DataFrame({"gene": [s.split(" (")[0] for s in order], "label": order,
                  "cluster": [cl[lab.index(s)] for s in order],
                  "module": [mod[lab.index(s)] or "" for s in order]}).to_csv(OUT / f"fib_coexpr_top20_{arm}.tsv", sep="\t", index=False)
    print(arm, {cl_name[k] or f"(unnamed {k})": v for k, v in genes_by_cl.items() if len(v) >= 2})
