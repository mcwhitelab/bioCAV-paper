#!/usr/bin/env python3
"""PREVIEW: clustered gene x gene co-expression heatmap of the top 20 + bottom 20
fibroblast/colorectum genes by paired-donor CAV r, per arm. Same statistic as the
fibroblast-subpop supplement panel B: Pearson r of log1p(CP10K) across the
paired-donor fibroblast cells. Average linkage on 1 - r.
Writes figures/previews/fib_coexpr_top20_<arm>.png and fib_coexpr_top20_<arm>.tsv"""
import sys
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
import matplotlib; matplotlib.use("Agg")
# editable text in the PDF (TrueType, not Type 3) and the house Arial-metric font used by figures.R
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42,
                            "font.family": "Liberation Sans"})
import matplotlib.pyplot as plt, seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import textwrap
from scipy.spatial.distance import squareform

LIB = Path("../cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d")
PAIR = "fibroblast__colorectum__normal_vs_colorectal_cancer"
OUT = Path("../../figures/previews")
ECM = {"MGP", "DCN", "OGN", "C3", "CCDC80", "MFAP5", "SFRP2", "DPT", "LUM"}
# arms: CAV arm names (rank by CAV r), or "de_stat" / "de_lfc" (rank by pseudobulk DESeq2
# min10: signed Wald statistic log2fc_mle/lfcSE, or raw log2fc_mle)
ARMS = sys.argv[1:] or ["np_conv4", "stored"]
CAV_REF = "conv1"   # CAV r shown next to DE-ranked genes: converged C=1 PCA-128 (was stored)
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
        cavr = pd.read_csv(f"gene_corr_{CAV_REF}/{PAIR}.tsv", sep="\t").set_index("gene").r
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
    # membership strength: mean co-expression with the other members of the same cluster.
    # A gene whose strength is < half its module's median loses the module label (e.g. RPL39
    # rides along with the plasma genes at 0.15 vs median 0.45) and moves to the edge of its
    # cluster as an unlabelled gene. Relative, so uniformly loose modules are not stripped.
    REL_MEMBER_R = 0.5
    Cv = C.values
    member_r = np.array([Cv[i, [j for j in range(len(cl)) if cl[j] == cl[i] and j != i]].mean()
                         if (cl == cl[i]).sum() > 1 else np.nan for i in range(len(cl))])
    med = {c: np.nanmedian(member_r[cl == c]) for c in set(cl)}
    weak = np.array([m is not None and r < REL_MEMBER_R * med[c] for m, r, c in zip(mod, member_r, cl)])
    for i in np.where(weak)[0]:
        print(f"  {arm}: {sel.gene_name.iloc[i]} dropped from '{mod[i]}' (mean r {member_r[i]:.2f}, module median {med[cl[i]]:.2f})")
        mod[i] = None
    new_cl = cl.copy(); nxt = cl.max() + 1
    for c in set(cl[weak]):
        w = [i for i in np.where(weak & (cl == c))[0]]
        # weak genes that co-express with each other (e.g. CD79A + MS4A1 split off a T-cell
        # block) form their own sub-block and get named; isolated ones become singletons
        if len(w) >= 2 and np.mean([Cv[i, j] for i in w for j in w if i != j]) >= 0.3:
            new_cl[w] = nxt; nm = name_cluster([sel.gene_name.iloc[i] for i in w])
            for i in w: mod[i] = nm
            nxt += 1
        else:
            for i in w: new_cl[i] = nxt; nxt += 1
        rest = [i for i in np.where(cl == c)[0] if i not in w]   # re-name what is left
        nm = name_cluster([sel.gene_name.iloc[i] for i in rest])
        for i in rest: mod[i] = nm
    cl = new_cl
    # --- order: dendrogram leaf order (fcluster clusters are contiguous in it), no trees drawn ---
    leaves = dendrogram(Z, no_plot=True)["leaves"]
    # within each original cluster run, push weak (unlabelled) members to the end
    base = fcluster(Z, t=0.85, criterion="distance")
    runs, cur = [], []
    for i in leaves:
        if cur and base[cur[-1]] != base[i]:
            runs.append(cur); cur = []
        cur.append(i)
    runs.append(cur)
    leaves = [i for r in runs for i in sorted(r, key=lambda i: weak[i])]
    blocks = []                                   # [name, first, last] runs of the same cluster
    for pos, i in enumerate(leaves):
        if blocks and cl[leaves[blocks[-1][1]]] == cl[i]:
            blocks[-1][2] = pos
        else:
            blocks.append([mod[i], pos, pos])
    merged = []                                   # adjacent unnamed runs collapse into one gap-free run
    for b in blocks:
        if merged and not b[0] and not merged[-1][0]:
            merged[-1][2] = b[2]
        else:
            merged.append(b)
    blocks = merged
    M = C.values[np.ix_(leaves, leaves)]
    n = len(leaves)
    fig, ax = plt.subplots(figsize=(8.2, 7.6))
    im = ax.imshow(M, cmap=cmap, vmin=-1, vmax=1, interpolation="nearest", aspect="equal")
    for _, a0, _ in blocks[1:]:                   # thin white gaps between blocks
        ax.axhline(a0 - 0.5, color="white", lw=1.6); ax.axvline(a0 - 0.5, color="white", lw=1.6)
    ax.set_xticks([]); ax.set_yticks(range(n))
    ax.set_yticklabels([lab[i] for i in leaves], fontsize=6.5)
    ax.yaxis.tick_right()
    for t, i in zip(ax.get_yticklabels(), leaves):   # label colour = direction, Fig 6 cancer/normal colours
        t.set_color("#c0392b" if sel.r.iloc[i] > 0 else "#3a6fad")
    ax.tick_params(axis="y", length=0, pad=2)
    for s in ax.spines.values(): s.set_visible(False)
    # module names above their blocks (named blocks only), quiet grey brackets
    for name, a0, a1 in blocks:
        if not name: continue
        y = -1.2
        ax.plot([a0 - 0.35, a1 + 0.35], [y, y], color="#666666", lw=0.8, clip_on=False)
        w = a1 - a0 + 1
        txt = name if w >= 6 else textwrap.fill(name, max(6, int(w * 2.2)), break_long_words=False)
        ax.text((a0 + a1) / 2, y - 0.6, txt, ha="center", va="bottom", fontsize=6.5 if w >= 3 else 5.5, color="#333333",
                linespacing=1.0, clip_on=False)
    ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(n - 0.5, -0.5)
    cax = fig.add_axes([0.06, 0.12, 0.018, 0.18])
    cb = fig.colorbar(im, cax=cax, ticks=[-1, 0, 1]); cb.ax.tick_params(labelsize=6)
    cb.set_label("co-expression r", fontsize=6.5); cax.yaxis.set_ticks_position("left"); cax.yaxis.set_label_position("left")
    unit = "CAV r" if key == "CAV r" else f"{key}; CAV r"
    ax.set_title(f"Fibroblast — colorectum: top/bottom 20 genes by {key} ({arm})\n"
                 f"Pearson co-expression across {len(cells)} paired-donor cells; label = gene ({unit}); red = up in cancer, blue = down",
                 fontsize=7.5, pad=46, color="#222222")
    fig.savefig(OUT / f"fib_coexpr_top20_{arm}.png", dpi=250, bbox_inches="tight")
    fig.savefig(OUT / f"fib_coexpr_top20_{arm}.pdf", bbox_inches="tight"); plt.close(fig)
    order = [lab[i] for i in leaves]
    # long-format matrix in display order, for the R supplementary figure
    SUPP = Path("../../figures/figure_data_staging/singlecell_heatmaps"); SUPP.mkdir(parents=True, exist_ok=True)
    gl = [sel.gene_name.iloc[i] for i in leaves]
    pd.DataFrame([(gl[x], gl[y], M[x, y]) for x in range(n) for y in range(n)],
                 columns=["gene1", "gene2", "r"]).to_csv(SUPP / f"fib_coexpr_{arm}_corr.csv", index=False)
    pd.DataFrame({"pos": range(1, n + 1), "gene": gl, "value": [sel.r.iloc[i] for i in leaves],
                  "cav_r": [sel.r.iloc[i] if key == "CAV r" else cavr.get(sel.gene.iloc[i], np.nan) for i in leaves],
                  "module": [mod[i] or "" for i in leaves], "key": key,
                  "log2fc": [sel.log2fc_mle.iloc[i] if "log2fc_mle" in sel else np.nan for i in leaves]}).to_csv(SUPP / f"fib_coexpr_{arm}_genes.csv", index=False)
    pd.DataFrame({"gene": [s.split(" (")[0] for s in order], "label": order,
                  "cluster": [cl[lab.index(s)] for s in order],
                  "module": [mod[lab.index(s)] or "" for s in order],
                  "member_mean_r": [round(member_r[lab.index(s)], 3) for s in order]}).to_csv(OUT / f"fib_coexpr_top20_{arm}.tsv", sep="\t", index=False)
    print(arm, {cl_name[k] or f"(unnamed {k})": v for k, v in genes_by_cl.items() if len(v) >= 2})
