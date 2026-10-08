#!/usr/bin/env python3
"""
PREVIEW (not in figures.R): does the per-term vs shared scaler choice matter
more for GO terms with more diverse CAV training positives?

Inputs: go_perterm_vs_shared_valauroc.tsv (654 terms; val-protein AUROC of the
original per-term-scaler CAV vs the retrained shared-scaler CAV, same val set)
and go_eval/diversity/term_diversity_{ont}.csv (families_per_100 of the CAV
training positives; DIAMOND >=30% id / >=50% cov components).
"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
D = pd.read_csv("go_perterm_vs_shared_valauroc.tsv", sep="\t")
div = pd.concat([pd.read_csv(f"../../go_eval/diversity/term_diversity_{o}.csv").assign(ns=o) for o in ("mf", "bp", "cc")])
M = D.merge(div[["go_term", "families_per_100", "n_train"]], left_on="term", right_on="go_term")
M["dAUC"] = M.auc_shared - M.auc_per_term
M.to_csv("go_perterm_vs_shared_with_diversity.tsv", sep="\t", index=False)
r1 = spearmanr(M.families_per_100, M.dAUC); r2 = spearmanr(M.families_per_100, M.cos)
print(f"terms with diversity: {len(M)} ({M.ns.value_counts().to_dict()})")
print(f"Spearman families_per_100 vs dAUC: rho={r1.statistic:.3f} p={r1.pvalue:.3g}; vs cosine: rho={r2.statistic:.3f} p={r2.pvalue:.3g}")
M["bin"] = pd.qcut(M.families_per_100, 4, labels=["least diverse", "Q2", "Q3", "most diverse"])
print(M.groupby("bin", observed=True).agg(n=("dAUC", "size"), fam100=("families_per_100", "median"),
      median_dAUC=("dAUC", "median"), frac_abs_gt_0_01=("dAUC", lambda x: (x.abs() > 0.01).mean()),
      median_cos=("cos", "median")).round(4).to_string())

COL = {"bp": "#0072B2", "cc": "#D55E00", "mf": "#009E73"}; NAME = {"bp": "BP", "cc": "CC", "mf": "MF"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, (a, b) = plt.subplots(1, 2, figsize=(5.5, 2.3), constrained_layout=True)
for ns, g in M.groupby("ns"):
    lab = f"{NAME[ns]} (n={len(g)})"
    a.scatter(g.families_per_100, g.dAUC, s=7, color=COL[ns], alpha=0.6, lw=0, label=lab)
    b.scatter(g.families_per_100, g.cos, s=7, color=COL[ns], alpha=0.6, lw=0, label=lab)
a.axhline(0, color="#8a8a8a", lw=0.6, ls="--")
a.set_xlabel("training-positive diversity\n(sequence families per 100 positives)"); a.set_ylabel("AUROC, shared − per-term scaler")
a.text(0.02, 0.98, f"Spearman ρ = {r1.statistic:.2f} (p = {r1.pvalue:.2g})", transform=a.transAxes, va="top", fontsize=6.5)
b.set_xlabel("training-positive diversity\n(sequence families per 100 positives)"); b.set_ylabel("cosine, per-term vs shared CAV")
b.text(0.02, 0.04, f"Spearman ρ = {r2.statistic:.2f} (p = {r2.pvalue:.2g})", transform=b.transAxes, va="bottom", fontsize=6.5)
b.legend(frameon=False, fontsize=6.5, loc="center right", markerscale=1.6)
for ax, l in ((a, "A"), (b, "B")): ax.text(-0.2, 1.04, l, transform=ax.transAxes, fontweight="bold", fontsize=8)
fig.savefig("preview_scaler_vs_diversity.png", dpi=300); fig.savefig("preview_scaler_vs_diversity.pdf")
