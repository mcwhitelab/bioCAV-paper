# Figure 6 — differential expression methods (notes for the Methods section)

> **SUPERSEDED 2026-09-29.** Figure 6's DE axis is now pseudobulk DESeq2 with a paired
> design (`~ donor + condition`, samples < 10 cells dropped), not the mixed model described
> here. The current description is `METHODS_FIGURE_STEPS.md` §6.3. This file is kept for
> provenance. Note that its statement that pseudobulk DESeq2 "failed" on neutrophil/breast
> came from a per-donor aggregation bug in `de_methods.py::run_deseq2_pair`, not a power
> limit (§6.3 item 18).

Scope: Figure 6 panels B/C/D (`figures/figures.R`, lines ~1075–1445). Panel A is an
externally supplied schematic. Each of B/C/D is a DE-vs-CAV scatter over a
transcriptional-continuum strip plot, for one cell-type × tissue × normal-vs-cancer
pair:

| Panel | Pair | Cells (normal / cancer) | Donors |
|---|---|---|---|
| B | neutrophil, breast, normal vs breast cancer | 66 / 231 | 3 (3 paired) |
| C | epithelial cell, lung, normal vs lung cancer (10x 3' v2 only) | 160 / 498 | 4 (paired only) |
| D | fibroblast, colorectum, normal vs colorectal cancer | 1000 / 1000 | 7 (7 paired) |

All counts are after restriction to donors contributing both conditions. Panel C
was re-fit on the paired-donor subset (2026-09-13); B and D were unaffected
(all donors already paired) but were re-run through the same code path anyway.

## 1. The DE method actually used

**Per-gene linear mixed model with a donor random intercept**, not pseudobulk —
despite the implementing file being named `de_methods.py`. The figure reads
`results/de_mixedlm/*.tsv`, written by `run_mixedlm_pair()`
(`bioCAV/specific_scripts/de_methods.py:292`).

- Model: `expr ~ condition`, `groups = donor_id`, fitted with
  `statsmodels.formula.api.mixedlm`, REML, L-BFGS optimizer.
- Response: counts library-size-normalised to 10,000 per cell (CP10K, *not* CPM —
  the docstring says CPM but `target_sum = 1e4`), then `log2(x + 1)`. Library size
  is computed over all genes before any gene filtering.
- Condition coded 0 = normal (baseline), 1 = cancer.
- Reported effect size: the `condition` fixed-effect coefficient, written to the
  `log2fc` column. It is a difference in mean log2-normalised expression, not a
  log ratio of means. (The module's Wilcoxon branch defines `log2fc` as a true log
  ratio — same column name, different definition. Only the mixed-model definition
  applies to Figure 6.)
- Gene pre-filter: genes detected (count > 0) in at least `max(1, 1% of the pair's
  cells)`. Typically ~11–14k of 61,497 genes survive per pair.
- Zero-variance genes are assigned `log2fc = 0`, `p = 1`; fits that raise are dropped.
- Pairs with fewer than 3 donors in either condition are skipped entirely.
- Multiple testing: Benjamini–Hochberg FDR within each pair, across all genes with
  a finite p-value. Two distinct fates, easy to conflate: genes whose fit **raised**
  get `pval = NaN` and are **dropped** by the `dropna(subset=["pval"])` at the end of
  `run_mixedlm_pair` (22 / 4 / 16 genes for B / C / D — the gap between the detection
  filter's "genes kept" and the written row count); genes that return a p-value but
  fail the optimiser's convergence criterion are **retained**, flagged in the
  `converged` column, and counted in the BH denominator (207 / 33 / 127 genes,
  1.9% / 0.3% / 1.0%). State it this way — a reviewer reading the code will find the
  `dropna`.
- Genes tested: 10,681 (neutrophil/breast), 11,140 (epithelial/lung10x),
  13,154 (fibroblast/colorectum). Non-converged fits retained: 207 (1.9%),
  33 (0.3%), 127 (1.0%) respectively.

### Why a mixed model rather than pseudobulk or per-cell Wilcoxon

Both alternatives exist in the same script and were tried:

- Pseudobulk DESeq2 (`--method deseq2`, via pydeseq2, one aggregate per donor;
  results in `results/de_vs_cav/`) **failed for neutrophil/breast** — too few donors
  (3) to fit.
- Per-cell Wilcoxon (`--method wilcoxon`, the module default) treats every cell as
  independent, which is pseudoreplication across donors and yields anti-conservative
  p-values.

The mixed model is the middle ground: it uses every cell (power) while modelling
donor as the unit of non-independence (calibration), and it is the one method that
runs on all three pairs, so the three panels stay comparable.

### Execution

Pairs were enumerated from `library_structure.json` and run as a 29-task SLURM array
(`single_cell/archive/mixedlm_scaleup/{COMMANDS.sh, array_mixedlm.sbatch}`), one
`de_methods.py --method mixedlm` invocation per pair.

The lung10x pair is a special case: "lung10x" is not a real `tissue` value but a
naming convention for the assay-restricted axes introduced to fix the lung
assay confound. `single_cell/scripts/lung10x_de_and_corr.py` selects those cells
explicitly (`cell_type == epithelial cell`, `tissue == lung`,
`assay == "10x 3' v2"`) and then calls the **same** `run_mixedlm_pair()` function
with the same defaults. Statistically identical; only cell selection differs.

## 2. The comparison axis (CAV side), for contrast

The scatter's y axis is the Pearson correlation, across single cells, between that
cell's score on the pair's L2 condition CAV axis and log-CP10K expression of the
gene (`cav_gene_correlation.correlate_scores_with_genes`, min 10 cells expressing).
These were recomputed **restricted to donors contributing both conditions**
(`gene_correlation_paired_donors/`, via `scripts/paired_donor_gene_corr.py`) because
the pooled, donor-blind correlation is vulnerable to a Simpson's-paradox confound
when single-condition donors are mixed in.

## 3. Caveats a reviewer will ask about — state these explicitly

1. ~~DE and CAV axes on different cells for panel C.~~ **Resolved 2026-09-13.** The
   mixed model is now re-fit on the paired-donor cells for every panel
   (`scripts/paired_donor_de_mixedlm.py` → `results/de_mixedlm_paired_donors/`),
   so both axes come from the same cells in B, C and D. Panel C went from 1,407
   cells / 8 donors to 658 cells / 4 donors. Verification: the re-run reproduces the
   old atlas-wide tables for B and D to optimiser noise (max |Δlog2FC| 2.3e-11 and
   2.5e-10, identical gene sets, zero convergence-flag flips), confirming the cell
   selection matches the original array job.
2. **No significance threshold enters the figure.** Point categories and gene labels
   use `|log2FC| ≥ 1.5` only; `pval`/`padj` are never read by `figures.R`.
3. ~~Spearman r computed on a selected subset.~~ **Fixed 2026-09-13.** The annotated
   r is now computed over every gene with both a DE estimate and a CAV correlation,
   with n printed in the subtitle: **B r = 0.60 (n = 10,681), C r = 0.73
   (n = 11,140), D r = 0.27 (n = 13,154)** — all with p far below any threshold
   (largest p = 6e-217, for D). The previous subset-based values (B 0.81 n = 25,
   C 0.59 n = 35, D 0.55 n = 10 p = 0.10) were the liability: panel D's rested on
   10 genes and was not significant. The same fix was applied to the orthogonalized-
   projection supplemental panel, which had the identical construction.
4. **Plotted gene universe** is the inner join of DE-tested and CAV-correlated genes,
   which is now also the n behind each reported r; genes below the log2FC threshold
   are drawn as a faint grey background cloud.
   Fibroblast/colorectum has **no** gene reaching |log2FC| ≥ 1.5 (max 1.27) — that
   is a real property of the data (population-average DE has little signal in this
   cell type) and belongs in the Results text, not hedged in the Methods.
5. **Dropped pairs.** Skin/melanoma pairs were excluded as unfixably assay-confounded
   (normal = 10x 3' v2, melanoma = Smart-seq2, no overlap); the original lung pairs
   carried a milder version (normal ~91% 10x, cancer ~50% 10x) and were replaced by
   the 10x-only `lung10x` retrain.
6. **Cell caps.** The atlas h5ad is capped at 1,000 cells per condition per group at
   build time — that is why fibroblast/colorectum reads exactly 1000/1000. The
   continuum strip additionally subsamples to 1,000 cells per pair for plotting only;
   DE and CAV statistics use the full available cells.

## 4. Provenance

- DE implementation: `bioCAV/specific_scripts/de_methods.py` (`run_mixedlm_pair`)
- DE execution (Figure 6): `bioCAV-paper/single_cell/scripts/paired_donor_de_mixedlm.py`
  → `results/de_mixedlm_paired_donors/`
- DE execution (atlas-wide, superseded for Fig 6): `single_cell/archive/mixedlm_scaleup/`
  (29-task array), `scripts/lung10x_de_and_corr.py` (lung10x, unpaired)
- CAV correlations: `bioCAV-paper/single_cell/scripts/paired_donor_gene_corr.py`
- Figure data export: `bioCAV-paper/single_cell/scripts/export_fig5_case_studies.py`
  (run by `08_export_fig5_data.sh`) → `figures/figure_data/sc_*.csv`
- Rendering: `bioCAV-paper/figures/figures.R`
- Not used by Figure 6: `results/de_vs_cav/` (DESeq2), the Wilcoxon branch, and the
  superseded exploratory drivers -- all now in `single_cell/archive/` (see its README).

## 5. Draft prose

Results and Methods paragraphs are in `FIG6_DE_TEXT.md`.
