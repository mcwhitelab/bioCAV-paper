# Figure 6 — differential expression: draft Results and Methods text

> **SUPERSEDED 2026-09-29.** Figure 6's DE axis is now pseudobulk DESeq2 with a paired
> design (`~ donor + condition`, samples < 10 cells dropped), not the mixed model described
> here. The current description is `METHODS_FIGURE_STEPS.md` §6.3. This file is kept for
> provenance. Note that its statement that pseudobulk DESeq2 "failed" on neutrophil/breast
> came from a per-donor aggregation bug in `de_methods.py::run_deseq2_pair`, not a power
> limit (§6.3 item 18).

## Results (short)

To ask whether the cancer-direction CAV axis recovers the same biology a
conventional expression comparison would, we computed, for each case-study pair,
a standard normal-versus-cancer differential expression test alongside the
per-cell CAV gene correlations. Because each pair contains cells from only a
handful of donors, and donor identity is the dominant source of non-independence
in these data, we used a per-gene linear mixed model with a random intercept per
donor rather than a per-cell test (which would treat every cell as an independent
replicate) or a donor-level pseudobulk test (which is severely underpowered at
three to seven donors). Both the DE test and the CAV correlations were restricted
to donors contributing both a normal and a cancer sample, so the two axes of each
comparison are computed on exactly the same cells.

Across all genes tested, the CAV correlation and the mixed-model log2 fold change
agree for the neutrophil/breast (Spearman r = 0.60, n = 10,681 genes) and
epithelial-cell/lung (r = 0.73, n = 11,140) pairs, with the strongest concordant
genes being the expected MHC class II program (HLA-DRA, HLA-DPA1, CD74) and
epithelial keratins (KRT7, KRT18, KRT19), respectively. In the
fibroblast/colorectum pair the agreement is much weaker (r = 0.27, n = 13,154),
and no gene in that pair reaches |log2FC| >= 1.5 — population-average differential
expression simply has little to report for colorectal fibroblasts, whereas the
per-cell CAV score still resolves an ordered axis and recovers genes such as
ADAMDEC1 (CAV r = 0.21) whose fold change is essentially zero. ADAMDEC1 marks a
specific mucosal fibroblast subpopulation, consistent with a shift in
subpopulation abundance that a population-average test averages away but a
per-cell projection detects.

## Methods (long)

**Differential expression.** For each case-study pair (cell type x tissue,
normal versus a matched cancer condition) we tested differential expression with
a per-gene linear mixed model with a random intercept for donor, implemented with
`statsmodels.formula.api.mixedlm` (model `expression ~ condition`,
`groups = donor_id`, REML, L-BFGS optimiser). Raw counts were normalised to a
library size of 10,000 counts per cell and log2(x + 1)-transformed; library size
was computed across all genes before any gene filtering. Condition was coded as
0 for normal (baseline) and 1 for the cancer condition, so the fitted fixed-effect
coefficient for condition is the difference in mean log2-normalised expression
between conditions; we report this coefficient as the log2 fold change. Genes
detected (non-zero count) in fewer than 1% of the pair's cells were excluded
before fitting, as were genes with zero variance across cells; p-values were
corrected for multiple testing within each pair by the Benjamini-Hochberg
procedure over all genes returning a finite p-value. Genes whose fit raised an
error or returned no p-value were dropped (22, 4 and 16 genes for
neutrophil/breast, epithelial-cell/lung and fibroblast/colorectum respectively);
genes whose fit returned a p-value but did not satisfy the optimiser's
convergence criterion were retained and are flagged in the output tables (207,
33 and 127 genes; 1.9%, 0.3% and 1.0%). This left 10,681 genes tested for
neutrophil/breast, 11,140 for epithelial-cell/lung and 13,154 for
fibroblast/colorectum.

We chose a mixed model over the two more common alternatives deliberately, and
implemented all three in the same code path so they could be compared. A per-cell
Wilcoxon rank-sum test treats each cell as an independent observation; because
cells from one donor share donor-level biological and technical variation, this
is pseudoreplication and yields anti-conservative p-values. Pseudobulk DESeq2
(pydeseq2, one aggregate per donor) is the conservative alternative, but with
3-7 donors per pair it is badly underpowered, and it failed outright for the
neutrophil/breast pair (three donors). The mixed model uses every cell while
modelling donor as the unit of non-independence, and was the only one of the
three that ran on all pairs, keeping the panels directly comparable.

**Donor pairing.** Each pair was restricted to donors contributing cells under
both conditions before any statistic was computed. This removes, at the design
level, the Simpson's-paradox confound that arises when single-condition donors
(whose score variation reflects donor identity rather than disease state) are
pooled with dual-condition donors. Three of three donors are paired for
neutrophil/breast (297 cells: 66 normal, 231 cancer), seven of seven for
fibroblast/colorectum (2,000 cells: 1,000 normal, 1,000 cancer), and four of
eight for epithelial-cell/lung (658 cells: 160 normal, 498 cancer). The same
paired-donor cell set was used for the differential expression test and for the
CAV gene correlations in every panel.

**Assay confounding.** Disease status is confounded with sequencing assay for
some tissue/condition combinations in this atlas. All skin-epidermis/melanoma
pairs are completely confounded (normal cells assayed exclusively with 10x 3' v2,
melanoma cells exclusively with Smart-seq2) and were excluded from the analysis.
The lung pairs carried a milder version of the same confound (~91% of normal
cells versus ~50% of cancer cells assayed with 10x 3' v2); for these, the CAV
axis was retrained from scratch on cells assayed with 10x 3' v2 only, and both
the differential expression test and the gene correlations were restricted to the
same assay-matched cells. Pairs labelled "lung10x" throughout refer to this
assay-restricted analysis.

**Comparison with CAV gene correlations.** The CAV axis value plotted against
log2 fold change is the Pearson correlation, across the pair's cells, between a
cell's score on the level-2 cancer-direction CAV axis and that cell's
log-normalised expression of the gene (genes expressed in fewer than 10 cells
excluded). The Spearman correlation annotated on each panel is computed over
every gene with both a differential expression estimate and a CAV correlation
(the intersection of the two gene sets; n is given in each panel), not over any
selected subset. The |log2FC| >= 1.5 threshold shown in the figure affects only
point styling and which genes are labelled; it does not enter the reported
correlation or any significance calculation. No significance threshold was
applied to the genes shown.

**Implementation.** Differential expression: `run_mixedlm_pair()` in
`bioCAV/specific_scripts/de_methods.py`, driven for the case-study pairs by
`bioCAV-paper/single_cell/scripts/paired_donor_de_mixedlm.py`. CAV gene
correlations: `bioCAV-paper/single_cell/scripts/paired_donor_gene_corr.py`.
Figure data export: `scripts/export_fig5_case_studies.py`; rendering:
`bioCAV-paper/figures/figures.R`.
