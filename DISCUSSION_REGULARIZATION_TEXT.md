# Discussion draft — choice of regularisation strength (C)

Drafted 2026-10-05; MYC paragraph added after Fig 5B/C and the PCA supplementary were promoted to C = 10⁻⁴ PCA-128. Numbers from METHODS_FIGURE_STEPS.md §1.6, §6.1 items 5/5b, §6.4, and
`single_cell/convergence_check/COMPARISON.md`, `motif_clustering/results/converged/`.

---

**Regularisation strength determines which part of a concept a CAV represents.** A CAV is
the weight vector of an L2-regularised logistic regression. Its regularisation strength is
therefore not a nuisance parameter: it decides what the vector encodes. Under strong
regularisation (C·n ≈ 1, i.e. C = 10⁻⁴ for our training sets of ~10⁴ embeddings), the
penalty is as large as the loss. The solution is then pulled toward the difference between
the positive and background class means, so the CAV captures the dominant contrast between
concept and background. Under weak regularisation (C = 1, C·n ≈ 10⁴), the fit can use the
covariance of the embeddings. It down-weights directions that vary widely across all inputs
and recruits lower-variance directions that also separate the classes, so several distinct
patterns of variation can contribute to a single vector. Neither setting is correct in
general. The right choice depends on whether the downstream question concerns the shared
signature of a concept or the structure within it.

**For protein annotation, the dominant contrast is the concept.** A Pfam family, GO term or
EC class is defined by what its members share, and each CAV is used directly as concept
versus background. Trained to convergence at C = 10⁻⁴, Pfam CAVs clustered more coherently
by clan than our original library (ARI 0.185 vs 0.163; median clan span 39 vs 258 clusters).
The same CAVs fit to convergence at C = 1 fragmented clan structure (ARI 0.078), consistent
with weakly regularised fits latching onto family-idiosyncratic directions. Moving the
GO/EC library to C = 10⁻⁴ raised top-ranked EC third-level predictions from 31% to 55% of
held-out enzymes and left GO rankings essentially unchanged. The original library had
obtained its useful structure implicitly, through the default stopping tolerance of the
optimiser, which halted fits well short of their optimum. Fitting every CAV to convergence
at an explicit C makes that regularisation a stated, reproducible choice rather than a side
effect of the solver.

**For single-cell states, the question is about structure within a concept.** Each
single-cell condition axis is the residual between two closely related CAVs, for example
cancer and normal colorectal fibroblasts, each trained against the same atlas-wide
background. At C = 10⁻⁴ both CAVs collapse onto the shared "fibroblast versus everything
else" direction (cosine between them 0.79 → 0.98). Their residual then reduces to the raw
mean shift between cancer and normal cells, which in this atlas largely reflects immune
contamination of the fibroblast annotation: correlation of the axis with PTPRC rose from
0.12 to 0.42, and with CD52 from 0.15 to 0.49. At C = 1 the same axis instead resolved
distinct co-expressed programs. These included an extracellular-matrix/quiescent-fibroblast
module, which tracked axis position in all seven donors, and an ADAMDEC1-centred mucosal
fibroblast module, which was clear in three of seven. Because these genes were selected by
their correlation with the axis, the modules show what the axis can resolve rather than an
independent validation. Strong regularisation did increase agreement with conventional
differential expression (fibroblast ρ 0.23 → 0.39), but this is expected by construction
for an axis that approaches the mean difference. It is not evidence of a better axis. Weak
regularisation is what lets multiple, partly independent biological programs contribute to
a cell-state direction.

**Dimensionality reduction is a second resolution setting.** The MYC interface CAV makes
this concrete. It has a single positive, the MYC bHLH domain. Under strong regularisation it
reduces to the standardised MYC window, so a proteome scan becomes a z-scored similarity
search. That search returns bHLH proteins almost exclusively (97 of the top 100), but it does
not distinguish proteins that bind MYC's CACGTG E-box from other bHLH factors: the AUROC for
CACGTG binders versus other bHLH proteins is 0.53. Projecting onto the leading principal
components of the background before fitting restricts the contrast to the embedding's major
axes of variation. Between 96 and 160 components, this separated CACGTG binders from other bHLH
proteins (AUROC 0.82–0.86, 0.85 at 128 components), and 16 of the top 20 hits were annotated
CACGTG binders. A second single-positive concept, the MAD1 bHLH segment, showed the same
plateau (AUROC 0.80–0.87). Too few components lost the family altogether: at 16 or fewer,
fewer than 10 of the top 100 hits were bHLH proteins. Rather than a single tuned value,
dimensionality therefore defines a range over which the CAV resolves a functional subclass
within a fold family.

**Practical guidance.** Because scikit-learn sums the logistic loss over samples, a fixed C
corresponds to different effective regularisation at different training-set sizes. A
cross-species SERK/CIK interface CAV trained on 32 sequences retained no discrimination at
C = 10⁻⁴ (AUROC 1.00 → 0.54). Our software therefore defaults to C = 1/n for user-defined
concepts, which reproduces the library setting for ~10⁴ samples. We recommend
strong regularisation (C·n ≈ 1) when a CAV is used to annotate or retrieve members of a
concept. We recommend weaker regularisation when CAVs are subtracted from one another or
used to resolve heterogeneity within a concept, together with dimensionality reduction to
keep the fit well conditioned. At C = 1, projecting onto 128 principal components was needed
for the optimiser to converge on the rank-deficient Geneformer embedding. At C = 10⁻⁴, PCA
had little effect (cosine 0.89 with and without it), because the penalty alone conditions
the problem. Both C and the number of retained components are best understood as resolution
settings, from the coarsest shared contrast to finer within-concept structure, and should
be reported alongside any CAV.
