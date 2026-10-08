# Single-cell CAVs: convergence and C=1e-4 check (2026-10-01)

Library: `cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d` (93 CAVs, Geneformer V2-104M,
shared scaler → shared PCA-128 → LR C=1, balanced, lbfgs tol 1e-4, max_iter 1000, n = 500 pos (cap) + 10,000 ref neg).
Comparison only: nothing under `cav_library/*/results`, `figures/` or `figure_data/` was touched.

Arms (all in the stored PCA-128 frame, same pos/neg):
- `stored` — the library's concept_v1.npy
- `orig`   — refit with the stored recipe (C=1, tol 1e-4)
- `conv1`  — C=1, tol 1e-10, max_iter 1e5
- `conv4`  — C=1e-4, tol 1e-10 (library rule; C·n ≈ 1.05, ≈ the new 1/n default 9.5e-5)

## 1. Did the current C=1 CAVs fail to converge? Essentially no.
`01_refit_cavs.py` → `refit_summary.tsv`
- `orig` reproduces stored: median cos 1.0000 (min 0.998). Stored fits stop after 22–171 lbfgs iterations (median 85).
- `conv1` (needs 74–640 iterations) vs stored: **median cos 0.9995**, 5th pct 0.921. Unlike Pfam no-PCA (cos ~0.7), PCA-128 makes the
  early-stopped fit almost the optimum, as in the Pfam PCA-128 check (0.97–0.99).
- Only exceptions are small/odd sets: all skin_epidermis/melanoma CAVs (0.85–0.93, already excluded from Fig 6) and
  neutrophil/breast normal (n_pos=53, 0.930) and cancer (0.967).
- `conv4` vs stored: median cos **0.63** (0.52–0.74). C=1e-4 is a substantively different estimator, not a convergence fix.
  Train AUROC drops 0.999 → 0.994.

## 2. Fig 6 case pairs, after the hierarchy (L0 mean → L1/L2 Gram-Schmidt residuals)
`02_build_coords.py` (gate: stored CAVs reproduce `results/hierarchy/cell_coordinates.tsv` to 5e-5 on all 101 axes),
`04_axis_metrics.py` → `axis_metrics.tsv`. Cells = paired donors only.

| pair | arm | cos(L2, stored L2) | r(cell L2 score, stored) | L2 AUROC N vs C (in-sample) |
|---|---|---|---|---|
| neutrophil/breast | conv1 | 0.952 | 0.982 | 0.870 (stored 0.869) |
| neutrophil/breast | conv4 | 0.564 | 0.847 | 0.861 |
| epithelial/lung10x | conv1 | 0.996 | 0.999 | 0.841 (stored 0.841) |
| epithelial/lung10x | conv4 | 0.506 | 0.732 | 0.904 |
| fibroblast/colorectum | conv1 | 0.999 | 1.000 | 0.585 (stored 0.585) |
| fibroblast/colorectum | conv4 | 0.351 | 0.514 | 0.644 |

The residual step amplifies small raw differences: even the `orig` refit (raw cos 0.9999) moves one non-case axis
(L2 epithelial/colorectum) by up to 0.5 (sd 1.3). For the three case pairs, `orig` = stored.

cos(raw CAV, balanced class-mean difference): stored 0.38–0.50, conv1 0.38–0.57, **conv4 0.79–0.91**. C=1e-4 moves well
towards the mean difference but is not a pure mean-difference axis (unlike MYC, n=1).

## 3. Gene level (paired-donor CAV r; `05_gene_metrics.py` → `gene_metrics.tsv`, `top_genes.tsv`, `fibroblast_callouts.tsv`)

| pair | arm | ρ(r, stored r) | top5/bot5 kept | ρ(r, DE log2fc_mle) | ρ(r, label r) |
|---|---|---|---|---|---|
| neutrophil/breast | stored | 1 | 5/5 | 0.750 | 0.663 |
| | conv1 | 0.986 | 5/4 | 0.758 | 0.668 |
| | conv4 | 0.912 | 4/3 | 0.756 | 0.648 |
| epithelial/lung10x | stored | 1 | 5/5 | 0.807 | 0.712 |
| | conv1 | 0.999 | 5/5 | 0.805 | 0.709 |
| | conv4 | 0.797 | 2/4 | 0.810 | 0.762 |
| fibroblast/colorectum | stored | 1 | 5/5 | 0.226 | 0.335 |
| | conv1 | 1.000 | 4/5 | 0.223 | 0.332 |
| | conv4 | 0.640 | 1/0 | 0.379 | 0.451 |

Top/bottom 5 under conv4:
- lung10x top: IGKC, SPP1, KRT7, KRT18, PAEP (stored: KRT7, KRT18, HMGB3, KRT19, RAB25); bottom adds CD68.
- fibroblast top: CXCR4, CD52, CD79A, CORO1A, LTB (lymphocyte genes); bottom: CST3, S100A6, IFITM3, TIMP1, CD63
  (stored top: IGKC, IGHA1, ADAMDEC1, CXCR4, RPL39; bottom: MGP, CCDC80, C3, DCN, OGN).

Fibroblast callouts (CAV r; label r ≈ 0 for all):
- **ADAMDEC1: stored 0.207, conv1 0.212, conv4 −0.033.** The ADAMDEC1 callout disappears under C=1e-4.
- ECM/quiescent program (MGP, DCN, OGN, C3, CCDC80, MFAP5, SFRP2, DPT, LUM): stored −0.30 to −0.37, conv4 −0.23 to −0.34.
  It survives but drops out of the bottom 5.

## Reading
- **Convergence:** the current single-cell CAVs are effectively converged. Converging at C=1 leaves lung10x and fibroblast
  unchanged: L2 cos ≥0.996, identical ρ and AUROC to about 0.005. The fibroblast top 5 is 4/5 (STMN2 in for RPL39), which is the
  refit noise floor, since `orig` also swaps that slot (CD3D). **neutrophil/breast is the one case pair that moves slightly.** It has
  the smallest positive set (n_pos 53 normal) and 3 donors: L2 cos 0.952, ρ_DE 0.750→0.758, max |r| 0.77→0.82, bottom 5 4/5 (ETS1 in
  for RPS29), same story. No retrain is needed for convergence reasons.
- **Fragility of the residual axes:** across all 101 axes, the `orig` refit (raw cos ≈0.9999) gives axis cos ≥0.9977. The worst is
  L2 epithelial/colorectum: cell-score r 0.999, but max |Δ| 0.5 on single cells. This is harmless at the level of the reported numbers.
- **C=1e-4** is a different axis (L2 cos 0.35–0.56 to current). It raises DE/label concordance (fibroblast ρ_DE 0.23→0.38) and the
  in-sample L2 AUROC. That is expected for an estimator pulled toward the class-mean difference (cos 0.79–0.91), so it is not
  independent evidence of improvement. The concrete change to the Fig 6 story is that **ADAMDEC1 is lost** in fibroblast
  (r 0.21 → −0.03). Immune/plasma genes sit at the extremes in every arm (stored fibroblast top: IGKC, IGHA1, CXCR4; stored lung10x
  bottom: FCER1G, TYROBP). Under conv4 they take over: fibroblast top 5 is all lymphocyte genes (CXCR4, CD52, CD79A, CORO1A, LTB), and
  the lung10x top gains IGKC/SPP1 while its bottom gains CD68. The fibroblast ECM program survives (r −0.23 to −0.34) but leaves the bottom 5.
- The only argument for C=1e-4 here is consistency with the library rule (METHODS §1.6, C·n ≈ 1). `specific_scripts/run_cav_pipeline.py`
  still passes regularization_C=1.0.

## METHODS mismatches (not edited)
- §1.6 says "Every CAV in the main figures is trained in the full 1152-d standardized space, except the MYC CAV". Fig 6 single-cell
  CAVs are 768-d Geneformer embeddings, shared scaler → shared PCA-128, C=1, lbfgs default tol (stopped at 22–171 iterations).
- §6.1 gives no estimator. Item 7 points at `cav_hierarchy.py`, but the current version projects in scaler space with no PCA;
  the live coordinates came from the earlier PCA version (+ `scripts/extend_hierarchy_lung10x.py`).

## 4. No-PCA arms (added 2026-10-01, for consistency with the protein libraries)
Same pos/neg, shared scaler frame (768-d), no PCA. Verified first that all 93 CAVs share one scaler (identical mean/scale,
n_samples_seen 10,000) and the same 10k negatives (row order differs). `06_refit_nopca.py` → `refit_nopca_summary.tsv`;
arms `np_orig` (C=1, tol 1e-4), `np_conv1` (C=1, converged), `np_conv4` (C=1e-4, converged = protein-library estimator).
PCA arms are mapped into the scaler frame (`pca.components_.T @ v`, `axes_<arm>_sc.npz`) for cross-frame cosines.

Raw CAVs (median over 93):
- no PCA at C=1 is moderately early-stopped: np_orig stops at 18–52 iterations and np_conv1 needs 61–446. cos(np_orig, np_conv1) has median 0.94 (0.71–0.97), milder than Pfam's ~0.7. CORRECTED 2026-10-02: this line earlier said "badly, like Pfam" on the basis of the cos-to-stored value, which measures an estimator difference.
  np_orig vs stored cos 0.31, np_conv1 vs stored 0.45.
- **np_conv4 vs PCA conv4: cos 0.89** (0.83–0.93). np_conv4 vs mean difference 0.86. At C=1e-4 the penalty dominates, so PCA
  barely matters. np_conv4 vs stored: 0.56.

Fig 6 L2 axes (`axis_metrics.tsv`, scaler frame):

| pair | np_conv4 vs PCA conv4 | np_conv4 vs stored | AUROC stored / conv4 / np_conv4 |
|---|---|---|---|
| neutrophil/breast | 0.84 | 0.46 | 0.869 / 0.861 / 0.889 |
| epithelial/lung10x | 0.84 | 0.43 | 0.841 / 0.904 / 0.913 |
| fibroblast/colorectum | 0.73 | 0.25 | 0.585 / 0.644 / 0.677 |

Gene level (`gene_metrics.tsv`, `top_genes.tsv`). np_conv4 reproduces the PCA conv4 picture almost gene for gene:
- neutrophil ρ_DE 0.759 (PCA conv4 0.756); lung10x 0.814 (0.810); fibroblast 0.391 (0.379).
- Top/bottom 5 are near-identical to PCA conv4. Fibroblast top: CXCR4, CD52, CORO1A, LTB, CD79A; bottom: CST3, IFITM3, S100A6, TIMP1, CD63.
  lung10x top: IGKC, SPP1, KRT7, KRT18, PAEP (identical).
- ADAMDEC1 r −0.029 (PCA conv4 −0.033, stored 0.207). ECM program r −0.22 to −0.34 (same as PCA conv4).
- The no-PCA C=1 arms are a third, different picture: fibroblast top H3-3A, IGKC, IGHA1; bottom HBB, TIMP1, MGP; ADAMDEC1 0.07.

**Reading:** the consistent estimator (shared scaler, no PCA, C=1e-4, converged) gives essentially the same Fig 6 as
PCA-128 + C=1e-4. The PCA choice stops mattering once C is set. The current Fig 6 gene story, including ADAMDEC1, is specific to
PCA-128 + C=1. The consistent estimator gains DE/label concordance (expected: it is near the mean difference) and loses ADAMDEC1.

**Immune-marker check (paired-donor CAV r, stored / PCA conv4 / np_conv4):**
- fibroblast/colorectum: PTPRC 0.12 / 0.44 / 0.42, CD52 0.15 / 0.50 / 0.49, MS4A1 0.13 / 0.41 / 0.39, CD79A 0.15 / 0.47 / 0.44.
  COL1A1 −0.19 / −0.32 / −0.31, PDGFRA −0.13 / −0.21 / −0.20.
  The C=1e-4 fibroblast L2 axis largely tracks immune-marker-expressing cells (contamination or doublets in the fibroblast
  annotation), and it points away from fibroblast identity. The stored axis has the same tilt, only weaker.
- epithelial/lung10x: LYZ −0.33 → −0.52, CD68 −0.31 → −0.52 (myeloid signal on the normal side). KRT7 is unchanged at 0.29.
- The ECM program (MGP, DCN, OGN …) keeps CAV r −0.22 to −0.34 with label r ≈ 0 under np_conv4. So the "CAV sees what DE misses"
  contrast survives; ADAMDEC1 is what is lost. Note that the fibroblast-identity genes now also go negative along this axis.
