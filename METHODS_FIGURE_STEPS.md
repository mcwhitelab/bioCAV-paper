# Extended Online Methods — step lists for figure production

Compiled 2026-09-18 from the driver scripts and `figures/figures.R`. This is a
step inventory for prose writing, not the prose itself. Every step names the
script that performs it so a reader can be pointed at the repository.

Conventions used below:
- `$BIOCAV_REPO` = the `bioCAV` method repository; `$BIOCAV_PAPER_ROOT` = `bioCAV-paper`
  (this analysis repository). Both, plus the external resource paths, are resolved by
  `config/paths.sh`, which every driver sources.
- Paths in **bold** are non-default and must be stated explicitly in the methods,
  because sourcing `config/paths.sh` as-is will *not* reproduce that step.
- `[TK]` marks a number or version the author still needs to supply.

---

## 0. Figure inventory and numbering

Figures were renumbered on 2026-07-27. Several file and variable names are one
generation stale; the mapping below is authoritative. Write the paper in the
**final** column.

| Final | Rendered file | Internal names still in the code |
|---|---|---|
| Fig 1 | *(none — schematic)* | — |
| Fig 2 | `figures/figures/cav_fig2.{pdf,png}` | `fig2` |
| Fig 3 | `cav_fig3.{pdf,png}` | `fig3`; panel objects named `p_2d`, `p_2e_solo`, `p_2f_solo` (Fig-2 legacy) |
| Fig 4 | `cav_fig4.{pdf,png}` | `fig4` |
| Fig 5 (SERK/CIK) | `cav_fig5.{pdf,png}` | data CSVs are `fig6_heatmap.csv`, `fig6_densities.csv`, `fig6_auroc.csv`; exporters `export_fig5E_data.py` |
| Fig 5 B/C (MYC top 20, standalone panel) | `cav_fig5BC_myc_top20.{pdf,png}` | `p_myc`; data `fig5BC_myc_top20.csv`; exporter `myc_mad/pca_sweep_supp/export_fig5BC_myc_data.py` (§7.5 step 12) |
| Fig 6 (single cell) | `cav_fig6.{pdf,png}` | drivers `08_export_fig5_data.sh`, `09_export_fig5_umap_data.sh`; exporters `export_fig5_case_studies.py`, `export_fig6_fib_modules.py` (D lower half) |
| Supp — projection comparison | `fig_supp_projection_comparison.{pdf,png}` | `proj_de_vs_cav.csv`; converged C = 1 CAVs since 2026-10-07 (§7.1 item 4) |
| Supp — fibroblast heatmaps (C = 1 / Wald DE / C = 10⁻⁴) | `fig_supp_fib_heatmaps.{pdf,png}` | data `supp_fib_coexpr_{conv1,de_stat,conv4}_{corr,genes}.csv`; exporter `single_cell/convergence_check/08_fib_module_heatmap.py conv1 de_stat conv4` (since 2026-10-07; §6.4 item 22b) |
| ~~Supp — fibroblast subpopulation~~ | DEMOTED 2026-10-07, superseded by `fig_supp_fib_heatmaps`; block gated by `SUPP_FIB_SUBPOP` in `figures.R`, old outputs in `figures/figures/demoted_2026-10-07/` | — |
| ~~Supp — SERK/CIK no-paralog LOO~~ | DEMOTED 2026-10-07; block gated by `SUPP_NOPARALOG_LOO` in `figures.R`, old outputs in `figures/figures/demoted_2026-10-07/` | `fig6_noparalog_scatter.csv` |
| Supp — PCA dimension (MYC/MAD1) | `fig_supp_pca_sweep.{pdf,png}` | `pca_sweep_cacgtg.csv` (C = 10⁻⁴, CACGTG; since 2026-10-05) |
| Pfam CAV circle pack (placement `[TK]`, §7.6) | `motif_clustering/results/converged/C0.0001_sharedscaler/threshold_cuts/diffusion_match4_nopromote/pfam_circlepack_clanid_diff_match4_nopromote.{pdf,png}` | not rendered by `figures.R`; C = 10⁻⁴ library, decided 2026-10-05 |
| Supp — EC evaluation | `fig_ec_eval.pdf` | built but not currently in the numbered set — confirm whether it ships |

**Figure 1** has no data pipeline: it is a hand-drawn schematic of the CAV
construction and scoring procedure, produced in Illustrator. No script generates it.

Other `draft_fig*.R` scripts in `figures/` (VAV reproduction standalone, Pfam clan
UMAP standalone, ADAMDEC1 module, Appendix-A motif supplemental, Fig-4B candidate
variants) are exploratory drafts. They are **not** invoked by `figures/00_make_figures.sh`
and are not paper figures; some are sourced into the session by `figures.R` only for
shared helper code. Do not write methods for them unless one is promoted.

**Single entry point.** `figures/00_make_figures.sh` sources `config/paths.sh` and
runs `Rscript figures/figures.R`, which reads only flat CSV/JSON from
`figures/figure_data/` and writes every panel. All computation is upstream; the R
script performs no statistics beyond the per-panel Spearman correlations noted in §7.

---

## 1. Shared preamble — protein-sequence CAV construction

Steps 1.1–1.6 are common to Figures 2, 3, 4 and 5. Describe once, then refer back.
(The single-cell library in §7 is an analogous but separate construction on a
different model and is described there.)

### 1.1 Protein language model
- ESM++ large (`ESMplusplus_large`), i.e. **ESM-C 600M**, 36 transformer layers,
  hidden size 1152, loaded from a local HuggingFace directory with
  `trust_remote_code=True`.
- License note: ESM-C 600M is MIT-licensed (Biohub, 2026). The `NOTICE.txt` shipped
  alongside the local weights is stale and states a Cambrian non-commercial licence;
  it does not apply.

### 1.2 Embedding
All embedding uses `mcwlab_utils/hf_embed_new.py` (**the copy at
`/groups/clairemcwhite/claire_workspace/github/mcwlab_utils/`, commit `f733f2cc3cb1`** —
the lab-shared copy at `/groups/clairemcwhite/mcwlab_utils` differs by 43 lines).
Standard invocation:

- layer selection `-l -11` (the 26th of 36 hidden states, matching the layer used in
  the preceding preprint), mean pooling (`--strategy mean` / `-ss mean`), batch size 1,
  `--max_length 2048` (1024 for the human proteome search space in §1.7).
- `--get_aa_embeddings` (per-residue) for anything that will be span-pooled;
  `--get_sequence_embeddings` for whole-sequence means.

### 1.3 Background/negative population
- 10,000 randomly drawn proteins, `neg_10000.fasta`.
- Embedded once with the §1.2 settings → `shared_protein/reference_population/neg_10000.pkl`
  (`shared_protein/00_get_negatives.sh`, one GPU job).
- This single set serves as the negative class for **every** sequence CAV in the paper.

### 1.4 Shared coordinate frame
- `$BIOCAV_REPO/scripts/fit_reference_scaler.py --neg neg_10000.pkl --out scaler_v1.pkl`
  fits one `StandardScaler` on the background population.
- Every CAV is trained and scored in this one frame, which is what makes scores from
  different concepts directly comparable. See §10.1 for the one place this was violated.

### 1.5 Span pooling
- `$BIOCAV_REPO/scripts/prepare_embeddings.py` mean-pools the per-residue embeddings
  over each concept's annotated span, producing the positive matrix.

### 1.6 CAV training
- `$BIOCAV_REPO/scripts/train_cav_from_embeddings.py` → `core/src/train_cav.py`, wrapped
  by `specific_scripts/train_cav_from_span.sh`, separates the span-pooled positives from
  the 10k background negatives in the shared scaler frame.
- **Estimator (current code, 2026-10-01):** scikit-learn `LogisticRegression`, L2 penalty,
  `solver='lbfgs'`, `class_weight='balanced'`, `random_state=42`, fit to convergence
  (`tol = 1e-10`, `max_iter = 100000`), **`C = 1/n_train` by default** (`default_C()` in
  `core/src/train_cav.py`; `--C` on `train_cav_from_embeddings.py` overrides it). An intercept
  is fitted but is **not** part of the concept vector. The report records `C`, `C_rule`
  (`1/n_train` or `user`), `tol`, `n_iter` and `converged`.
- **Why C scales with n.** sklearn minimises ½‖w‖² + C·Σᵢ lossᵢ: the loss is *summed*, so the
  effective regularisation is ~1/(C·n), and balanced class weights still sum to n. A fixed C
  therefore means different regularisation at different training-set sizes. C = 1/n keeps
  C·n = 1. With the shared 10k background, 1/n ≈ 10⁻⁴, which is the value used for the
  libraries. On a small custom set the library value is wrong: on the SERK/CIK design
  (n = 32), C = 10⁻⁴ shrinks the CAV to the class-mean difference and cross-species AUROC
  falls from 1.000 to 0.535 (`interface_search/fastas/pairs/lrr_serkcik/converged_C1e-4/COMPARISON.md`).
  Matching C·n is a sensible default, not a guarantee: on that design C = 1/32 ties on
  cross-species AUROC (1.000) but is lower on LOO (0.770 vs 0.821) and homology-filtered k = 1
  (0.816 vs 0.847). For custom sets, report C and consider choosing it by inner CV.
- **Libraries (decided 2026-09-29, shared scaler):** GO, EC and Pfam were retrained with a
  **fixed C = 10⁻⁴** (not `1/n`; for n ≈ 10,000–11,000 the two differ by <10%), converged
  (`tol = 1e-10`, `max_iter = 100000`):
  - Pfam: `/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_converged_C0.0001_sharedscaler`
    (`distribution/retrain_pfam_converged/`). It clusters Pfam clans better than the old library
    (ARI at 1,326 clusters 0.185 vs 0.163; clan contiguity 64 vs 11).
  - GO/EC: `/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_sharedscaler`
    (`distribution/retrain_go_ec_C1e-4/`; 7,631/7,631, all converged, n_iter median 28, max 39;
    EC 3.4.22.28 trained on positives truncated to 2,046 aa, §1.7).
  - The **old libraries** used C = 1 with sklearn's default `tol = 1e-4`, `max_iter = 1000`,
    and stopped after 12–50 iterations. Their directions therefore depend on the stopping
    rule (tol 1e-5 moves Pfam CAVs to cos 0.46–0.89). The old Pfam library also used one
    scaler per motif, not the shared one. Median cosine between a C = 10⁻⁴ GO/EC CAV and its
    old counterpart is about 0.45, so the swap changes directions substantially.
  - **Figures 2 and 3 use the C = 10⁻⁴ GO/EC library since 2026-10-01** (every panel; §2 and §3
    data status). For Figures 2 and 3 the methods sentence below applies.
  - **The circle pack (§7.6) uses the C = 10⁻⁴ Pfam library since 2026-10-05.**
  - **Still on the old libraries:** Fig 4 (Pfam), the thresholds and
    the Zenodo records. Wherever those appear, the methods must describe the old estimator.
  - Methods sentence (applies to Figs 2 and 3; extend to Pfam once Fig 4 / §7.6 are
    regenerated): "A linear L2-regularised logistic regression (L-BFGS,
    class_weight = "balanced") was trained for each concept to separate positive from negative
    embeddings in the shared reference scaler frame, fit to convergence (tolerance 10⁻¹⁰). The
    inverse regularisation strength was C = 10⁻⁴ for the GO, EC and Pfam libraries
    (≈ 1/n for n ≈ 10⁴ training embeddings). Because scikit-learn sums the loss over samples,
    C was scaled to the training-set size for CAVs trained on custom sets (default C = 1/n)."
- **Interface CAVs are separate designs.** The SERK/CIK CAV (Fig 5E, §5.4) uses C = 1 and a
  scaler fit on its own 32 training genes. Converged C = 1 reproduces it (cos 0.99997, identical
  metrics). The MYC CAV (Fig 5B/C, §7.5 step 12) uses **C = 10⁻⁴, fit to convergence, PCA-128**
  (promoted 2026-10-05; previously C = 1, early-stopped, PCA-128). With one positive and the 10k
  negatives (C·n ≈ 1) the fit is close to the class-mean difference in the 128-PC frame
  (cos 0.995). Without PCA it equals the standardised MYC window (cos 0.9998), i.e. a z-scored
  similarity search, which is why the PCA step is kept. Trial and comparisons:
  `interface_search/fastas/pairs/myc_mad/converged_C1e-4/COMPARISON.md`.
- **Concept vector:** `coef_[0]` divided by its L2 norm (`cav / (||cav|| + 1e-10)`),
  saved as `concept_v1.npy`. It is therefore a unit direction; the classifier's
  intercept and probability calibration are discarded.
- Training diagnostics: stratified k-fold cross-validated AUROC is recorded per concept
  alongside the vector.
- **No dimensionality reduction by default.** Every CAV in the main figures is trained in the
  full 1152-d standardized space, **except the MYC CAV in Fig 5 B/C, which uses 128 PCA
  components** (§7.5 step 12). The single-cell library of Fig 6 is a separate construction (768-d Geneformer, PCA-128, C = 1; §6.1 items 5 and 5b). PCA before the logistic regression is a user-selectable option,
  not part of the default pipeline. Its effect is characterised in the supplemental sweep (§7.5). Earlier versions of the pipeline used a PCA-128 step; the current
  `train_cav_from_embeddings.py` has no `--pca-pkl` flag (`$BIOCAV_REPO/docs/internal_cluster_paths.md`).
- **Scoring** a new window = standardise with `scaler_v1.pkl`, then take the dot product
  with the unit concept vector — a signed projection onto the concept direction, not a
  probability. No figure panel z-scores CAV scores (the Fig 2 pos/neg z-score density was
  dropped 2026-10-06; the Fig 6 module scores z-score gene expression, not CAV scores).
- Note `train_cav_from_span.sh` deletes the pooled positive embeddings at its step 5;
  retraining a concept therefore requires re-embedding.

### 1.6b The LLR confidence statistic
Several panels gate on an LLR; define it once.

- For each concept, `eval_go_temporal.py` / `eval_ec.py` fit two univariate Gaussians to
  that concept's CAV projection scores on its **held-out proFAB test lists** (§1.8):
  N(`pos_mean`, `pos_std`) over `random_positive_test_max1000.span` and
  N(`neg_mean`, `neg_std`) over `random_negative_test_max1000.span` (`POS_SPAN` / `NEG_SPAN`
  in both scripts). *(Corrected 2026-10-01: this item previously said training positives and
  the 10k background negatives. The code uses neither; see also the docstring note in §2.1.)*
- The Gaussians are fit with whichever CAV is being evaluated, so the LLR gate changes when the
  library changes. For the C = 10⁻⁴ library it was recomputed from the cached test embeddings.
  Recomputing it with the old CAVs reproduces the old shipped LLR to file rounding (max
  |Δ| 5e-7 GO, 5e-5 EC).
- For a query score *x*, **LLR = log N(x; pos_mean, pos_std) − log N(x; neg_mean, neg_std)**
  (`scipy.stats.norm.logpdf`), i.e. a natural-log likelihood ratio of the positive over
  the background model. LLR > 0 means the score is more likely under the concept's own
  positive distribution than under background.
- Thresholds used: **LLR ≤ 0** is the "CAV made no confident call" rule in Fig 2B and
  Fig 3A; `--llr-threshold 2.3` (≈ ln 10, i.e. 10:1 odds in favour of the positive model)
  is the stricter threshold passed to `compare_tool_temporal.py` and `compare_ec_tools.py`
  for the tool-comparison recall figures. The script default is 1.0; 2.3 is set
  explicitly in `go_eval/03_compare_to_outside_tools.sh` and `ec_eval/04_compare_ec_tools.sh`.

### 1.7 Concept libraries used in the paper

| Library | n concepts | Source labels | Built by |
|---|---|---|---|
| GO | 6,677 | proFAB GO sets, parts 1–10, `random` positive train split | `shared_protein/02_make_cav_commands.sh` → `go_eval/COMMANDS_files/cav_GO_COMMANDS*.sh`, run as SLURM arrays |
| EC | 952 | proFAB EC sets, parts 1–4, `random` positive train split | same, → `ec_eval/COMMANDS_files/cav_EC_COMMANDS*.sh`, `ec_eval/SBATCH_files/array_EC*.sbatch` |
| Pfam motif | 27,472 (27,459 trained, 13 failed) | Pfam domain spans | external run; **`/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all`**, layer 25 only |

GO and EC above are the **old** libraries. **Figures 2 and 3 use the C = 10⁻⁴ library** (§1.6):
**6,678 GO and 953 EC CAVs** at `/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_sharedscaler`.
That is the old 6,677 + 952 plus two terms the old library lacked:
- **GO:0032436:** no old CAV existed.
- **EC 3.4.22.28 (picornain 3C):** all 76 positives are picornavirus polyproteins of
  2,052–2,432 aa, which the 2,048-aa embedding limit drops. Here each was truncated to its
  first 2,046 aa (`hf_embed -t 2046`), which keeps the 3C protease region at about residues
  1,500–1,750 (`distribution/retrain_go_ec_C1e-4/03_ec_3-4-22-28_truncated.sh`). State the
  truncation in the methods.
- **Positive counts:** the C = 10⁻⁴ positives are the same proteins as the old library's.
  n_positive matches for all 7,629 shared terms.
- **Pipeline wiring:** the figure pipelines read CAVs through a symlink overlay of proFAB
  (`/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay`,
  `distribution/retrain_go_ec_C1e-4/04_build_overlay.py`) that differs only in `concept_v1.npy`.

Dataset preparation for GO/EC:
1. `shared_protein/01_download_profab_datasets.sh` — download and extract the proFAB
   GO (10 parts) and EC (4 parts) benchmark archives from `profab.kansil.org`.
2. `shared_protein/format_proFAB.sh` → `specific_scripts/format_proFAB.py --max-n 1000 --seed 42`
   — convert proFAB's index files into accession lists (`.span`), subsampling to at most
   1,000 per list. See §1.8 for exactly which of proFAB's splits this produces and which
   one is used.
3. `shared_protein/02_make_cav_commands.sh` — emit one `train_cav_from_span.sh` command
   per concept; execute as SLURM job arrays.

### 1.8 Which proFAB split is used — be specific, proFAB ships three

proFAB provides **three independent split schemes** per concept, and the methods must
name one. Each scheme is stored as a pair of file types:

- `{group}_{polarity}_gaac.txt` — one row per protein, column 1 the UniProt accession
  followed by its precomputed GAAC features. This file is the accession dictionary.
- `{group}_{polarity}_{split}_indices.txt` — row numbers indexing into that gaac file.

with `group ∈ {random, target, temporal}`, `polarity ∈ {positive, negative}` and
`split ∈ {train, test, validation}` (validation exists only for `temporal`).

Coverage on disk is uneven, which is one reason the choice matters:

| Scheme | GO concepts | EC concepts |
|---|---|---|
| `random_` (random train/test partition) | 6,678 / 6,678 | 953 / 971 |
| `target_` (target-based partition) | 4,577 | 291 |
| `temporal_` (time-based, with a validation split) | 1,597 | 294 |

`format_proFAB.py` auto-detects every group present and writes
`{group}_{polarity}_{split}_max1000.span` — a plain one-accession-per-line file — for
all of them (`--max-n 1000`, `--seed 42`, so the subsample is reproducible). All three
schemes therefore exist on disk. Only one is used:

- **CAV training positives: `random_positive_train_max1000.span`**, and the sequences
  fetched for it, `random_positive_train_max1000.span.fasta`. This is the *only* proFAB
  span file any script in the repository consumes for training — verified by grep across
  all `.sh` and `.py` in the paper repository, which returns 23 references, every one of
  them to `random_positive_train_max1000`, and none to any `target_` or `temporal_`
  span. The trained CAV is written beside it in `random_positive_train_max1000_cav/`.
- **CAV training negatives are not from proFAB at all** — they are the single shared 10k
  background population of §1.3. proFAB's own `*_negative_train_*` lists are unused.
- **Evaluation negatives: `random_negative_test_max1000.span`** (hard-coded as `NEG_SPAN`
  in both `eval_go_temporal.py` and `eval_ec.py`), i.e. held out from the training split.
- **Evaluation positives:** for EC, `random_positive_test_max1000.span` (hard-coded
  `POS_SPAN`). For the GO temporal evaluation in Figure 2, the positives are **not** a
  proFAB split — they are a separately downloaded QuickGO annotation set (§2.1).

**A naming collision worth one sentence in the paper.** The Figure 2 evaluation is
described as "temporal", and proFAB also ships a `temporal_` split. They are unrelated.
The paper's temporal hold-out is the QuickGO set of annotations dated after 2021-12-31;
proFAB's `temporal_` files were generated by `format_proFAB.py` but never read. Readers
familiar with proFAB will otherwise assume the `temporal_` split was used.

Caution for the methods text: the smaller 7,692-motif Pfam library
(`$AHMAD_CAV_OUTPUTS_DIR`) is the **default** value in `config/paths.sh`, but no
current paper figure uses it — the 27,472-motif library above is used throughout, and
its path is the commented-out line in `paths.sh`. State the 27k path explicitly.

---

## 2. Figure 2 — GO evaluation (temporal hold-out vs. DeepGo-SE)

Layout (double column since 2026-09-29; relettered 2026-10-06): A blank (reserved, full
width); B rank-tier composition by sequence-family count, full width (§2.4 B, "current
panel B"); C AUC violins | D AUPR violins. Output 7.2 × 7.9 in.

**Panel relettering (2026-10-06).** The former panel B, a pooled positive/negative density of
z-scored CAV scores, was dropped to keep the methods lean: it was the only z-scored panel,
it pooled raw scores across GO terms (unlike the per-term AUC/AUPR panels), and its
positive set (Fig 2 rank pairs) and negative set (proFAB test negatives) matched neither
of the other panels. The panels after it moved up one letter: old C → **B**, old D/E →
**C/D**. Every Fig 2 reference in this file now uses the new letters, including dated
notes written before the change. Script and file names (`export_fig2C_families.py`,
`p_2a`/`p_2b`/`p_2c`) keep the old letters. Pre-change copies:
`figures/figure_data/_backup_pre_fig2_dropB_2026-10-06/`. The source CSVs
(`figure_data/temporal_pos_neg_density_*.csv`) are left in place but no longer read.

**Data status (2026-10-01): all Fig 2 panels use the C = 10⁻⁴ converged GO library** (§1.6, §1.7).
- **Unchanged:** the known-positive-excluded pairs, the term set, the candidate filters and
  DeepGo-SE.
- **Changed:**
  - the CAV vectors;
  - the LLR gate, recomputed with the new CAVs (§1.6b);
  - one extra Fig 2B candidate, GO:0032436 (BP 4,298).
- **Driver:** `figures/rerun_fig23_C1e-4.sh`, which runs each panel's own pipeline against
  the proFAB overlay. Its steps:
  - `06_eval_exclude_known_positives.sh` (B, D, E);
  - `go_eval/rank_C1e-4/rerank_go.py` (C);
  - BLAST extended to GO:0032436 by `go_eval/rank_C1e-4/extend_blast.sh`. All 280 of its
    training positives are already in the original DIAMOND database, which was searched with
    `--max-target-seqs 0`. So the original hits are complete for this term, and extending
    only adds its member rows. No new search is needed, and this is exactly equivalent to
    having had the term in the original search.
- **Regression check:** `LIB=old` reproduces the previous live files exactly (ranks
  identical, LLR to file rounding, BLAST ranks identical to `noobsolete/`).
- **Archive:** the old-library data is in
  `figures/figure_data/archive_fig23_pre_C1e-4_2026-10-01/` (README) and the old renders are
  in `figures/previews/archive_rendered_pre_C1e-4_2026-10-01/`.
- **Term set:** GO:0032436 is a candidate only. The evaluated-term set stays pinned to
  `static/fig2_go_terms_{ont}.txt` (136 / 256 / 162), as for the other CAVs trained after
  the original eval.
- **Numbers:** those in §2.4 not marked "C = 10⁻⁴" describe the old library.

**Data status (2026-09-29, superseded by the library swap above): all Fig 2 panels use the known-positive-excluded pairs.**
`figures/figure_data/temporal_*` and `go_specificity_ranks_{ont}.csv` were replaced by
the staged `figure_data_staging/noknown/` files. `go_specificity_ranks_allterms_{ont}.csv`
now comes from `go_eval/blast_baseline/noobsolete/` (strict `_pf` filter, obsolete terms
excluded). The pre-exclusion files are archived in
`figure_data/archive_fig2_pre_noknown_2026-09-29/` (with a README describing how to
restore them). Numbers in §2.1 and §2.4 C that are marked "superseded" describe the
archived data.

### 2.1 Validation set — the post-2021 GO annotations

**Source files.** Three aspect-specific QuickGO GPAD downloads, all pulled
**2026-05-27**, in `go_eval/static/go_annotations/`:

| File | QuickGO aspect | Download limit | Rows | Original download |
|---|---|---|---|---|
| `mf.gpad` | molecular_function | 500,000 | 498,025 | `QuickGO/QuickGO-annotations-20260527_mf.gpad` |
| `bp.gpad` | biological_process | 115,000 | 113,139 | `QuickGO/QuickGO-annotations-1779883181548-20260527.gpad` |
| `cc.gpad` | cellular_component | 270,000 | 262,241 | `QuickGO/QuickGO-annotations-1779891376684-20260527.gpad` |

(The working copies are byte-identical to the QuickGO originals, verified by md5. The
`.txt` and `.tsv` files in `QuickGO/` are the same annotations in QuickGO's other export
formats and are not read by any script.) Each was downloaded with the QuickGO filter
`evidenceCode=ECO:0000314&evidenceCodeUsage=descendants` — see the caveat below, it is
not what it sounds like. Annotation dates in the files span 1999-01-01 to 2026-04-28; the
GPAD headers record ECO release 2026-04-05 and GO release 2026-04-27.

**Gold-standard construction** (`load_gold_standard` in
`$BIOCAV_REPO/specific_scripts/eval_go_temporal.py`), in order:
1. Keep only `DB == "UniProtKB"` rows (drops IntAct CPX-XXXX complexes and other DBs,
   which cannot be fetched from UniProt).
2. **`--date-cutoff 20211231`** — keep annotations dated strictly after 2021-12-31.
   This is the temporal hold-out: MF 334,988 rows (38,502 proteins, 3,422 terms);
   BP 31,462 (13,896 proteins, 4,567 terms); CC 78,837 (40,495 proteins, 1,799 terms).
3. Reduce to unique (protein, GO term) pairs.
4. Normalise accessions to bare UniProt form — strip a `UniProtKB:` prefix, then strip
   any `-PRO_XXXXXXX` chain or `-2` isoform suffix.
5. **`--min-n 5`** — keep only GO terms with at least 5 unique post-cutoff proteins.
6. Keep only GO terms that have a trained CAV (i.e. a proFAB concept directory).

**What survives** — this is the actual evaluation set, and the attrition is severe:
*(Superseded: Fig 2 now uses the known-positive-excluded run — see §2 data status, 2026-09-29.)*

| | MF | BP | CC |
|---|---|---|---|
| Validation proteins (`val_proteins.span`) | 2,292 | 4,019 | 6,100 |
| GO terms evaluated | 226 | 395 | 258 |
| (protein, term) pairs scored | 2,472 | 4,819 | 7,118 |

**Caveat that must be stated — the MF download is dominated by protein binding.**
ECO:0000314 is IDA, but `evidenceCodeUsage=descendants` also pulls in ECO:0000353 (IPI,
inferred from physical interaction), which is 385,001 of the 498,025 MF rows. As a
result **GO:0005515 "protein binding" alone is 89.6% of all post-2021 MF annotations**
(300,103 of 334,988). It contributes nothing to the evaluation only because proFAB has
no `GO_0005515` set and therefore no CAV, so step 6 removes it. That accident — not a
design choice — is why MF falls from 3,422 candidate terms to 226 evaluated. Say so
explicitly; a reviewer who recomputes the gold standard without step 6 will get a
completely different, protein-binding-dominated benchmark.

**Per-protein artefacts**, cached once per ontology in
`go_eval/results/temporal_eval_{mf,bp,cc}/eval_embeddings/`:
- `val_proteins.span` — the accession list written by the script.
- `val_proteins.span.fasta` — sequences retrieved from UniProt.
- `val_proteins.span.fasta.pkl` (+ `.seqnames`, `.description`) — ESM-C embeddings, §1.2
  settings. This is the file `compare_tool_temporal.py` reads via `--val-pkl`.
- `val_proteins.span_esm.pkl` — an earlier embedding cache, superseded.

**Per-term background files.** For each evaluated GO term, the LLR of §1.6b needs two
Gaussians, fitted from the proFAB `random` split's held-out lists (§1.8):
- `random_positive_test_max1000.span` → the positive distribution (`pos_mean`, `pos_std`).
- `random_negative_test_max1000.span` → the background distribution (`neg_mean`, `neg_std`).

Note for anyone reading the code: `eval_go_temporal.py`'s docstring says the positive
Gaussian is fitted on "the per-term training positives", but `POS_SPAN` is the **test**
positives. The code is right and the docstring is wrong — describe the test positives.

**DeepGo-SE predictions on the same proteins**, in `go_eval/outside_tools/`:
`val_proteins.span_preds_{bp,cc,mf}.tsv` (positives) and
`neg_test_deepgose_preds_{bp,cc,mf}.tsv` (the proFAB negatives), concatenated by
`02_format_outside_tool_preds.sh` into `deepgose_preds_{bp,cc,mf}.tsv`. The
`val_proteins_{bp,cc,mf}.deepgo.csv` files are the dropped DeepGOWeb predictions (§2.3
item 8) and are not used.

7. GO ontology structure from `go_eval/static/go_annotations/go.obo`, GO release
   2026-04-27 per the GPAD headers.
8. Negatives for the AUC/AUPR calculations are the proFAB per-concept
   `random_negative_test_max1000.span` proteins — the held-out negatives of the same
   `random` split the CAV positives come from (§1.8).

### 2.2 CAV evaluation
4. `eval_go_temporal.py` scores every validation protein against every trained GO CAV,
   standardising with the shared `scaler_v1.pkl`, → `eval_temporal_results.tsv`.
5. `go_eval/01_summarize_temporal_eval.sh` → `specific_scripts/summarize_temporal_eval.py`,
   run once per ontology (`--label mf|bp|cc`), writes the per-term summary and, via
   `--figure-data-dir`, exports `temporal_pos_neg_density_*.csv` (no longer used: its panel was dropped
   2026-10-06),
   `temporal_pr_curves_*.csv`, `temporal_per_term_summary_*.csv`,
   `temporal_protein_pairs_*.csv`.

### 2.3 Comparison tool
6. **DeepGo-SE** is the sole comparator, run locally on our own hardware (not via a
   web service). Full configuration:
   - **Container:** `/xdisk/clairemcwhite/deepgose_latest.sif` (the DeepGO2/DeepGO-SE
     image, pulled 2026-06-09), executed with `apptainer exec --nv`.
   - **Model data:** `/xdisk/clairemcwhite/data` bound to `/workspace/deepgo2/data`,
     supplying `go.obo`, `go-plus.norm` and the per-ontology `terms.pkl` and model
     checkpoints.
   - **Command, complete:**
     `python predict.py -if <fasta> --device cuda --batch-size 1`.
     Those are the only three flags passed; **every other parameter is the DeepGO-SE
     default.**
   - **Prediction threshold: 0.1, the `predict.py` default** (`--threshold/-t`,
     `default=0.1`), not overridden. Verified two ways: the click option in
     `predict.py`, and empirically — the minimum score in every output file, positives
     and negatives across all three ontologies, is exactly 0.100, with zero rows below
     it. *(Threshold 0 was the separate DeepGOWeb API run, §2.3 item 8 — a different
     tool, dropped from the paper. Don't carry that number over to DeepGo-SE.)*
   - **Architecture, as run:** per-ontology ensembles of `DeepGOModel` over 2560-d ESM2
     features extracted by the container's own `extract_esm` — 6 checkpoints each
     (`deepgozero_esm_plus_{n}.th`), averaged: MF {0,1,2,5,6,8}, BP {2,5,6,7,8,9},
     CC {1,3,4,5,6,7}. All three ontologies are predicted in one pass per FASTA.
   - **Inputs:** the validation proteins
     (`go_eval/SBATCH_files/deepgose_val_proteins.sbatch`, over
     `temporal_eval_*/eval_embeddings/val_proteins.span.fasta`) and the proFAB negative
     test proteins (`go_eval/SBATCH_files/deepgose_negatives.sbatch`, over every
     `random_negative_test_max1000.span.fasta`). Both scoring runs are needed: the
     AUC/AUPR comparison requires the tool scored on the same negatives as the CAVs.
     These two files are repository copies of the scripts that actually ran, which
     lived at `/xdisk/clairemcwhite/test_deepgose_val.sbatch` and
     `/xdisk/clairemcwhite/test_deepgose.sbatch` and hard-coded a previous repository
     layout; the copies take their paths from `config/paths.sh` and carry a provenance
     header. They have not been re-run — the original outputs are md5-identical to the
     files shipped in `go_eval/outside_tools/`. Container, model-data and torch-cache
     locations can be overridden with `DEEPGOSE_SIF`, `DEEPGOSE_DATA` and
     `TORCH_CACHE`.
   - `go_eval/02_format_outside_tool_preds.sh` gunzips and concatenates positives and
     negatives into `outside_tools/deepgose_preds_{bp,cc,mf}.tsv` (3 columns:
     protein, GO term, score).
7. **Consequence of the 0.1 threshold, state it once.** DeepGo-SE emits only
   (protein, term) pairs scoring ≥ 0.1, so its output is sparse. In
   `compare_tool_temporal.py` any (protein, term) pair with no emitted row is
   **filled with 0.0** before AUC/AUPR (`merged["tool_score"].fillna(0.0)`, and
   `pid_to_score.get(p, 0.0)` for negatives). Because positives *and* negatives were
   both run through the same threshold, this zero-fill is symmetric and does not bias
   the comparison — unlike the DeepGOWeb run, which had no negatives and was dropped
   partly for that reason (§2.3 item 9). Still worth one sentence, since a reviewer will
   notice the tool has no scores below 0.1.
8. `go_eval/03_compare_to_outside_tools.sh` → `specific_scripts/compare_tool_temporal.py`,
   per ontology, with `--llr-threshold 2.3` and the shared scaler. Exports
   `temporal_tool_comparison_*.csv` and `go_specificity_ranks_*.csv` to `figure_data/`.
9. **DeepGOWeb/DeepGOPlus was evaluated and then dropped** (2026-07-12) and must not
   appear in the figure or its legend. Reasons, in order of weight: (i) DeepGOPlus
   combines a CNN with a DIAMOND sequence-similarity component that transfers
   annotations from homologous training proteins, so it is not a
   sequence-representation-only method and is not comparable on the question the
   figure asks; (ii) batch API scores (v1.0.27, threshold 0) sat on a different scale
   from the web interface and could not be reconciled; (iii) the batch covered only
   validation proteins, with no negatives, so its AUC/AUPR were zero-fill inflated.
   The recipe is retained, commented out, in `03_compare_to_outside_tools.sh`, and its
   outputs are quarantined in `figure_data/deepgoweb/`.

### 2.4 Panels
- **Current panel B (2026-09-29; was C until 2026-10-06)** (`p_2a`, data `figure_data/go_rank_families_{ont}.csv`,
  exporter `go_eval/blast_baseline/export_fig2C_families.py`): 100% stacked rank-tier bars
  for CAV, DeepGo-SE and a DIAMOND nearest-neighbour baseline, grouped by ontology (rows:
  MF / BP / CC blocks, one bar per method). Columns split by the number of sequence
  families among the true term's CAV training positives: <20 / 20–59 / ≥60. The cutoffs
  come from the screen in `go_eval/blast_baseline/family_cutoff_screen.py`
  (`preview_family_screen.png`): all methods drop steeply up to about 15–20 families,
  plateau at 20–59, and the CAV lead over BLAST disappears from about 60 families.
  - Pairs: known-positive-excluded, 1,223 / 2,739 / 4,052 (MF / BP / CC). One exact
    duplicate (protein, term) row per ontology is dropped, from two annotation IDs that
    normalise to one accession, e.g. P0DTC2 and its -PRO_ chain.
  - Families: `go_eval/diversity/term_diversity.py`, i.e. DIAMOND all-vs-all of the
    term's `random_positive_train_max1000` FASTA, with an edge at ≥30% identity over
    ≥50% of the shorter sequence, and connected components counted as families.
  - Ranks: `*_rank_ont_filt_pf` over non-obsolete same-ontology CAVs
    (C = 10⁻⁴ library: MF 1,114 / BP 4,298 / CC 660; old library: BP 4,297), minus the
    true term's ancestors and the protein's GPAD and proFAB-positive terms (both
    ancestor-closed).
  - BLAST: best DIAMOND bitscore (`--more-sensitive`, e ≤ 1e-3, self-hits removed) to
    each term's CAV training positives, tie-averaged rank (`blast_ranks.py`). It is ranked
    over the same candidates as CAV and DeepGo-SE, including GO:0032436 for the C = 10⁻⁴
    library (§2 data status).
  - "No prediction": CAV LLR ≤ 0; DeepGo-SE term not emitted (score < 0.1); BLAST no hit
    to the true term's training positives.
  - **Result, C = 10⁻⁴ library (current Fig 2B)** (Rank 1 / top 10 %, gated), from
    `figure_data/go_rank_families_*.csv`:

  | | Families | Pairs | CAV | DeepGo-SE | BLAST | CAV no prediction |
  |---|---|---|---|---|---|---|
  | MF | all | 1,223 | 17.8 / 41.2 | 10.1 / 22.7 | 11.4 / 34.3 | 33.4 |
  | MF | <20 | 262 | 37.8 / 66.0 | 18.7 / 37.4 | 21.0 / 50.4 | 26.3 |
  | MF | 20–59 | 745 | 14.4 / 32.1 | 7.9 / 19.3 | 8.6 / 28.6 | 38.1 |
  | MF | ≥60 | 216 | 5.6 / 42.6 | 6.9 / 16.7 | 9.3 / 34.7 | 25.5 |
  | BP | all | 2,739 | 5.5 / 20.9 | 1.5 / 6.2 | 2.9 / 12.5 | 38.1 |
  | BP | <20 | 241 | 12.0 / 28.6 | 0.0 / 1.7 | 0.4 / 15.4 | 47.3 |
  | BP | 20–59 | 1,728 | 6.0 / 24.7 | 1.3 / 6.7 | 3.7 / 13.4 | 36.2 |
  | BP | ≥60 | 770 | 2.3 / 10.0 | 2.6 / 6.8 | 1.8 / 9.6 | 39.6 |
  | CC | all | 4,052 | 8.0 / 26.5 | 3.5 / 15.7 | 8.3 / 29.2 | 33.3 |
  | CC | <20 | 315 | 12.4 / 31.7 | 7.3 / 30.2 | 5.4 / 20.0 | 57.8 |
  | CC | 20–59 | 1,685 | 9.2 / 30.1 | 1.1 / 15.8 | 8.5 / 26.2 | 31.9 |
  | CC | ≥60 | 2,052 | 6.3 / 22.7 | 4.9 / 13.4 | 8.6 / 33.1 | 30.8 |

  - **What the swap did to panel B.**
    - **Totals:** the overall rates barely move (old → new rank 1 / top 10: MF 16.4 / 42.2 →
      17.8 / 41.2, BP 5.7 / 20.5 → 5.5 / 20.9, CC 8.1 / 28.6 → 8.0 / 26.5).
    - **Individual pairs:** these churn heavily. Only 5–17% of pairs keep the same CAV rank.
      Rank-1 pairs lost vs gained: MF 61 / 79, BP 80 / 73, CC 131 / 124; McNemar p 0.15–0.71.
      These counts compare the shipped term set against the C = 10⁻⁴ CAVs on that same set;
      adding GO:0032436 changes CAV ranks only for BP pairs.
    - **The exception is BP <20 families:** rank 1 falls 18.7 → 12.0% (20 lost, 4 gained,
      p = 0.002, n = 241).
    - **No prediction:** CAV "no prediction" falls about 10 points in every ontology (MF 43.8 →
      33.4, BP 47.3 → 38.1, CC 43.6 → 33.3%). Most of the newly gated pairs land at
      Rank > 10, so the change is more coverage, not better ranking.
    - **Methods sentence:** CAV leads DeepGo-SE on rank 1 and top 10 in every ontology. CAV
      leads BLAST in MF and BP, but **not in CC** (8.0 / 26.5 vs 8.3 / 29.2), so don't write
      that CAV leads BLAST throughout.
    - **Callouts:** re-picked on the C = 10⁻⁴ library in
      `go_eval/blast_baseline/callouts_C1e-4/`, with the old selection rules and homology
      recipe recovered and calibrated (see its README). Strict, no-detectable-homology
      callouts: MF 6, BP 11, CC 34. Of the old MF ones only Q2RVS1 (EncFtn, ferroxidase)
      survives. The old `callouts/` and `bp_cav_only_rank1_cases.csv` describe the old
      library; don't cite them.
  - Result, old library (archived; was the Fig 2 rank panel until 2026-10-01):

  | | Families | Pairs | CAV | DeepGo-SE | BLAST |
  |---|---|---|---|---|---|
  | MF | <20 | 262 | 37.8 / 64.9 | 18.7 / 37.4 | 21.0 / 50.4 |
  | MF | 20–59 | 745 | 12.2 / 36.1 | 7.9 / 19.3 | 8.6 / 28.6 |
  | MF | ≥60 | 216 | 4.6 / 35.6 | 6.9 / 16.7 | 9.3 / 34.7 |
  | BP | <20 | 241 | 18.7 / 36.5 | 0.0 / 1.7 | 0.4 / 15.4 |
  | BP | 20–59 | 1,728 | 6.2 / 23.5 | 1.3 / 6.7 | 3.7 / 13.4 |
  | BP | ≥60 | 770 | 0.6 / 8.8 | 2.6 / 6.8 | 1.8 / 9.6 |
  | CC | <20 | 315 | 15.2 / 27.6 | 7.3 / 30.2 | 5.4 / 20.0 |
  | CC | 20–59 | 1,685 | 10.3 / 35.4 | 1.1 / 15.8 | 8.5 / 26.2 |
  | CC | ≥60 | 2,052 | 5.3 / 23.2 | 4.9 / 13.4 | 8.6 / 33.1 |
- **Previous rank panel** (archived data): for each (validation protein, true GO term) pair, the rank of the true
  term among the candidate GO CAVs of **its own ontology**, binned as Rank 1 / 2–3 /
  4–10 / >10 / no prediction, shown as 100% stacked bars faceted by method
  (CAV, DeepGo-SE) × ontology. n pairs: MF 2,472, BP 4,819, CC 7,118.
  *(Superseded: Fig 2 now uses the known-positive-excluded run — see §2 data status, 2026-09-29.)*

  **Computed by** `go_eval/rank_all_terms/compute_ranks_all_terms.py` →
  `figure_data/go_specificity_ranks_allterms_{mf,bp,cc}.csv`, column
  `cav_rank_ont_filt` / `tool_rank_ont_filt` (2026-09-26; replaces the
  `make_rank_scatter()` ranks in `go_specificity_ranks_*.csv`). Steps:
  1. **Candidate set = all trained, non-obsolete GO CAVs in the pair's ontology**,
     namespace taken from `go.obo` (2026-04-27 release): **MF 1,114, BP 4,297, CC 660**
     (6,071 of 6,677). Excluded: 484 CAV terms marked `is_obsolete: true` (MF 47, BP 390,
     CC 47) and 122 whose IDs are not in this `go.obo` at all. Report the per-ontology
     counts in the legend. *(Corrected 2026-09-29: this item previously gave MF 1,161 /
     BP 4,687 / CC 707 and said the 122 absent IDs were "obsolete". Obsolete terms keep
     their `namespace:` line in `go.obo`, so all 484 were candidates until
     `parse_obo(drop_obsolete=True)` was added. The live and staged
     `go_specificity_ranks_allterms_*.csv` still include them until
     `compute_ranks_all_terms.py` is re-run; see "Obsolete terms" below.)*
  2. **Filtered rank.** Removed from each pair's candidates: (i) the true term's
     `is_a`/`part_of` ancestors, and (ii) every other GO term the protein is annotated
     with in the three GPADs at any date (pre- and post-cutoff), closed over ancestors.
     A correct-but-less-specific ancestor, or another correct annotation, outranking the
     true term is not an error. Median terms removed per pair: MF 33, BP 38, CC 24.
  3. **CAV score** = `scaler_v1.pkl`-standardised sequence embedding · unit concept
     vector (raw projection, not LLR or z), using the **retrained** GO CAVs on disk
     (all 6,677 in the shared frame; §10.1).
  4. **DeepGo-SE score** = its emitted score, 0.0 where it emitted nothing. Each
     validation FASTA was predicted in all three ontologies, so predictions for all three
     are pooled per protein before restricting to the same candidate set as the CAV —
     both methods are ranked over identical candidates.
  5. Rank = `1 + #{candidates scoring strictly above the true term}`, so ties take the
     optimistic rank.
  6. "No prediction" is unchanged: LLR ≤ 0 for CAV, term not emitted (score < 0.1) for
     DeepGo-SE — two different definitions, state both. The pair list and both gates are
     carried over from `go_specificity_ranks_{ont}.csv`, so **the LLR gate still comes
     from the July CAVs** and refreshes when `eval_go_temporal.py` is re-run.
     *(Superseded: Fig 2 now uses the known-positive-excluded run — see §2 data status, 2026-09-29.)*

  *(Superseded: Fig 2 now uses the known-positive-excluded run — see §2 data status, 2026-09-29.)*
  Result (Rank 1 / top 10): CAV MF 32.5% / 67.8%, BP 12.7% / 36.2%, CC 23.1% / 52.8%;
  DeepGo-SE MF 16.9% / 65.9%, BP 0.6% / 5.4%, CC 3.6% / 48.3%.

  **Why the previous definition was replaced.** `make_rank_scatter()` filtered candidates
  with `if go_id not in go_terms: continue`, where `go_terms` is only the terms evaluated
  in that ontology's run — a denominator of **226 / 395 / 258** (`n_go_terms` column), in
  which the true term competed only against other terms that happened to reach the gold
  standard. Rank-1 figures from that version (e.g. DeepGo-SE MF 74%) are **not**
  comparable to the current panel.

  **Why same-ontology rather than all 6,677, and why filtered.** Pooling all three
  ontologies compares scores across MF/BP/CC; without filtering, DeepGo-SE is outranked
  by the true term's own ancestors by construction, because its outputs are
  hierarchy-consistent (ancestor ≥ descendant). Even after filtering, many terms that
  outrank the true term for DeepGo-SE are broad, high-prior terms (e.g. "membrane",
  "catalytic activity") that are likely true but absent from the IDA-descendant GPADs.
  This is a property of ranking calibrated per-term probabilities across terms, and a
  reviewer may raise it; state it in one sentence. The unfiltered and all-6,677
  variants are kept as columns in the same CSV (`*_rank_ont`, `*_rank_all`,
  `*_rank_all_filt`), and `figures/draft_fig2C_allterms_variants.pdf` shows all of them
  side by side.

  **Obsolete terms (2026-09-29).** `parse_obo()` in `compute_ranks_all_terms.py` now sets
  the namespace of obsolete terms to None (`drop_obsolete=True`, the default), which
  removes them from every per-ontology mask. `blast_ranks.py` and
  `go_eval/pca_test/pca_test_cc.py` import the same `parse_obo`, so they inherit the
  change. The pooled `*_rank_all*` columns are also restricted to the 6,071 live terms
  (`n_terms_all` = 6,071), and `n_dropped_filt*` counts only live terms. Why it matters:
  obsolete terms have no `is_a`/`part_of` parents, so step 2 cannot remove them. Many are
  broad grouping terms (e.g. "obsolete cellular metabolic process") that DeepGo-SE,
  trained on an earlier GO release, scores highly for most proteins. Left in, they were
  DeepGo-SE's top-scoring competitors for most BP pairs. Checked: the edited scripts
  reproduce the wrapper run `go_eval/blast_baseline/noobsolete/run_noobsolete.py` exactly
  for every `*_ont*` column and for `blast_ranks_*.csv`.

  Effect (gated Rank 1 / top 10 %, `_ont_filt_pf`, known-positive-excluded pairs;
  obsolete included → excluded):

  | | CAV | DeepGo-SE | BLAST (mid-tie) |
  |---|---|---|---|
  | MF | 15.8 / 42.2 → 16.3 / 42.3 | 8.4 / 22.5 → 10.0 / 22.7 | 8.3 / 33.5 → 11.3 / 34.3 |
  | BP | 5.4 / 19.8 → 5.7 / 20.6 | 0.4 / 3.9 → 1.5 / 6.2 | 1.9 / 11.5 → 2.9 / 12.5 |
  | CC | 7.3 / 27.6 → 8.1 / 28.6 | 2.7 / 15.3 → 3.5 / 15.7 | 7.5 / 29.0 → 8.3 / 29.2 |

  **Methods text (draft):**
  > Candidate terms were all trained GO CAVs in the pair's ontology, excluding terms
  > marked obsolete in the GO release used for evaluation (47 MF, 390 BP, 47 CC).
  > Obsolete terms have no `is_a` or `part_of` parents, so the ancestor filter cannot
  > remove them. Many are broad grouping terms (e.g. "cellular metabolic process") that
  > DeepGo-SE, trained on an earlier release in which they were still active, scores
  > highly for most proteins. Left in, they outrank the true term for reasons unrelated
  > to its prediction. Including them lowers every method's rank-1 rate and does not
  > change the comparison (CAV / DeepGo-SE rank 1: MF 15.8 / 8.4% included vs
  > 16.3 / 10.0% excluded; BP 5.4 / 0.4% vs 5.7 / 1.5%; CC 7.3 / 2.7% vs 8.1 / 3.5%).
- **C**, **D** (`p_2b`, `p_2c`; were D, E until 2026-10-06): per-term AUC and AUPR (validation positives vs. proFAB
  test negatives), dodged violins by ontology, median quantile drawn, y fixed to [0,1].

  **Methods text (2026-09-29, reviewer comment "AUROC/AUPR Methods"; numbers are for the
  known-positive-excluded run, now in Fig 2):**
  > Figure 2C and D report AUROC and AUPR separately for each GO term; each term
  > contributes one point to the violins (MF 136, BP 256, CC 162 terms), and no scores are
  > pooled across terms. For a given term, the positives are the validation proteins with a
  > UniProtKB annotation to that term dated after 31 December 2021, excluding any protein
  > that proFAB already lists as a positive for that term in either its training or its
  > test split. proFAB positives are drawn from Swiss-Prot releases up to 2020-01 and
  > include annotations to descendant terms, so a post-cutoff annotation date alone does
  > not make a protein–term pair new. Terms with fewer than five remaining positives were
  > dropped. After retrieval losses, terms have 3–409 positives (median 7–9). The negatives
  > are that term's held-out proFAB negative test proteins (9–77 per term, median 26–29).
  > These are disjoint from the positives and were never used to train the CAV. The CAV
  > score is the projection of the standardised ESM-C embedding onto the unit concept
  > vector, using the shared reference-population scaler. We use this raw projection, not
  > the LLR or a z-score. The DeepGo-SE score is its predicted probability for the term.
  > DeepGo-SE was run on both the positive and negative proteins, and any (protein, term)
  > pair it did not emit (score < 0.1) was assigned 0. AUROC was computed with
  > scikit-learn `roc_auc_score`. AUPR is average precision (`average_precision_score`), a
  > step-wise sum and not trapezoidal area. Both methods were scored on identical positive
  > and negative sets. The chance level for AUPR is each term's positive fraction, not
  > 0.5: median 0.21 (MF), 0.23 (BP) and 0.28 (CC), range 0.06–0.91.

  Code: CAV columns `auc_/aupr_val_vs_test_neg` from `summarize_temporal_eval.py`
  (`compute_go_term_auc/aupr`); DeepGo-SE `tool_auc/tool_aupr` from
  `compare_tool_temporal.py`. Shipped (July) numbers were independently recomputed from
  `eval_temporal_results.tsv` + `{go}_test_neg_scores.tsv` + `deepgose_preds_*.tsv`:
  max |ΔAUC| = 0.0000.

  **Known-positive leakage and its fix (2026-09-29).** The temporal gold standard
  filtered annotation *records* by date, not whether the (protein, term) pair was new. A
  protein already known for a term can gain a new post-2021 record, e.g. from a new
  paper. Of the shipped pairs, 915 MF / 1,566 BP / 2,624 CC (33–37%) had a protein in
  that term's CAV training positives (`random_positive_train_max1000.span`). Of those,
  347 / 456 / 1,445 already had an IDA/IPI annotation to the same term before 2021 in our
  GPADs. The rest come from proFAB's true-path propagation or other evidence codes.
  A further 107 / 164 / 280 were in proFAB test positives, which set the LLR Gaussians.

  Fix: `eval_go_temporal.py --exclude-known-positives` drops pairs whose protein is in
  the term's proFAB random train or test positive span, then re-applies `--min-n 5`.
  `--go-terms-file static/fig2_go_terms_{ont}.txt` holds the term set to the shipped
  226/395/258. Without it, 175+ GO CAVs trained in June 2026, after the original eval,
  enlarge the benchmark to 401/650/402 terms whose proteins mostly have no embeddings or
  DeepGo-SE predictions. Run: `go_eval/06_eval_exclude_known_positives.sh`. It writes to
  `go_eval/results/temporal_eval_{ont}_noknown/` and `figures/figure_data_staging/noknown/`,
  and uses the CAVs on disk (retrained 2026-09-16). Promoted to
  `figures/figure_data/` on 2026-09-29 (see §2 data status).

  | | Terms | Pairs | CAV AUROC shipped → July CAVs excl. → retrained excl. | DeepGo-SE AUROC | CAV AUPR shipped → excl. | DeepGo-SE AUPR |
  |---|---|---|---|---|---|---|
  | MF | 226 → 136 | 2,472 → 1,224 | 0.968 → 0.908 → 0.908 | 0.700 → 0.571 | 0.957 → 0.829 | 0.579 → 0.374 |
  | BP | 395 → 256 | 4,819 → 2,740 | 0.920 → 0.815 → 0.815 | 0.545 → 0.500 | 0.875 → 0.713 | 0.409 → 0.318 |
  | CC | 258 → 162 | 7,118 → 4,053 | 0.961 → 0.867 → 0.866 | 0.688 → 0.556 | 0.966 → 0.840 | 0.708 → 0.463 |

  (Medians over terms. The July-vs-retrained column shows the scaler-frame retrain has no
  effect on these terms; the whole drop is the exclusion.) On the old library CAV AUROC
  exceeds DeepGo-SE on 82% / 89% / 92% of terms.

  **C = 10⁻⁴ library (current Fig 2C/D, 2026-10-01).** Same terms, positives and negatives;
  only the CAV changes, so the DeepGo-SE columns are identical.

  | | Terms | CAV AUROC old → C = 10⁻⁴ | CAV AUPR old → C = 10⁻⁴ | DeepGo-SE AUROC / AUPR | CAV > DeepGo-SE (AUROC / AUPR) |
  |---|---|---|---|---|---|
  | MF | 136 | 0.908 → 0.938 | 0.829 → 0.882 | 0.571 / 0.374 | 87% / 87% of terms |
  | BP | 256 | 0.815 → 0.838 | 0.713 → 0.733 | 0.500 / 0.318 | 89% / 93% |
  | CC | 162 | 0.866 → 0.897 | 0.840 → 0.857 | 0.556 / 0.463 | 94% / 94% |

  The methods-text paragraph above needs no change for the swap. Its term counts, positive
  and negative counts and AUPR chance levels do not depend on the CAV.

  Panel B, filtered same-ontology rank as the figure draws it (Rank 1 / top 10, with
  "no prediction" pairs counted as misses), shipped → excluded:
  CAV MF 32.4 / 66.3 → 15.0 / 41.7, BP 12.2 / 35.1 → 4.9 / 18.7, CC 22.9 / 51.9 → 6.7 / 26.8;
  DeepGo-SE MF 16.5 / 41.7 → 7.9 / 22.1, BP 0.6 / 4.9 → 0.4 / 3.0, CC 3.5 / 26.9 → 2.6 / 14.3.
  CAV leads on both measures in all three ontologies. **Caution:** the console summary
  of `compute_ranks_all_terms.py` and the older "Result" line above report DeepGo-SE
  top-10 *without* the no-prediction gate (e.g. MF 65.9 shipped / 55.1 excluded), which
  counts unemitted terms as ranked. That is not what panel B draws, and mixing the two
  conventions wrongly makes DeepGo-SE look ahead on top-10 in MF and CC.

  Why the panel B rates drop: the removed pairs were the CAV's easiest. On the known-positive pairs the
  CAV had Rank 1 52.7 / 23.4 / 45.2% and top 10 96.4 / 61.5 / 86.7% (MF/BP/CC), with
  no-prediction ≤1.4%. That is what you expect for proteins the CAV was trained on.
  DeepGo-SE also did best on those pairs (they are well-studied proteins), but gained less.

  **proFAB-known terms in the rank filter (checked 2026-09-29).** The `_filt` rank drops
  only terms found in our IDA/IPI-descendant GPADs. After exclusion, 60% of kept pairs
  involve a protein that proFAB lists as a positive for another CAV term the filter
  misses (median 8). `compute_ranks_all_terms.py` now also writes `*_rank_ont_filt_pf`,
  which additionally drops those terms and their ancestors (median 4 / 13 / 1 more per
  pair, MF / BP / CC). It barely matters. Gated Rank 1 / top 10, `_filt` → `_filt_pf`:
  CAV MF 15.0/41.7 → 15.8/42.1, BP 4.9/18.7 → 5.4/19.7, CC 6.7/26.8 → 7.3/27.7;
  DeepGo-SE MF 7.9/22.1 → 8.4/22.5, BP 0.4/3.0 → 0.4/3.9, CC 2.6/14.3 → 2.7/15.3.
  Panel B and the diversity previews still use `_filt`. The LLR "no prediction" gate is also refreshed, so it now comes from
  the retrained CAVs, but that makes no difference: CAV no-prediction (LLR ≤ 0) is
  24.8 / 30.2 / 26.0% shipped → 43.6 / 47.4 / 44.0% with the July gate on the excluded
  pairs → 43.8 / 47.3 / 43.6% with the retrained gate, so the jump is entirely the
  exclusion. Preview: `figures/figure_data_staging/noknown/preview/cav_fig2.{png,pdf}`
  (rendered by `preview_fig2.R`, figures.R lines 1–647 with DATA/OUT redirected).
  DeepGo-SE's own training-set overlap with these proteins is
  unchecked; the exclusion removes the same pairs for both methods.

---

## 3. Figure 3 — EC evaluation and cross-ontology validation

Layout, two rows: A EC rank composition | B EC hierarchy depth; C combined CAV space |
D GO–EC cosine | E rank of annotated EC partner. Output 7.2 × 4.1 in
(double-column width).

**Bottom-row reorder (2026-10-06).** The UMAP moved ahead of the GO–EC concordance pair:
old C (cosine) → **D**, old D (EC-partner rank) → **E**, old E (combined CAV space) → **C**.
Every Fig 3 reference in this file uses the new letters, including dated notes written
before the change. Section numbers keep their pipeline order (§3.3 cosine/rank, §3.4 UMAP).
Pre-change copies: `figures/figure_data/_backup_pre_fig3B_P0DUE1_2026-10-06/`.

**Data status (2026-10-01): every Fig 3 panel uses the C = 10⁻⁴ GO/EC library** (§1.6, §1.7).
- **Library:** 6,678 GO + 953 EC CAVs, including GO:0032436 and EC 3.4.22.28.
- **Driver:** `figures/rerun_fig23_C1e-4.sh`, with each panel's own pipeline pointed at the
  proFAB overlay:
  - `ec_eval/rank_C1e-4/rerank_ec.py` (A);
  - `ec_hierarchy_decay.py` (B);
  - `go_ec_cosine_validation.py` plus the 20,000-value background sample, `rng(0)`, as in
    `plot_go_ec_cosine_validation.py` (C/D);
  - `combined_cav_umap_3d.py --dims 2` (E).
- **Regression check:** `LIB=old` reproduces the previous B and C/D files byte-for-byte and
  the A ranks exactly. E has the same concept set, but UMAP coordinates are not
  bit-reproducible across runs even with `random_state=42`.
- **BLAST:** extended to EC 3.4.22.28 by `go_eval/rank_C1e-4/extend_blast.sh`. All 76 of its
  positives are already in the original EC DIAMOND database, so only member rows are added.
  No EC-Bench test protein hits them.
- **Archive:** old data in `figures/figure_data/archive_fig23_pre_C1e-4_2026-10-01/`.
- **Numbers:** those below not marked "C = 10⁻⁴" describe the old library.

### 3.1 Panel A — EC specificity vs. CLEAN-SupCon

**Current panel A (2026-09-29): leak-free, with a DIAMOND nearest-neighbour baseline.**
- **Data:** `figure_data/ec_ranks_samelevel_noknown.csv`, written by
  `ec_eval/blast_baseline/ec_ranks_noknown.py` → `export_fig3A.py`.
- **Pairs:** those whose protein is already a proFAB positive of the true EC (train or
  test span of that EC's CAV) are dropped, 209 → 138 (L3 161 → 107, L4 48 → 31). The old
  `ec_train_overlap_pairs.tsv` (32 pairs) covered only exact-EC pairs of the gold
  standard and missed most level-3 pairs.
- **UniProt fragments dropped (2026-10-06):** 138 → **124** (L3 107 → 103, L4 31 → 21).
  14 pairs whose protein is a UniProt "fragment" entry are removed: 16–106 aa partial
  sequences, mostly snake-venom PLA2s and proteases plus C0HLV0, C0HLQ6 and C0HLQ4. The
  list covers all 468 EC-Bench test proteins (23 fragments) and comes from
  `ec_eval/static/ec_annotations/EC-Bench/fetch_uniprot_fragments.py` →
  `uniprot_fragments.tsv`. It is applied in `ec_eval/rank_C1e-4/rerank_ec.py` after the
  known-positive filter. The surviving 124 rows are identical to the corresponding 138-pair
  rows. Fragments are only removed from Fig 3A; CAV training sets and the UMAP (Fig 3C) are
  unchanged. Fig 3B still shows C0HLQ6 (§3.2 item 8, open).
- **Candidates:** same-level EC CAVs, minus the protein's other EC-Bench ECs and its
  other proFAB-positive ECs at that level.
- **BLAST:** best DIAMOND bitscore (2.0.9, `--more-sensitive`, e ≤ 1e-3, self-hits
  removed) to each EC CAV's training positives (210,130 proteins), tie-averaged rank.
  "No prediction" means no hit to the true EC's training positives.
- **CAV "No prediction":** LLR ≤ 0. For the C = 10⁻⁴ library the gate is recomputed from
  the new CAVs (§1.6b). The old library's gate was carried from the shipped EC eval.
- **CLEAN-SupCon:** EC-Bench `results_test_100_go_3.csv` (unchanged).
- **Result, C = 10⁻⁴ library, fragments dropped (current Fig 3A, 2026-10-06):**

  | | CAV Rank 1 / top 3 | CAV no prediction | BLAST Rank 1 / top 3 | CLEAN correct |
  |---|---|---|---|---|
  | L3 (103) | 55.3 / 75.7% | 13.6% | 75.7 / 78.6% | 71.8% |
  | L4 (21) | 90.5 / 100% | 0% | 100 / 100% | 90.5% |

  - **All four level-4 CAV "no prediction" pairs were fragments**, as were the other
    level-4 misses apart from the two below.
  - **Level 4, CAV wrong and BLAST right: 2 pairs, both rank 2.** A0A345BJN1 (*Cochliobolus
    lunatus*, squalestatin cluster) and P0DUM0 (*Aspergillus hancockii*, hancockiamide
    cluster) are fungal phenylalanine ammonia-lyases (4.3.1.24, LLR ≈ 324). Both are
    outranked by the histidine ammonia-lyase CAV (4.3.1.3; scores 45 vs 39).
    - **Same family:** PAL and HAL are both MIO-dependent ammonia-lyases (Pfam PF00221 in
      all training positives of both CAVs).
    - **Training sets differ by taxon:** the PAL CAV's 54 positives are 46 plant, 6 fungal
      and 2 bacterial (median 711 aa). The HAL CAV's 181 are 171 bacterial (median 509 aa).
    - **So** the PAL direction is mostly a plant-PAL direction. Fungal PALs, which are
      closer to the HAL-like ancestral enzymes, project at least as strongly onto HAL.
    - **Why BLAST gets them:** it finds close fungal relatives among the 6 (89% and 63%
      identity).
    - **Fig 3B caveat:** A0A345BJN1 is also the top Fig 3B curve. Its level-4 PAL score
      (39.3) is lower than its sibling-EC HAL score (45.1), which that panel does not show.
- **Result, C = 10⁻⁴ library, before the fragment filter (archived 2026-10-06):**

  | | CAV Rank 1 / top 3 | CAV no prediction | BLAST Rank 1 / top 3 | CLEAN correct |
  |---|---|---|---|---|
  | L3 (107) | 55.1 / 74.8% | 15.0% | 76.6 / 79.4% | 72.0% |
  | L4 (31) | 80.6 / 87.1% | 12.9% | 93.5 / 93.5% | 80.6% |

  - **Level 3:** this is the one real gain from the swap, 30.8 → 55.1% rank 1. 31 pairs moved
    to rank 1 and 5 moved off (McNemar p = 1e-5). CAV still trails CLEAN and DIAMOND here.
  - **Level 4:** rank 1 is unchanged. Top 3 falls by 2 pairs and "no prediction" rises from
    1 to 4 of 31 pairs.
- **Result, old library (archived):**

  | | CAV Rank 1 / top 3 | BLAST Rank 1 / top 3 | CLEAN correct |
  |---|---|---|---|
  | L3 (107) | 30.8 / 66.4% | 76.6 / 79.4% | 72.0% |
  | L4 (31) | 80.6 / 93.5% | 93.5 / 93.5% | 80.6% |

- **Old file:** `ec_specificity_ranks_samelevel.csv` is kept in `figure_data/` but no
  longer plotted.

**Panel B leak (found 2026-09-29; FIXED the same day, see §3.2 step 8).** Two of the five example proteins in
`ec_hierarchy_decay.csv` are training positives of CAVs whose scores the panel plots:
A0A452E9Y6 is in the training span of 1.11, 1.11.1 and 1.11.1.7, and P0DPS1 is in the
training span of 2.7.7.7. The other three (C0HLQ6, P0DUB9, A0A345BJN1) are clean at
every level. Swap in leak-free examples before submission.

**Provenance of the EC gold standard.** `test_ec.csv` is a **verbatim upstream EC-Bench
artifact**, not something derived in this project. It is md5-identical
(`afeb9909f8b7edf0a0cfd2ca682a66ad`) to `data/test_ec.csv` in the EC-Bench checkout at
`/xdisk/clairemcwhite/clairemcwhite/EC-Bench` (github.com/dsaeedeh/EC-Bench, commit
`5e98a2c`), where it is git-tracked and unmodified. The same holds for `test_ec.fasta`,
`price.csv`, `price-149.csv` and the `results_test_*_go_3.csv` tool-prediction tables —
all copied byte-for-byte into `ec_eval/static/ec_annotations/EC-Bench/`. Nothing in this
repository writes any of them; `01_reformat_gold_standards.sh` only reads `test_ec.csv`
and emits the long-form `test_ec_long.tsv`. (Within EC-Bench itself `test_ec.csv` is
generated by their own pipeline — `data_preprocessing.sh` then `create_data.sh` →
`code/create_train_test_data.py` — so it is a processed benchmark file, not a raw
database download.)

**EC-Bench is itself a temporal benchmark — state this correctly.** Its splits are
time-based (EC-Bench README, "Summary of Dataset Splits"):

| EC-Bench split | Source | Release |
|---|---|---|
| Pretraining | TrEMBL | 2018-02 |
| Training | Swiss-Prot | 2018-02 |
| Test (`test_ec.csv`) | Swiss-Prot | 2023-01 |
| Challenging test (`price-149.csv`) | CLEAN paper, manual curation | — |

confirmed by the pipeline invocation in `data_preprocessing.sh`
(`--train_name uniprot_sprot2018_02.tsv --test_name uniprot_sprot2023_01.tsv`). Test
proteins are additionally mmseqs2-filtered against the training set at 30% and 100%
sequence-identity thresholds. The thresholds define two *training* regimes, not two test
sets: `results_test_30_go_3.csv` and `results_test_100_go_3.csv` carry the same 468 test
proteins, and `test_ec.csv` is that same 468. **We use the 100%-threshold tool results**
(`04_compare_ec_tools.sh --tool-predictions results_test_100_go_3.csv`).

**But it is not a temporal hold-out with respect to *our* training data, and that is the
distinction the methods must draw.** EC-Bench's 2018-02 → 2023-01 design is relative to
EC-Bench's own training set. Our EC CAVs are trained on the proFAB `random` positive
train split (§1.8), a 2021 snapshot with no relationship to EC-Bench's 2018-02 cut, so
EC-Bench's similarity filtering gives us no guarantee at all. Specifically:
- No script in our EC path applies any date filter. Grepping `eval_ec.py`,
  `compare_ec_tools.py`, `summarize_ec_eval.py`, `reformat_ec_goldstandard.py`,
  `ec_hierarchy_decay.py` and all five `ec_eval/0*.sh` drivers for date/temporal/cutoff
  logic returns nothing — there is no `--date-cutoff` equivalent to §2.1's.
- `test_ec.csv` carries no dates anyway: its three columns are `id`, `seq`, `ec_number`.
- What protects the EC evaluation instead is a **direct leakage audit**:
  `ec_eval/00_check_train_overlap.sh` checks every evaluated (protein, EC) pair against
  that EC's own proFAB training FASTA; **32 pairs** overlapped and were excluded via
  `--exclude ec_train_overlap_pairs.tsv` (step 2 below).
- proFAB also ships its own `temporal_` split for 294 of the 971 EC concepts, and
  `format_proFAB.py` wrote `.span` files for it, but no script reads them (§1.8).

So the accurate sentence is: *the EC test proteins come from a benchmark built as a
2018→2023 temporal split, but because our CAVs are trained on an independent 2021
resource, train/test separation is enforced here by explicit overlap removal rather than
by date.* Do not write "no temporal split" — EC-Bench has one — and do not write that our
EC evaluation is temporally held out from our training data, because it is not.

1. Gold standard: **EC-Bench** `test_ec.csv` → 509 (protein, EC) pairs, reformatted to
   long form (`protein_id`, `ec_number_raw`, `ec_number`, `ec_cav_id`) by
   `ec_eval/01_reformat_gold_standards.sh` → `specific_scripts/reformat_ec_goldstandard.py`.
   Sequences and their ESM-C embeddings are cached beside it as `test_ec.fasta` /
   `test_ec.fasta.pkl`.
   The Price-149 set was processed but **deliberately not used**: its labels disagree
   with structural matches too often. Say this, and say it was excluded a priori of the
   comparison, not after seeing results.
2. Training-overlap audit: `ec_eval/00_check_train_overlap.sh` walks every evaluated
   (protein, EC) pair and checks whether the protein appears in that EC's own training
   FASTA; overlapping pairs are written to
   `static/ec_annotations/EC-Bench/ec_train_overlap_pairs.tsv` and excluded downstream
   via `--exclude`. **32 (protein, EC) pairs** were found overlapping and removed.
3. `ec_eval/02_eval_ec.sh` → `specific_scripts/eval_ec.py`, shared scaler, `--min-n 1`
   (contrast `--min-n 5` on the GO side — the EC gold standard is far smaller, so no
   term-size floor is imposed).
   Internally this also scores proFAB's held-out `random_positive_test_max1000.span` and
   `random_negative_test_max1000.span` per EC concept (§1.8) to fit the per-concept
   positive/background Gaussians the LLR is computed against (§1.6b).
4. `ec_eval/03_summarize_ec_eval.sh` → `summarize_ec_eval.py` with the overlap exclusion.
5. `ec_eval/04_compare_ec_tools.sh` → `specific_scripts/compare_ec_tools.py`, comparator
   predictions from EC-Bench's published results table
   (`results_test_100_go_3.csv`: stacked, majority-voting-learnt/regular,
   enzbert-learnt, deepectransformer-regular, CLEAN-SupCon, …), `--llr-threshold 2.3`,
   shared scaler. Exports `ec_specificity_ranks.csv`, `ec_tool_comparison_summary.csv`,
   `ec_tool_comparison_by_level.csv`, `ec_per_term_recall.csv`, `ec_recall_vs_llr.csv`.
6. Panel A renders CAV's rank tiers alongside CLEAN-SupCon's Correct/Incorrect. Note
   the asymmetry explicitly: CLEAN-SupCon emits a single predicted EC per protein, not
   a scored candidate list, so it has no rank of the true term. n = 209 (protein, EC)
   pairs after exclusions; the tool-comparison summary covers 316 pairs at full coverage.

   **CLEAN-SupCon was not run by us — do not describe it as a local run.** No CLEAN code
   is invoked anywhere in `bioCAV` or `bioCAV-paper`. The numbers are the `clean-supcon`
   column of EC-Bench's published prediction table
   (`static/ec_annotations/EC-Bench/results_test_100_go_3.csv`), read generically by
   `compare_ec_tools.py`, which treats every non-`{id, ec_number}` column as a tool and
   normalises the name (`col.replace("-", "_")`) to give the `clean_supcon_correct`
   field that `figures.R` reads. The same table supplies the other ~13 comparators
   (stacked, majority-voting, enzbert, deepectransformer, ecrecer, blastp, …).
   Our only choices are therefore: which file (`results_test_100_go_3.csv`, the
   100%-identity training regime — `results_test_30_go_3.csv` is the alternative),
   exact-match as the correctness criterion, and the `--exclude` overlap list.

   For completeness, EC-Bench's own CLEAN configuration, from their vendored checkout
   (`EC-Bench/CLEAN/`) — attribute this to EC-Bench, not to us:
   - CLEAN's supervised-contrastive variant (`app/train-supconH.py`; the alternative in
     the same repo is `train-triplet.py`), on ESM-1b embeddings via
     `facebookresearch/esm`.
   - Inference: `CLEAN_infer_fasta.py` → `infer_maxsep(train_data, test_data,
     report_metrics=True, pretrained=False, gmm='CLEAN/data/pretrained/gmm_ensumble.pkl')`
     — i.e. **max-separation** inference (not the p-value mode), with `pretrained=False`
     so EC-Bench's own retrained weights are used rather than CLEAN's released ones.
   - Environment: Python 3.10.4, PyTorch 1.11.0, CUDA toolkit 11.3.
   - **Verify before citing:** `CLEAN_infer_fasta.py` hard-codes
     `train_data = 'train_ec_30'`, i.e. the 30%-identity training set, while we read the
     `_100_` results file. Either the script was edited per run or the naming does not
     mean what it appears to; the entrypoint's `--fasta_data init` is likewise a
     placeholder. Resolve against the EC-Bench paper before making any claim about which
     training regime CLEAN-SupCon used. `[TK]`
6b. **CAV rank = rank among the EC CAVs of the true EC's own level** (2026-09-26),
   computed by `ec_eval/rank_same_level/compute_ec_ranks_same_level.py` →
   `figure_data/ec_specificity_ranks_samelevel.csv`, column `cav_rank_level_filt`;
   replaces `cav_rank` in `ec_specificity_ranks.csv`. Steps:
   1. **Candidate set = trained EC CAVs at the same EC level: 153 for level 3, 733 for
      level 4** with the C = 10⁻⁴ library (of 953: level 1 = 6, level 2 = 61; the old
      library had 732 at level 4, of 952, lacking 3.4.22.28). Report these in the legend.
      Ancestors (e.g. 1.1 and 1 for 1.1.1) and descendants are at other levels and are
      therefore excluded by construction.
   2. **Filtered rank.** The protein's other true ECs in `test_ec_long.tsv`, truncated
      to the pair's level, are also removed from the candidates (affects 12 of 209 pairs).
   3. CAV score = `scaler_v1.pkl`-standardised embedding · unit concept vector; rank =
      `1 + #{candidates scoring strictly above the true EC}`. "No prediction" = LLR ≤ 0,
      unchanged. Pair list, `llr` and `clean_supcon_correct` are carried over unchanged
      from `ec_specificity_ranks.csv`; CLEAN-SupCon's bars do not change.
   4. **Embeddings: `static/ec_annotations/EC-Bench/test_ec.fasta.pkl`**, not the
      `--val-pkl` in `04_compare_ec_tools.sh`. The file that driver names
      (`results/ec_eval/ecbenchtest/eval_embeddings/val_proteins.span.fasta.pkl`) now
      holds only 43 proteins and could not have produced the shipped 209-pair ranks;
      `test_ec.fasta.pkl` reproduces the shipped all-952 ranks exactly (209/209 pairs).
      Fix the driver before any re-run.

   Result (Rank 1 / top 3 / top 10): level 3 (n = 161) 36.6% / 74.5% / 91.3%;
   level 4 (n = 48) 87.5% / 95.8% / 95.8%. Previously 9.3% / 31.1% / 61.5% and
   70.8% / 95.8% / 95.8%.

   **Why same-level — and state the reason correctly.** The improvement is **not**
   mainly from removing ancestors or true terms. Decomposing what outranked the true EC
   under the old all-952 pooling: for level-3 pairs, 79% of outranking CAVs were
   *unrelated ECs at a different level* (almost all level 4), 17% unrelated same-level,
   3% descendants, <1% ancestors or other true ECs. Raw projection scores are not
   calibrated across EC levels — deeper, more specific CAVs give larger projections
   (panel B shows the same climb within one lineage) — so pooling levels lets unrelated
   level-4 CAVs outrank a correct level-3 call on scale alone. For level-4 pairs, 72%
   of outranking CAVs were already same-level, 22% other-level, 5% the protein's other
   true ECs. The methods sentence is: *ranked among CAVs of the same EC level, because
   projection scores are not comparable across levels; the protein's other true ECs
   were excluded.* Note the candidate set **shrank** (952 → 153 / 732), unlike Fig 2B
   where it grew.

### 3.2 Panel B — score sharpening along the EC hierarchy
7. `hierarchy_decay/00_run.sh` → `hierarchy_decay/ec_hierarchy_decay.py`, with the
   EC-Bench validation embeddings, the long-form gold standard, all four proFAB EC
   parts and the shared scaler.
8. Protein selection (revised 2026-09-29): gold-standard proteins carrying a level-4 EC
   label whose entire 1→2→3→4 lineage has a trained CAV (55 proteins), restricted to
   those that are **leak-free at every level**, i.e. in neither the train nor the test
   span of any CAV in their lineage (20 proteins: class 1: 1, class 2: 3, class 3: 14,
   class 4: 2). Rule: one protein per level-2 EC class where possible, first accession
   alphabetically. Class 1 has only one leak-free protein.
   Shown: C0HLQ6 (1.2.1.12), P0DQP1 (2.1.1.72), P0DUB9 (3.1.1.3), P0DUE1 (3.2.2.6),
   A0A345BJN1 (4.3.1.24).
   - **Replaced:** A0A452E9Y6 (1.11.1.7, a training positive of the 1.11, 1.11.1 and
     1.11.1.7 CAVs) and P0DPS1 (2.7.7.7, a training positive of 2.7.7.7).
   - **C0HLV0 replaced by P0DUE1 (2026-10-06).** C0HLV0 (endoglucanase gh5, *Fomes meliae*,
     3.2.1.4), the class 3.2 pick under the rule, is a 27-aa UniProt fragment
     (`QAVLSETGGGDTASCDEYLYQELAYVK`, protein-level peptide evidence, PubMed 33892036): about
     1/10 of the mass of the published enzyme, with no Pfam domain. Its curve was flat
     (C = 10⁻⁴: −0.9 / 2.0 / 2.8 / 2.5 at levels 1–4), which reflects the missing sequence,
     not the hierarchy. P0DUE1 (TIR-STING, *Magallana gigas*, NAD⁺ hydrolase 3.2.2.6,
     415 aa) is the **only** other leak-free class 3.2 protein, so the replacement is forced
     by the rule once severe fragments are excluded. Its curve was computed before the swap
     was decided (3.2 / 5.6 / 10.2 / 14.4), but the choice did not depend on it.
     The other four proteins' scores are unchanged (max |Δ| 0).
   - **Open:** C0HLQ6 is also a UniProt fragment (16 aa of GAPDH, *Chattonella marina*), but
     it is the only leak-free class 1 protein and its curve sharpens (7.2 → 17.5); state
     the fragment in the legend or methods. The class 3.1 pick does not follow the written
     rule: the first leak-free 3.1 accession alphabetically is A0A345BJN2 (3.1.1.1, 218 aa;
     curve 5.9 / 7.5 / 13.8 / 21.4), not P0DUB9 (3.1.1.3). Either change the pick or the
     rule text.
   - **C = 10⁻⁴ scale:** projection scores are larger than on the old library. At level 4
     they reach 14.4–39.3 for all five (17.5–39.3 for the four kept from the old set, vs 4.1–10.4 before), 1.7–6.6×
     depending on the protein. The y axis changes accordingly; nothing is rescaled. The
     ordering of the curves also changes: P0DQP1 and P0DUB9 now end above C0HLQ6.
   - **Old data:** `hierarchy_decay/results/archive_pre_leakfree_2026-09-29/`.
9. For each protein, the CAV projection score against its level-1, -2, -3 and -4 CAVs
   is plotted as a connected line. This is a projection score, not a similarity metric.
10. The parallel GO version (`go_hierarchy_decay.py`, using `is_a` ancestor chains;
    candidates found by `find_go_chains.py`) exists and is documented in `00_run.sh` but
    is **not** in the figure.

### 3.3 Panels D and E — GO/EC cross-ontology agreement
11. `go_ec_crossval/00_run_validation.sh` → `go_ec_cosine_validation.py`, using the
    Gene Ontology Consortium's curated `ec2go` mapping
    (`go_ec_crossval/static/ec2go.txt`, snapshot checked in; `[TK]` snapshot date).
12. For each GO term with a curated EC partner and a trained CAV on both sides, compute
    the cosine similarity between the two concept vectors, and the rank of the annotated
    partner among all trained EC CAVs (953 with the C = 10⁻⁴ library; 952 old). n = 239 GO
    terms.
    - **Result, C = 10⁻⁴ (current Fig 3D/E):** matched-pair median cosine 0.831, vs a
      background median of 0.021 (sd 0.090). The partner is rank 1 for 88.3% of terms and in
      the top 5 for 98.7%, with 1 term ranked > 10.
    - **Old library:** matched median 0.601, background 0.026 (sd 0.043), rank 1 88.3%, top 5
      98.3%, 2 terms ranked > 10.
    - **Interpretation:** the swap separates matched from background pairs much more clearly
      but leaves rank 1 unchanged. The added EC 3.4.22.28 does not create a new pair: ec2go
      maps it to GO:0004197, whose curated partner remains the level-3 EC 3.4.22, still at
      rank 1.
13. Background distribution: cosine similarity of each GO CAV to all other EC CAVs,
    subsampled to `go_ec_cosine_background_sample.csv`.
14. The essential control to state: **GO and EC CAVs are trained completely
    independently** — different label sources, different positive sets, no shared
    supervision — so agreement is not built in.
15. Panel D is the two overlaid cosine histograms with medians; panel E is the rank
    histogram, binned at >20, annotated with the rank-1 percentage.
16. Supplemental table (`make_rank_vs_training_size_table.py`, `render_supp_table.py`):
    training-set size of the annotated EC partner by rank bucket, showing the rare
    rank > 10 failures are the most under-trained EC CAVs rather than biological
    mismatches.

### 3.4 Panel C — combined CAV space
17. `motif_clustering/combined_cav_umap_3d.py` (driven from `go_eval/05_visualize_umap.sh`)
    embeds all EC + GO concept vectors in one UMAP, annotating each with its category
    using `enzyme.dat`/`enzclass.txt` for EC and `go.obo` for GO namespace.
    GO namespace membership is supplied as plain include-lists built once by
    `go_eval/scripts/build_go_namespace_lists.py` (GO terms are not split into
    MF/BP/CC subdirectories on disk).
18. Coordinates exported with `--coords-csv-out` → `figure_data/go_ec_combined_umap.csv`.
    n = 7,509 concept vectors with the C = 10⁻⁴ library: EC 953, MF 1,161, BP 4,688,
    CC 707 (old library 7,507: EC 952, BP 4,687). For the new library the namespace lists
    are rebuilt from the overlay (`go_eval/results/go_namespace_lists_C1e-4/`) so that
    GO:0032436 is placed in BP.
19. Rendered as a rasterised scatter (`ggrastr`, 150 dpi) coloured by category.

---

## 4. Figure 4 — motif CAV direction space and per-residue localisation

Layout: two columns, A on the left, B over C on the right. Output 6.33 × 2.67 in,
PNG at 600 dpi. (Earlier drafts had two blank spacer panels in a single row; both were
removed 2026-09-04.)

**Layer-index convention — state this once and use it consistently.** A CAV filename
`L{k}` indexes `hidden_states[k]` directly, where k = 0 is the embedding layer and
k = 1…36 are the transformer layers. The paper indexes from 1 counting the embedding as
layer 1, so `display_layer = k + 1`. Thus `L25` is displayed as layer 26.

### 4.1 Panel A — 3-D UMAP of the motif CAV library
1. **Panel A and panel C use different motif libraries — say so explicitly.** Panel A
   plots 7,688 of the 7,692 concept vectors of the *earlier* Pfam library
   (`$AHMAD_CAV_OUTPUTS_DIR`, all layers), embedded in 3-D by UMAP. Panel C uses the
   27,472-motif library at layer 25. The asymmetry is a provenance artifact, not a
   choice: panel A's coordinates are decoded from a precomputed HTML written 2026-07-12
   by `motif_clustering/01_visualize_umap.sh`, which was archived on 2026-09-16 as stale
   once the 27k library became the standard. Either state both libraries in the legend or
   regenerate panel A on the 27k library. (This also corrects
   `distribution/DISTRIBUTION_PLAN.md`, which states that no script in the paper repo
   still references the 7,692 library — Figure 4A references it through the archived HTML.)
2. Coordinates are **not recomputed** at figure time: they are decoded from the embedded
   Plotly payload of the precomputed
   `motif_clustering/results/figures/pfam_umap_3d.html`
   (the default approximate-distance variant; the `_sim` cophenetic variant collapses
   to a near-uniform disk with no visible structure).
3. `figures/figure_data/pfam_clan_tsne/render_pfam_clan_umap_3d.py` renders a single
   fixed view, azimuth 330°, elevation 20°, coloured by Pfam clan
   (`Pfam-A.clans.tsv`), no title and no legend — clan labels are added manually in
   Illustrator afterwards. Non-highlighted points are drawn at α = 0.40, highlighted
   clans at α = 0.85.

### 4.2 Panel B — VAV_DROME layerwise localisation
4. Reproduces Figure 2 of the preceding preprint (arXiv 2511.21614v1) with freshly
   trained CAVs. Four Pfam domains: PF00621 (RhoGEF), PF00130 (C1), PF00017 (SH2),
   PF00018 (SH3).
5. CAVs trained on **all 36 layers** for those four motifs
   (`figure_data/vav_motif_repro/train_all_layers.log`; ESM++ large, batch size 4,
   `class_weight='balanced'`, shared 10k negatives as `neg_10000.jsonl`).
6. `score_vav_all_layers.py`: embed the query (Q9NHV9 / VAV_DROME) once per layer;
   slide a window of the motif's median positive-domain length at stride 1; score each
   window (scaler transform then dot with the CAV, matching
   `run_test_detection_full_embed.py::score_motif_from_hidden`); assign each residue the
   mean score over all windows containing it. Output `vav_layerwise_scores.json`.
7. Rendering: one line per layer, all 36 layers drawn, the four domain-matching curves
   highlighted; domain names are placed in a fixed right-hand column at the end of each
   highlighted line rather than in a legend, matching the reference figure.

### 4.3 Panel C — Plexin A1, a higher-domain-count example
8. `score_bigreceptors.py` — same window/score/assign procedure, **layer 25 only**,
   against the **27,472-motif library at
   `/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs`** (not the
   `paths.sh` default). Candidates MET, LDLR and Plexin A1 were scored; Plexin A1
   (8 domain types) is the one shown.
9. Window length per motif is the motif's median positive-domain length, hard-coded in
   the script (e.g. PF01403 173, PF07714 254, PF08337 249 aa).
10. B and C are assembled as a single four-row `plot_grid` (curve / track / curve / track)
    rather than two nested grids, so cowplot can align their panel edges — nested grids
    are opaque to `align="v"` and the two y-axis labels have different widths.
11. Two further worked examples (SET1_SCHPO, PLDZ_DICDI) are in the Appendix-A motif
    supplemental, `draft_fig_supp_appendixA.R`.

---

## 5. Figure 5 — SERK vs. CIK LRR-RK co-receptor specificity

Layout, three rows: [A B] / [C (0.6) | D (0.4)] / [E full width]. Output 5.5 × 7.7 in.
Working directory for all of this: `interface_search/fastas/pairs/lrr_serkcik/`.

**Panels A, B and C are not yet built.** A = structure overview (Illustrator),
B = SERK-binding interface (Illustrator), C = TBD. They render as dashed placeholder
boxes in the current PDF. Per the figure plan, B/C are intended to carry the MYC case,
which has not been built. **The MYC case shows the top 20 of the k = 128 CAV from the §7.5
sweep** (exact PCA; `interface_search/fastas/pairs/myc_mad/MYC_HUMAN.fasta_top20_pcaK128.tsv`).

### 5.1 Training and test sets
1. Positive class: Arabidopsis LRR-RKs whose co-receptor is a SERK; negative class:
   those partnering CIK. Training set = 35 Arabidopsis genes (17 SERK-partner,
   18 CIK-partner). Of the 18 negatives only 5 are true CIK founders (CLV1, RPK2,
   BAM1–3); the other 13 are hard negatives folded in during a v2 retrain (CLV2,
   PRK1–6, SRF2–8). This composition matters for interpreting §5.4 and must be stated.
2. Test set = 36 non-Arabidopsis genes (24 SERK-partner, 12 CIK-partner): 8 curated
   independently (rice, tomato, maize, canola) plus 28 pulled by homology search
   against the training genes.
3. Concept definition: a 150 aa SERK-binding interface window per gene
   (`serkcik.span`, `homolog_usable.json` coordinates).

### 5.2 Representation — the one thing that must be matched
4. **In-context pooled embeddings**: embed the full-length protein, then mean-pool the
   residues of the 150 aa window. This is what the CAV is trained and scored on.
5. An isolated-window representation (embedding the 150 aa window as a standalone
   sequence) also exists in this directory and gives materially different baselines.
   The figure uses in-context for **both** the CAV and the baseline. Say so — the
   earlier version of this panel compared a CAV on in-context embeddings against a
   cosine baseline on isolated-window embeddings, which is not a fair comparison.

### 5.3 Panel D — nearest-homolog failure cases
6. `make_hardneg_blast_heatmap.py` / `make_2class_blast_heatmap.py`: all-vs-all BLASTp
   over the 150 aa windows of the 35 Arabidopsis genes → bitscore matrix
   (`serkcik_hardneg_blast_matrix.pkl`), exported to `figure_data/fig6_heatmap.csv`.
7. The panel shows four "no close paralog" query genes (BRI1, RLP23, TDR, PSKR1)
   against their top SERK and CIK hits. Column labels are coloured by partner class.
8. Honest framing required by the analysis notes (`BASELINE_FRAMING.md`):
   - 2 of the 4 are miscalled by top-hit transfer (BRI1: BAM3 92.0 over FLS2 88.2;
     PSKR1: BAM2 89.4 over BRI1 86.7). RLP23 and TDR are correct but by only 1.2 and
     3.1 bitscore points on rows with SD 13–15.
   - **These four are selected hard cases, not a failure rate.** Across all 35 training
     genes BLAST 1-NN is correct 31/35 with a median margin of 3.63 SD, and only 3/35
     sit within ±0.25 SD of the boundary. The legend must say the panel is a hard-case
     illustration.
   - Panel D is BLAST; panel E's baseline is cosine. They are not the same nearest
     neighbour, so D does not "explain" E. Only PSKR1 fails under both; BRI1 is a
     BLAST-only failure (cosine calls it correctly). The corresponding cosine statistic
     is far worse: median margin 0.12 SD, 19/35 within ±0.25 SD of the boundary.

### 5.4 Panel E — cross-species transfer
9. `export_fig5E_data.py` regenerates `figure_data/fig6_densities.csv` and
   `fig6_auroc.csv`. Train on the 35 Arabidopsis genes, test on the 36 non-Arabidopsis
   genes; both scores computed on the same in-context embeddings.
10. Two scores, plotted as paired signed-score densities with rugs:
    - **cosine 1-NN margin** = cosine to the nearest SERK-partner training gene minus
      cosine to the nearest CIK-partner training gene. Stated as a similarity margin so
      that, as with the CAV score, positive = more SERK-like. Nearest-neighbour label
      transfer is the standard homology baseline in this literature, and a top-hit BLAST
      transfer *is* 1-NN; earlier class-mean and mean-to-positive-only framings were
      discarded as not comparable to a CAV trained on both classes.
    - **CAV score**.
11. Reported: AUROC 0.986 (cosine 1-NN) vs. 1.000 (CAV).
12. Claims this panel supports, and does not:
    - **Supported:** the CAV matches homology transfer while training only on Arabidopsis.
    - **Not supported:** "homology fails here". BLAST 1-NN margin AUROC is 1.000 on this
      split; the baselines also do not collapse on the independently curated 8-gene
      subset (cosine 0.938 / BLAST 1.000 vs. CAV 1.000). The cross-species task is
      ortholog label transfer, which is what BLAST does by construction.
    - The 0.986 is **one gene**: all four discordant pairs out of 288 involve
      `P0DL10_TD1` (maize THICK TASSEL DWARF1), margin +0.00014. Worth naming.
13. Robustness (subsampling **without replacement**, drop k of 35, 300 draws; a
    with-replacement bootstrap was run first and discarded as invalid — it removes
    candidate neighbours from a max-based 1-NN score while a fitted CAV still uses every
    retained gene, structurally penalising cosine). At k = 3: cosine median 0.986, CAV
    median 1.000; the CAV reaches perfect separation in 74% of draws vs. 0% for cosine;
    intervals overlap throughout. The defensible claim is **robustness to training-set
    composition**, not accuracy.
14. Homology-filtered evaluation (`homology_filtered_eval.py`): for each test gene, drop
    its k closest Arabidopsis training genes by BLAST bitscore (embedding-independent),
    applying the exclusion to **both** methods — the cosine neighbour pool shrinks *and*
    the CAV is refit. At k = 1 with the negative class matched to the test negatives
    (the 5 true CIK founders): cosine 0.816, CAV 0.865, paired-bootstrap gap +0.047,
    95% CI [−0.031, +0.144]. Past k = 1 both methods collapse below 0.5. Report as
    suggestive, not significant.

### 5.5 Supplemental — within-Arabidopsis leave-one-out
15. `export_fig5E_noparalog.py` → `figure_data/fig6_noparalog_scatter.csv`;
    rendered as `fig_supp_noparalog_loo.{pdf,png}`.
16. A different experiment from panel E: leave-one-out over the same 35 Arabidopsis
    genes, same in-context embeddings, CAV refit per fold
    (`class_weight='balanced'`, unit weight), so the two AUROCs are directly comparable:
    cosine 1-NN 0.876 vs. CAV 0.850 (BLAST 1-NN 0.794). At n = 35 these spreads are
    within noise — report as comparable, not as a win for either method.
17. The scatter decomposes the margin and shows all 35 genes sitting at cosine > 0.985
    to both classes; panel D's four no-close-paralog genes are labelled and all four
    land on the y = x boundary, i.e. their nearest-neighbour class calls are ties.

---

## 6. Figure 6 — single-cell case studies

Layout: A blank bar (externally supplied schematic) over three case columns B/C/D,
each a DE-vs-CAV scatter above a transcriptional-continuum strip plot. Output
7.5 × 5.82 in.

**The DE axis changed on 2026-09-29** from a per-cell donor-random-intercept mixed model
to **pseudobulk DESeq2 with a paired design** (§6.3). §6.3 below is the current
authoritative description. `FIG6_DE_METHODS.md` and `FIG6_DE_TEXT.md` describe the
superseded mixed-model version and are kept for provenance only.

### 6.1 Single-cell CAV library
1. Model: **Geneformer-V2-104M** (a different model from the protein work), with
   `token_dictionary_gc104M.pkl`.
2. Background: 10,000 cells sampled uniformly at random (without replacement,
   `np.random.default_rng(42)`) from **all human cells in the CELLxGENE Census**
   (`census_version="stable"`), with no filtering on tissue, disease, assay, dataset or
   `is_primary_data`. They are **not** drawn from the experiment dataset (item 3).
   `single_cell/00_get_negatives.sh` → `specific_scripts/build_reference_population.py`
   → `reference_population/reference_joinids.txt` (the 10,000 `soma_joinid`s),
   `reference.h5ad`, `embeddings/cells.pkl` and `global_scaler_v1.pkl` (the single-cell
   analogue of §1.3–1.4). Built 2026-07-08 (scaler 2026-07-09).
   Composition of the saved sample (from `reference.h5ad`): 798 datasets, 62
   `tissue_general` values (brain 3,103; blood 1,995; eye 803; lung 641), 162 disease
   labels (73% `normal`), 32 assays (mostly 10x 3′ v3/v2). 3,952 cells have
   `is_primary_data=False`. Only 29 cells come from dataset `b617ee1b-…`, and 6 of
   those are among the 64,783 cells downloaded for the library, so the overlap is
   negligible. The background includes blood cells even though the library excludes
   blood/PBMC/lymph node; this is expected for a Census-wide random background.
   TODO: record the resolved Census release date. "stable" changes over time and is
   not reproducible as written.
   Suggested methods text: *"The reference population comprised 10,000 cells sampled
   uniformly at random (seed 42) from the full human CELLxGENE Census (release
   <DATE>; 798 datasets, 62 tissues), with no filtering by tissue, disease, or assay."*
3. Dataset: CELLxGENE Census **dataset** (not collection) `b617ee1b-f8c8-4de9-b82b-e803ab93550d`
   (`build_cav_library_agent.py` filters on `dataset_id == '…'`).
4. Library definition: `single_cell/01_download_data.sh` →
   `specific_scripts/build_cav_library_agent.py`, an LLM-driven metadata-inspection step
   (Gemini API) that enumerates every (cell_type × tissue × disease_state) combination
   with a normal/healthy or tumour/cancer disease state, naming one CAV per combination,
   max 500 cells per group, skipping groups under 50 cells, excluding blood, PBMC and
   lymph node samples, and never pooling across tissue or disease state. The prompt is
   the `--prompt` argument in `single_cell/01_download_data.sh`, not text inside
   `build_cav_library_agent.py`. Quote it verbatim in the methods; it is the actual
   specification of the library.
5. Training: `single_cell/02_make_cavs.sh` → `specific_scripts/run_cav_pipeline.py`,
   reference n = 10,000, positives n = 500 per concept, shared `global_scaler_v1.pkl`.
   **Estimator (differs from the protein libraries in §1.6):** Geneformer 768-d embeddings →
   shared scaler → shared PCA-128 → `LogisticRegression(C = 1, class_weight = "balanced",
   solver = "lbfgs")` with sklearn's default `tol = 1e-4`, `max_iter = 1000`. The fits stop
   after 22–171 iterations but are effectively converged: refit at `tol = 1e-10`, the median
   cos to the stored CAV is 0.9995 (5th percentile 0.92). The exceptions are the excluded
   skin/melanoma CAVs (0.85–0.93) and neutrophil/breast normal (n_pos = 53, 0.93). The Fig 6
   L2 axes are unchanged (cos ≥ 0.95), so the early-stopping problem of §1.6 does not apply
   here. **Fig 6 uses the converged refit** (since 2026-10-07): same pos/neg, same shared
   scaler and PCA-128, C = 1, `tol = 1e-10`, `max_iter = 1e5` (`single_cell/convergence_check/
   01_refit_cavs.py`, arm `conv1`), hierarchy rebuilt by `02_build_coords.py` →
   `coords_conv1.tsv`, gene correlations by `03_gene_corr.sh` → `gene_corr_conv1/`, exported by
   `export_fig5_case_studies.py --coords convergence_check/coords_conv1.tsv --gene-corr-dir
   convergence_check/gene_corr_conv1` and `export_fig6_fib_modules.py --coords=...coords_conv1.tsv`.
   Against the stored CAVs: gene r Pearson 0.991 / 0.9996 / 0.9998 (neutrophil / lung10x /
   fibroblast), fibroblast L2 cell score r 0.9995; continuum genes change by one each in
   neutrophil (RPS29 → ETS1) and fibroblast (RPL39 → STMN2). Previous data:
   `figure_data/_backup_pre_fig6_conv1_2026-10-07/`. Without PCA (768-d), C = 1 at the default tol is moderately early-stopped: cos to
   its own converged fit has median 0.94 (0.71–0.97), after 18–52 vs 61–446 iterations.
   *Why PCA converges and no-PCA does not:* lbfgs stops on a gradient tolerance, and the
   distance left to the optimum is roughly gradient / smallest Hessian curvature. The
   Geneformer embedding is rank-deficient (feature-covariance eigenvalues down to 1e-14),
   and the training set is linearly separable at 768-d. Along the near-null directions only
   the penalty curves the loss: the smallest Hessian eigenvalue is 1/(C·Σw) ≈ 1e-4 at C = 1,
   so the weights there can still be far from optimal when the stop triggers. PCA-128 keeps
   only components with variance ≥ 0.9, and the smallest eigenvalue rises tenfold (≈ 1e-3).
   At C = 10⁻⁴ it is ≈ 1, so the fit converges in ~30 iterations in either frame
   (fibroblast/colorectum cancer CAV).
   The Hessian condition numbers are ~1e3 at C = 1 in both frames (and ~11 at C = 10⁻⁴). (The
   cos 0.31 between no-PCA C = 1 and the stored CAV is an estimator difference, not early
   stopping.)
5b. **Why the single-cell library keeps C = 1 rather than the library value C = 10⁻⁴**
   (comparison 2026-10-01/02, `single_cell/convergence_check/COMPARISON.md`, scripts 01–09;
   nothing promoted).
   - *What is being estimated.* Each single-cell CAV separates one cell group (e.g. cancer
     fibroblasts) from the Census-wide random background, not normal from cancer. Most of
     that difference is "fibroblast vs everything else", which the normal-fibroblast CAV
     shares. The L2 condition axis (item 7) is the small residual left after removing the
     shared cell-type (L0) and tissue (L1) directions.
   - *What C decides.* C sets which part of the group-vs-background difference the CAV keeps.
     At C = 1 the fit is weakly regularised, so it uses the covariance structure. It
     down-weights directions that vary a lot in every cell type and keeps the specific,
     low-variance directions that separate the groups. At C = 10⁻⁴, with C·n ≈ 1 (n = 10,500),
     the penalty is as large as the loss and the fit is pulled toward the raw class-mean
     difference. That difference is dominated by the embedding's highest-variance
     directions, i.e. coarse cell-identity axes such as immune vs non-immune.
     Fibroblast/colorectum, PCA-128 frame (`convergence_check/why_C.tsv`); variances are
     measured on the 10k background and expressed relative to the mean eigenvalue:

     | | current (C = 1) | C = 10⁻⁴ |
     |---|---|---|
     | cos(cancer CAV, raw class-mean difference) | 0.49 | 0.86 |
     | background variance along the cancer CAV | 1.3 | 2.3 |
     | background variance along the L2 axis | 0.64 | 1.71 |

   - *Consequence for L2.* At C = 10⁻⁴ the normal and cancer CAVs collapse onto the same
     fibroblast-vs-Census direction: cos(normal CAV, cancer CAV) goes from 0.79 to 0.98, and
     the L2 residual norm before rescaling falls from 0.60 to 0.22. The rescaled residual is
     then close to the mean shift between cancer and normal positives (cos 0.78, vs 0.13 at
     C = 1). The same holds for epithelial/lung10x (0.55 vs 0.15) and neutrophil/breast
     (0.72 vs 0.26).
   - *What that mean shift is.* In this atlas it is largely immune contamination of the
     fibroblast annotation. Among fibroblast/colorectum cells, ≥ 2 of PTPRC/CD3E/CD52 are
     detected in 27% of cancer vs 16% of normal cells, and MS4A1 in 5.7% vs 0.9%. Under
     C = 10⁻⁴ the fibroblast L2 axis correlates with PTPRC (r 0.12 → 0.42) and CD52
     (0.15 → 0.49), and against COL1A1 (−0.19 → −0.31).
   - *Gene modules* (clustered co-expression of the top/bottom 20 genes,
     `convergence_check/08_fib_module_heatmap.py`, `figures/previews/fib_coexpr_top20_*`):
     - C = 1: four modules, namely T cell, plasma, a 20-gene ECM/quiescent-fibroblast block
       and an ADAMDEC1/CXCL14/CCL13/BMP5/ADAM28 mucosal-fibroblast block.
     - C = 10⁻⁴ (with or without PCA): two blocks, lymphocytes vs generic non-lymphocyte
       genes. ADAMDEC1 r goes 0.21 → −0.03.
     - Top/bottom 20 by pseudobulk DE: small non-fibroblast modules (myeloid, epithelial,
       endothelial, plasma, NK) plus one SFRP1/FBLN1 fibroblast module, mostly with CAV r ≈ 0.
   - *PCA is not the driver.* At C = 10⁻⁴, the no-PCA and PCA-128 CAVs agree (raw cos 0.89,
     L2 cos 0.73–0.84) and give near-identical gene lists. At C = 10⁻⁴ the outcome is set by C.
   - *Why the protein libraries differ.* There each CAV is used directly as concept vs
     background, the class mean is a good concept direction, and nothing is subtracted
     afterwards. The single-cell analysis instead takes the difference of two nearly identical
     CAVs, which needs the specific part that strong shrinkage discards.
   - *Caveat.* This is a mechanism, not a proof that C = 1 is optimal. The support for C = 1
     is empirical: the module structure above, and the ECM program whose CAV r is −0.30 to
     −0.37 while label r ≈ 0.
   - Higher DE concordance under C = 10⁻⁴ (fibroblast ρ 0.23 → 0.39, §6.4 item 21) is
     expected for a near-mean-difference axis and is not evidence of a better axis.
6. The resulting library contains **93 single-cell CAVs**. Library structure inferred by
   `specific_scripts/analyze_cav_library.py`
   (`03b_get_lib_structure.sh`) → `library_structure.json`.
7. Hierarchical decomposition: `single_cell/04_build_hierarchy.sh` →
   `specific_scripts/cav_hierarchy.py`, producing per-cell coordinates on three nested
   axes. Two level orders were built: `group-context-condition`
   (L0 = cell type, L1 = tissue) and `context-group-condition` (L0 = tissue,
   L1 = cell type). In both, **L2 is the condition axis with tissue and cell type
   projected away**, and L2 spans the same subspace either way — only L0/L1 differ in
   which variable is removed first. All figure statistics use L2.
8. Gene correlates: `single_cell/05_get_disease_genelists.sh` →
   `specific_scripts/cav_gene_correlation.py`, correlating each cell's L2 score with
   per-gene expression after library-size normalisation and `log1p`. State the
   normalisation: without it a gene's correlation can be dominated by per-cell library
   size rather than its own biology (worked example in the script header —
   MALAT1/fibroblast/melanoma, raw-count r = +0.46 vs. normalised r = −0.62, the latter
   agreeing with the DE direction). Minimum 10 cells expressing a gene.

### 6.2 Case selection and the two confound fixes
9. Three cases shown, one per column:

   | Panel | Pair | Paired-donor cells (normal / cancer) | Paired donors | Used in pseudobulk DE (≥ 10 cells/sample) |
   |---|---|---|---|---|
   | B | neutrophil, breast | 66 / 231 | 3 (all donors paired) | 2 donors, 46 / 230 cells |
   | C | epithelial cell, lung — **10x 3′ v2 cells only** | 160 / 498 | 4 (of 8) | 4 donors, 160 / 498 cells |
   | D | fibroblast, colorectum | 1,000 / 1,000 | 7 (all donors paired) | 6 donors, 997 / 980 cells |

   The CAV gene correlations use every paired-donor cell (column 3). Panel D's
   1,000 / 1,000 is the atlas's build-time cap of 1,000 cells per condition per group.
   (`sc_pair_meta.csv` records 514 / 486 for D. That is the continuum strip's
   **plotting-only** subsample of 1,000 of the 2,000 cells, `max_cells=1000` in
   `cav_continuum_viz.get_pair_data`, not the analysed cell count. An earlier version of
   this note got that backwards.)
10. **Assay confound.** Skin/melanoma pairs were dropped as unfixably confounded
    (normal = 10x 3′ v2, melanoma = Smart-seq2, zero overlap). The original lung pairs
    carried a milder version (normal ~91% 10x, cancer ~50% 10x) and were replaced by a
    from-scratch CAV retrain restricted to `assay == "10x 3' v2"`
    (`scripts/train_10x_lung_cavs.py`, `extend_hierarchy_lung10x.py`,
    `lung10x_de_and_corr.py`). "lung10x" is a naming convention for that
    assay-restricted axis, not a `tissue` value in the atlas.
11. **Paired-donor restriction.** Both axes are computed on cells from donors
    contributing *both* conditions (`scripts/paired_donor_gene_corr.py` for the CAV
    correlations, `scripts/paired_donor_pseudobulk.py` for DE). Donor-blind pooling is
    vulnerable to a Simpson's-paradox confound when single-condition donors are mixed
    in. Panel C went from 1,407 cells / 8 donors to 658 cells / 4 donors. **One
    asymmetry to state:** the DE additionally drops pseudobulk samples of < 10 cells
    (item 13), which removes one donor each from B (AziziBRCA_1: 20 normal + 1 cancer
    cell) and D (QianCRC_CRC_6: 3 normal + 20 cancer cells). The CAV correlations and
    continuum strips keep those donors. So the two scatter axes share donors and cells
    except for these 21 (B) and 23 (D) cells; C is identical on both axes.

### 6.3 Differential expression axis (x) — pseudobulk DESeq2, paired design
Standard multi-donor single-cell DE protocol: pseudobulk aggregation, then a bulk
RNA-seq count GLM with donor as the unit of replication (Squair et al. 2021 *Nat
Commun* 12:5692; Crowell et al. 2020 *Nat Commun* 11:6077 / muscat). It was chosen as a
standard test that uses the donor as the unit of replication. That makes it more
conservative than any per-cell test, the earlier mixed model included: those treat cells
as the unit of evidence and call far more genes (lung: 1,867 genes at padj < 0.05 under
the mixed model vs 145 here). Between the two standard pseudobulk engines, neither is
uniformly stricter (item 16). The plotted statistic, Spearman ρ over log2FC, is
engine-independent (ρ ≥ 0.96 between engines); the engine choice only changes which
points are coloured as DE.

12. **Aggregation.** `scripts/paired_donor_pseudobulk.py`: raw integer UMI counts
    (`cells.h5ad` `.X`, verified integer) are **summed** per donor × condition over the
    paired-donor cells from §6.2. The result is one normal and one cancer pseudobulk
    sample per donor, written to `results/de_pseudobulk_paired_donors/pseudobulk_counts/`
    with each sample's `n_cells`.
13. **Sample QC.** Pseudobulk samples built from < 10 cells are dropped (muscat's default
    minimum). A donor left with only one condition is then dropped too, because the
    design is paired. Effect: B keeps 2 of 3 donors (4 samples), C 4 of 4 (8 samples),
    D 6 of 7 (12 samples). Smallest retained samples: B 10 normal / 78 cancer cells,
    C 25 / 71, D 49 / 57.
    Sensitivity (`min1`, no minimum): B then keeps a 1-cell, 609-count sample. Because
    `filterByExpr` scales its CPM cutoff by the median library size, B's tested genes
    shrink from 1,326 to 350 (D: 5,429 → 5,084; C unchanged). The Spearman ρ with CAV r
    is essentially unchanged (B 0.79, C 0.80, D 0.22 shrunk-LFC). Library sizes in B span
    ~700-fold, outside what `filterByExpr` assumes.
14. **Gene filter.** `edgeR::filterByExpr(counts, group = condition)` with defaults
    (min.count 10, min.total.count 15, min.prop 0.7). The same filtered matrix goes to
    both engines. Genes tested: **1,326 / 1,446 / 5,429** for B / C / D. Genes failing the
    filter are **untested, not DE-null**; any "CAV finds genes DE misses" statement is
    restricted to tested genes.
15. **Model (primary, plotted).** `DESeq2` 1.46.0, `design = ~ donor + condition`, normal
    as the reference level, default size factors, parametric dispersion fit (local fit
    substituted automatically where the parametric trend fails), Wald test on
    `condition_cancer_vs_normal`, `results(alpha = 0.05)` with default independent
    filtering and Cook's-distance handling, BH-adjusted within each pair. Residual df =
    samples − donors − 1: **1 / 3 / 5** for B / C / D. The x axis is DESeq2's
    **unshrunk (MLE) log2FoldChange**, a true log2 ratio of normalised means (column
    `log2fc_mle`; the TSV's `log2fc` column is the apeglm-shrunk value and is not
    plotted).
    `lfcShrink(type = "apeglm")` (apeglm 1.28.0) is saved for robustness. It is not
    plotted because with 1–5 residual df it collapses almost every gene to ~0, which
    shows DE as less informative than it is.
16. **Robustness engine.** `edgeR` 4.4.0 quasi-likelihood on the same filtered counts and
    design: TMM normalisation, `estimateDisp`, `glmQLFit(robust = TRUE)`,
    `glmQLFTest(coef = "conditioncancer")`. LFC agreement with DESeq2 (Spearman) is
    0.957 / 0.999 / 0.995. The **number** of DE genes is engine-dependent (DESeq2 vs edgeR
    at padj < 0.05, |LFC| ≥ 1: B 30 vs 97, C 143 vs 35, D 2 vs 0), so the Results should
    not lean on a DE-gene count. R 4.4.3, Bioconductor 3.20; env
    `/xdisk/clairemcwhite/clairemcwhite/envs/pseudobulk_de`; `sessionInfo()` in
    `results/de_pseudobulk_paired_donors/min10/sessionInfo.txt`.
17. **DE call used for point styling:** padj < 0.05 and |log2FC| ≥ 1. B has NA padj for
    720 of its 1,326 genes (DESeq2 independent filtering at 1 residual df). Those genes
    were tested, stay in the Spearman ρ, and count as not DE. DE genes: B 30 (all up in
    cancer), C 143 (43 up / 100 down), D 2 (IGLC3 up, FCER1G down).
18. **Why not the old DESeq2 code path.** `bioCAV/specific_scripts/de_methods.py::
    run_deseq2_pair` aggregates **per donor only** and labels each aggregate with its
    first cell's condition. With paired donors this merges each donor's normal and
    cancer cells into one sample. The "pseudobulk fails on neutrophil/breast" statement
    in the superseded docs was that bug, not a power limit, and the tables in
    `results/de_vs_cav/` are invalid. The shared function is **not yet fixed** but no
    longer feeds any figure.
19. Superseded (kept for provenance, not used by Fig 6): per-cell `statsmodels` mixed
    model with a donor random intercept on log2(CP10K+1)
    (`scripts/paired_donor_de_mixedlm.py` → `results/de_mixedlm_paired_donors/`), whose
    "log2fc" was a difference in mean log expression, not a log ratio. On the same genes
    its agreement with CAV r was slightly *higher* (ρ 0.81 / 0.83 / 0.27) because it is
    nearly the same per-cell statistic as a label correlation (ρ 0.94–0.95 with per-cell
    r(expression, condition label)). The pseudobulk comparison is the harder baseline.
    **The supplements still use it**: `fig_supp_projection_comparison` (§7.1) and the
    subpopulation supplement read `de_mixedlm_paired_donors`.

Suggested methods text: *"Differential expression was tested by pseudobulk analysis with
a paired design. For each case-study pair, raw UMI counts were summed per donor and
condition over cells from donors contributing both conditions. Pseudobulk samples
comprising fewer than 10 cells were removed, together with the donor's other sample,
leaving 2, 4 and 6 donors for the neutrophil/breast, epithelial/lung and
fibroblast/colorectum pairs. (CAV gene correlations used all paired donors, 3, 4 and 7
respectively, so the two analyses differ by 21 and 23 cells in the neutrophil/breast and
fibroblast/colorectum pairs.) Genes were filtered with edgeR filterByExpr (1,326, 1,446 and
5,429 genes retained), and differential expression was tested with DESeq2 (v1.46) using
the design ~ donor + condition, with normal tissue as the reference level (Wald test,
Benjamini–Hochberg correction within each pair). We report DESeq2's maximum-likelihood
log2 fold change; genes with adjusted P < 0.05 and |log2FC| ≥ 1 were considered
differentially expressed. Results were robust to the choice of engine (edgeR
quasi-likelihood; Spearman correlation of log2FC with DESeq2 ≥ 0.96), to apeglm shrinkage
of fold changes, and to removing the minimum-cell threshold."*

### 6.4 Export and rendering
20. `single_cell/08_export_fig5_data.sh` → `scripts/export_fig5_case_studies.py` writes
    `figure_data/sc_pair_meta.csv`, `sc_de_vs_cav.csv`, `sc_continuum_cells.csv`,
    `sc_continuum_genes.csv`. Cell/gene selection reuses `cav_continuum_viz.py`'s logic.
    The DE step must be run first: `scripts/paired_donor_pseudobulk.py`, then
    `Rscript scripts/paired_donor_pseudobulk_de.R 10` in the pseudobulk_de env (both from
    `single_cell/`). `sc_de_vs_cav.csv` carries `log2fc` (= DESeq2 MLE), `padj`,
    `n_donors` and `cav_r`.
21. Scatter (upper): every **tested** gene (DE-tested ∩ CAV-correlated) is plotted. The
    statistic is a **Spearman correlation over all of them**. Since 2026-10-05 it is
    reported here, not in a panel subtitle. On the converged C = 1 CAVs (conv1, the live Fig 6 since
    2026-10-07): B ρ = 0.76 (n = 1,326, 2 donors, p = 3e-248), C ρ = 0.81
    (n = 1,446, 4 donors, p < 1e-300), D ρ = 0.22 (n = 5,429, 6 donors, p = 2e-62). (Stored,
    unconverged CAVs gave 0.75 / 0.81 / 0.23.)
    These are **not comparable** to the superseded mixed-model values (0.60 / 0.73 / 0.27),
    which were computed over 10–13k genes. On the same gene set the mixed model gives
    0.81 / 0.83 / 0.27. **No DE call is drawn (since 2026-10-05):** every tested gene is one
    grey cloud and only each case's continuum genes (top/bottom 5 by CAV r) are highlighted
    and labelled.
    - Before 2026-10-05, DE genes (padj < 0.05, |log2FC| ≥ 1) were shown in blue. That
      colour clashed with the blue "normal" condition colour, and the text never used the DE
      call; the effect size can be read from the x axis.
    - Earlier still, the top 8 DE genes were also labelled.
    - DE gene counts at that threshold, for the text if needed: B 30, C 143, D 2. At
      padj < 0.01 and |log2FC| ≥ 2: 14 / 64 / 1. For reference, on the library-consistent single-cell CAVs
    (shared scaler, no PCA, C = 10⁻⁴, converged; `single_cell/convergence_check/COMPARISON.md`)
    the same Spearman over all tested genes is B ρ = 0.76 (n = 1,326, 2 donors, p = 4e-249),
    C ρ = 0.81 (n = 1,446, 4 donors, p < 1e-300), D ρ = 0.39 (n = 5,429, 6 donors, p = 6e-198).
    The gene sets and donor counts are identical to the current figure. The rise in D is
    expected by construction, because the C = 10⁻⁴ CAVs lie close to the class-mean
    difference (cos ≈ 0.86); it is not independent evidence of a better axis. The
    per-panel tables are `figures/previews/fig6_C1e-4_{stored,np_conv4}_stats.tsv`.
22. Continuum strip (lower half of **B and C only** since 2026-10-05; D's lower half is the
    fibroblast-program panel, item 22b): cells ordered along the L2 condition axis, drawn as
    **condition-split strips** (normal and cancer as separate rows), not as interleaved
    jitter, with the top continuum genes' scaled expression shown as grey-ramp rows
    beneath. Subsampled to 1,000 cells per pair **for plotting only**; all statistics use
    the full available cells. State in the Results, not hedged in Methods, that
    fibroblast/colorectum has only 2 DE genes (IGLC3, FCER1G) of 5,429 tested. The
    quiescent-fibroblast continuum genes (MGP, DCN, OGN, C3, CCDC80) and ADAMDEC1 all have
    padj ≈ 1 (MLE |log2FC| ≤ 0.46). Normal and cancer cells overlap substantially on the L2
    axis in that case.
22b. **Fig 6D lower half: fibroblast program panel (LIVE in `figures.R` since 2026-10-05; it
    replaced D's continuum strips).**
    - Data: `single_cell/scripts/export_fig6_fib_modules.py` (from `single_cell/`) →
      `figure_data/sc_fib_module_{scores.csv,stats.tsv}`.
    - Rendering: `sc_build_fib_modules()` in the Fig 6 block of `figures.R`.
    - Two stacked plots, ECM/quiescent above ADAMDEC1+, show module score vs L2 as per-condition
      means in L2 deciles with 95% CI ribbons. Each title gives the module's gene list.
    - The x axis reads "L2 score (cancer-CAV direction)".
    - The original preview driver is `figures/draft_fig6_fib_modules.R`, and the staging copy of
      the data is `figure_data_staging/singlecell_modules/` (identical values).
    - *Programs* are the named blocks of the clustered co-expression heatmap of the top/bottom 20
      genes by CAV r, using the converged C = 1 CAVs (`convergence_check/08_fib_module_heatmap.py conv1`;
      **supplement `fig_supp_fib_heatmaps` (LIVE in `figures.R` since 2026-10-07)**: top row A = conv1,
      B = top/bottom 20 by pseudobulk Wald statistic (labels gene (log2FC; CAV r from A)); bottom row
      C = C = 10⁻⁴ PCA-128 converged (conv4), so A and C differ only in C. Gene labels 6.2 pt; colour =
      sign of the ranking statistic. Preview driver was `figures/draft_supp_fib_heatmaps.R`):
      - ECM / quiescent fibroblast: the **9 annotated ECM genes** MGP, DCN, OGN, C3, CCDC80, MFAP5,
        SFRP2, DPT, LUM (all inside the heatmap's 20-gene ECM block, which also contains PSAP,
        GPNMB, FBN1, IGFBP6, PRELP, GSN, SERPING1, EFEMP1, C1R, CST3, TIMP1). The 20-gene version
        gives r −0.43 / −0.36 / −0.50 (all / normal / cancer).
      - ADAMDEC1+ mucosal fibroblast, 7 genes: BMP5, ADAM28, CXCL14, CCL13, ADAMDEC1, STMN2, HAPLN1
        (HAPLN1 enters the top 20 under the converged fit; the stored-CAV block had 6 genes).
      - NSG1 (ADAMDEC1 block) and RPL39 (plasma block) are excluded by the membership rule: a gene is dropped when its mean
        co-expression with the module is < 0.5 × the module median.
    - *Score* = mean over the module genes of z-scored log1p(CP10K), across the 2,000 paired-donor
      fibroblasts. The panel plots per-condition means in L2 deciles with 95% CI ribbons.
    - *Statistics* (Pearson r of module score with the L2 score):

      | | all cells | normal (n = 1,000) | cancer (n = 1,000) | per donor (7 donors) |
      |---|---|---|---|---|
      | ECM / quiescent (9 genes) | −0.41 | −0.35 | −0.47 | −0.32 to −0.46, all p < 0.05 |
      | ADAMDEC1+ | +0.26 | +0.32 | +0.20 | +0.24 to +0.42 in 3 donors (CRC_1, 4, 7); ≈ 0 (−0.07 to +0.08, n.s.) in 4 |

      The mean score barely differs between conditions (normal − cancer: ECM +0.05 SD,
      ADAMDEC1 +0.03 SD). The programs therefore track position along the L2 axis, not the
      condition label. In the top L2 decile the ADAMDEC1 program is higher in normal than in
      cancer cells.
    - *Caveat for the text:* the module genes were selected by abs(CAV r), so the
      score-vs-L2 correlation is partly by construction. The non-circular findings are:
      (i) the selected genes form coherent co-expression modules;
      (ii) the trend holds within each condition and, for ECM, in every donor;
      (iii) the programs are DE-null (pseudobulk padj ≈ 1, §6.4 item 22).
    - Wording to avoid: "loss of an ECM program in cancer" or "increased ADAMDEC1 in cancer".
      Use instead: "along the L2 axis, cells lose a coordinated ECM/quiescent program and, in a
      subset of donors, gain an ADAMDEC1+/CXCL14+ mucosal-fibroblast program."
22c. **How much does the fibroblast/colorectum L2 axis separate the conditions?** (2026-10-05,
    `single_cell/convergence_check/12_fib_axis_separation.py conv1`; paired-donor cells, converged C = 1 CAVs;
    the stored CAVs gave the same numbers to within 0.02: held-out 0.562, shift 0.31 SD, L2 means 1.24 / 1.52, per-donor max 0.75).
    - *Construction.* The "cancer" end is defined by construction. L2 is the cancer-fibroblast CAV
      (500 cancer fibroblasts vs the Census background), minus its projection on L0 (mean
      fibroblast baseline) and L1 (colorectal normal-fibroblast residual). The arrow points
      along the remaining part of the cancer CAV; it does not come from separating normal from
      cancer fibroblasts.
    - *Separation is weak:*
      - AUROC 0.585 on all 2,000 cells. 1,000 of these were CAV training positives (500 per
        condition); on the other 1,000 the AUROC is 0.561.
      - The mean shift is 0.30 SD (L2 normal 1.16 ± 0.89, cancer 1.44 ± 0.91).
      - Per donor, AUROC is 0.51–0.73: 0.51–0.60 in 4 donors, 0.68–0.73 in 3.
    - *Null comparison.* Random unit directions orthogonalised against the same L0/L1 basis
      (taking the better of the two orientations) give a median AUROC of 0.525, 95th percentile
      0.565 and 99th percentile 0.584. The observed 0.585 exceeds 99.25% of them. The held-out
      0.561 is at the null's 95th percentile.
    - *Ceiling.* A logistic regression trained directly on normal vs cancer fibroblast embeddings
      (PCA-128 frame, leave-one-donor-out) reaches only AUROC 0.522. In this dataset the Geneformer
      embeddings of colorectal fibroblasts carry almost no donor-generalisable normal/cancer
      signal.
    - *Consequence for Fig 6D.* The fibroblast L2 axis is a direction of within-population
      heterogeneity (ECM/quiescent vs ADAMDEC1+/immune-associated programs) that is present in
      both conditions, with a slight cancer-ward tilt. It is not a normal→cancer separator. The
      "← normal  cancer →" axis label overstates this for D. For comparison, the in-sample L2
      AUROCs are 0.87 (neutrophil/breast) and 0.84 (epithelial/lung10x), §6.4 item 22b and
      `convergence_check/axis_metrics.tsv`.
22d. **Decision 2026-10-05: option 2. Fibroblast/colorectum stays as the third case, reframed as
    the deliberate contrast case.** It shows no disease separation and no DE (items 22c and 22e),
    yet the CAV axis exposes coherent fibroblast programs that DE cannot see.
    - The case-selection sheet is `figures/previews/fig6_case_options.*`, item 22e.
    - D's lower half is no longer the continuum strips. D's x axis reads "L2 score (cancer-CAV
      direction)", not "← normal  cancer →".
    - Full-figure mocks: `figures/draft_fig6_optionD.R [lines|blocks]` →
      `previews/fig6_optionD_{lines,blocks}.*`.
      - *lines*: module score vs L2 (item 22b).
      - *blocks*: a shrunken co-expression heatmap of the 27 module genes
        (`draft_fig6_fib_module_summary.R`, option C).
    - **Promoted 2026-10-05 (lines variant)**, item 22b. Pre-change backup:
      `figures/figure_data/_backup_pre_fig6D_modules_2026-10-05/` (figures.R, cav_fig6.{pdf,png}).
22e. **Case-selection sheet, all 10 usable pairs** (`single_cell/convergence_check/14_allpair_separation.py`).
    - DE for the 7 non-figure pairs comes from copies of the pseudobulk pipeline run into
      `convergence_check/de_pseudobulk_allpairs/` (13a/13b). They reproduce the live min10
      results for the existing pairs.
    - Held-out L2 AUROC vs the 95th percentile of random directions (n DE genes):

      | pair | held-out L2 AUROC | random 95th pct | DE genes |
      |---|---|---|---|
      | neutrophil / breast | 0.87 | 0.75 | 30 |
      | epithelial / lung10x | in-sample 0.84 only (all normal cells are CAV positives) | 0.75 | 143 |
      | neutrophil / lung | 0.67 | 0.81 | 133 |
      | B cell / colorectum | 0.65 | 0.63 | 4 |
      | T cell / colorectum | 0.64 | 0.59 | 7 |
      | epithelial / colorectum | 0.58 | 0.58 | 1 |
      | fibroblast / colorectum | 0.56 | 0.57 | 2 |
      | mast / colorectum | 0.55 | 0.58 | 0 |
      | endothelial / colorectum | 0.53 | 0.55 | 2 |
      | mononuclear phagocyte / colorectum | 0.53 | 0.57 | 2 |

    - Every colorectal pair shows weak separation and almost no DE.
23. `single_cell/09_export_fig5_umap_data.sh` → `export_fig5_umap_data.py` exports an
    overview t-SNE (raw CAVs plus hierarchy levels L0/L1/L2, max 20,000 cells, restricted
    by `--clean-pairs-csv` to cell_type × tissue combinations where both conditions are
    ≥ 80% the same dominant assay, and by `--paired-donors-only`). **These panels were
    removed from the figure** — they showed near-complete normal/cancer overlap at every
    level and carried nothing the case columns do not carry better. Panel A's bar is now
    blank for an externally supplied schematic. The CSV
    (`sc_umap_overview.csv`, 63,928 rows) is still exported; decide whether it becomes a
    supplemental or is dropped.

---

## 7. Supplemental figures

### 7.1 Orthogonalised-projection comparison (`fig_supp_projection_comparison`)
1. `single_cell/10_export_projection_comparison_supp_data.sh` →
   `scripts/export_projection_comparison_supp_data.py` → `figure_data/proj_de_vs_cav.csv`.
2. Compares the DE-vs-CAV relationship across three stages of the condition axis: raw
   condition CAV, L0 removed, and L2 (tissue + cell type removed). Three stages, not
   four — the script header explains why.
3. Same Spearman-over-all-genes fix as §6.4 item 21; the earlier subset-based version had
   the identical liability. **DE axis here is still the superseded mixed model**
   (`de_mixedlm_paired_donors`), not Fig 6's pseudobulk DESeq2; see §6.3 item 19 and the
   open item below.
4. **CAVs: converged C = 1 PCA-128 since 2026-10-07**, matching Fig 6 (`--cav-dir
   convergence_check/cavs_conv1 --coords convergence_check/coords_conv1.tsv`; the L2 stage
   reproduces `coords_conv1.tsv` to 2e-15). Spearman ρ (DE vs CAV), stored → converged:
   neutrophil raw 0.05 → 0.07, L0-removed 0.56 → 0.59, L2 0.60 → 0.61; lung10x 0.66 → 0.67,
   0.75 → 0.75, 0.73 → 0.73; fibroblast 0.16 → 0.17, 0.07 → 0.07, 0.27 → 0.27. Previous data in
   `figure_data/_backup_pre_fig6_conv1_2026-10-07/`.

### 7.2 Subpopulation-abundance shift (`fig_supp_fibroblast_subpop`) — DEMOTED 2026-10-07, superseded by `fig_supp_fib_heatmaps` (§6.4 item 22b)
4. `single_cell/scripts/export_subpop_supp_data.py` → `subpop_gene_scatter.csv`,
   `subpop_gene_scatter_background.csv`, `subpop_gene_corr.csv`, `subpop_cell_scores.csv`.
5. Worked example: the quiescent-fibroblast module in fibroblast/colorectum (MGP, DCN,
   OGN, C3, CCDC80) is DE-null (padj > 0.5 for all five under the mixed model this
   supplement still uses; padj ≈ 1 for all five under Fig 6's pseudobulk DESeq2, so the
   claim holds under the stricter test too) but a strong CAV correlate
   (r = −0.33 to −0.37, p ≈ 1e-53 to 1e-65) — consistent with a subpopulation whose
   relative abundance shifts with disease rather than uniform per-cell downregulation.
   A group-mean DE test is powered to detect the latter and not the former.

### 7.3 SERK/CIK leave-one-out — see §5.5.

### 7.4 EC evaluation panel (`fig_ec_eval.pdf`)
6. Built from `ec_per_term_recall.csv` (per-EC-term recall distribution, one curve per
   tool, restricted to tools with full overlapping coverage with the CAV EC set),
   `ec_tool_comparison_summary.csv` and `ec_recall_vs_llr.csv`. Confirm whether this
   ships as a supplemental or was superseded by Figure 3A.

### 7.5 Number of PCA components as a training option (`fig_supp_pca_sweep`)
Working directory: `interface_search/fastas/pairs/myc_mad/`. Two single-positive concepts:
**MYC 366–420** (P01106, 55 aa; the bHLH domain alone, UniProt 369–421, excluding the leucine
zipper at 428–449) and **MAD1 61–94** (Q05195, 34 aa, inside the MAD1 bHLH 56–108). Each has one
span-pooled positive plus the §1.3 10k negatives in the §1.4 shared scaler frame
(`MYC_HUMAN.fasta_concept/`, `MAD1_HUMAN.fasta_concept/`).

7. **Sweep CAVs:** `pca_sweep_mad1/build_pca_sweep.py CONCEPT_DIR OUT_PREFIX` (CPU, seconds).
   Recipe:
   - standardize with the shared scaler;
   - fit exact PCA (`svd_solver='full'`) on the 10k standardized negatives, keeping the top
     k ∈ {8, 16, 32, 64, 128, 256, 512} components;
   - fit balanced L2 logistic regression (C = 1, lbfgs, seed 42, sklearn default `tol = 1e-4`,
     i.e. the **old** early-stopped estimator, not the current §1.6 default) in the
     k-dimensional PC space;
   - fold the coefficient back to 1152-d (`components_.T @ coef`) and unit-normalize.

   k = 1152 skips PCA and reproduces the default concept vector (cosine 1.000).
   - The PCA basis is fit on the shared negatives only, so it is identical for every concept:
     the frame stays shared, and scoring is unchanged (standardize, then dot product).
   - Validation: the script reproduces the original MYC sweep (`pca_sweep/myc_pca_sweep_cavs.npz`)
     at cosine 1.000 for all k.
8. **Proteome scan:** `pca_sweep*/score_pca_sweep.py` via `run_pca_sweep.sbatch`, one job per
   concept (300 GB, about 5 min; the 94 GB `uniprot_human_all.fasta.pkl` does not fit in 64 GB).
   - Scoring is the same as `query_proteins_by_cav.py --sliding-window`: per-residue embeddings
     are mean-pooled over a window of the concept's span length (55 / 34) at stride 1, and the
     protein score is the maximum over its windows. It is vectorized with cumulative sums, and
     all k are scored in one pass.
   - Validation: the MYC k = 1152 column and the folded May PCA-128 CAV reproduce the production
     rankings `MYC_HUMAN.fasta_query_results_sharedscaler.tsv` and `..._mayPCA128_folded.tsv`
     (Spearman 1.000).
   - 19,997 reviewed human proteins.
9. **Family labels:**
   - bHLH = the PF00010 list in `human_hlh_bzip_accessions.json` (InterPro API, reviewed
     human, 2026-09-28) plus UniProt reviewed human entries with IPR011598, IPR036638 or
     PS50888 (`human_bhlh_interpro_IPR011598_IPR036638_PS50888.tsv`, 115 accessions in the
     union).
   - PF00010 alone misses the bHLH-PAS proteins (SIM1/2, NPAS1/3, HIF1A/3A) and must not be
     used by itself.
   - bZIP = PF00170 / PF07716 (41). No protein is in both sets.
10. **Export and render — PROMOTED 2026-10-05: C = 10⁻⁴ converged sweep, judged by CACGTG.**
    - CAVs: `converged_C1e-4/build_sweep_converged.py` (step 7 recipe with balanced LR C = 10⁻⁴,
      `tol = 1e-10`, `max_iter = 1e5`), k = 8, 16, 32, 64, 128, 256, 512, 1152 plus the fine grid
      80–176 (`*_fine_cavs.npz`); scored with the step 8 scorer (`converged_C1e-4/score_*.py`).
    - Labels: CACGTG (canonical E-box) binder / other motif / unlabelled per protein,
      `converged_C1e-4/ebox_label/` (`label_cacgtg.py`: HOCOMOCO v12 core + JASPAR 2024 vertebrate,
      rule fixed before ranking; `augment_labels.py` adds `literature_overrides.tsv` and inference
      for proteins with no database motif → label sets db / db+lit / db+lit+inf;
      `eval_augmented.py` → `cacgtg_eval_augmented.tsv`).
    - Export: `pca_sweep_supp/export_pca_sweep_cacgtg_data.py` → `figure_data/pca_sweep_cacgtg.csv`
      → `figures.R` section "Supplemental — number of PCA components" → `fig_supp_pca_sweep.{pdf,png}`
      (7 × 2.6 in; identical to the approved preview `previews/supp_pca_sweep_C1e-4.png`).
    - Panels (label set db+lit+inf; dashed line in A = db only): A = within-bHLH AUROC, CACGTG
      binders vs other labelled bHLH proteins over the scan ranking; B = CACGTG binders in the top
      100; C = bHLH proteins in the top 100. Shaded k = 96–160, dotted k = 128 (Fig 5 B/C CAV).
    - Numbers (db+lit+inf; db in brackets):

      | | k ≤ 16 | k = 32 | k = 64 | k = 96–160 | k = 128 | k = 256 | no PCA |
      |---|---|---|---|---|---|---|---|
      | MYC AUROC | 0.62–0.72 | 0.61 | 0.81 | 0.82–0.86 | 0.85 [0.84] | 0.62 | 0.53 [0.55] |
      | MAD1 AUROC | 0.61–0.68 | 0.56 | 0.71 | 0.80–0.87 | 0.84 [0.83] | 0.62 | 0.63 [0.65] |
      | MYC binders / bHLH, top 100 | 1–6 / 1–9 | 17 / 50 | 28 / 86 | 29–30 / 100 | 30 / 100 | 24 / 100 | 21 / 97 |
      | MAD1 binders / bHLH, top 100 | 0–7 / 0–7 | 16 / 35 | 25 / 82 | 31 / 99–100 | 31 / 100 | 28 / 100 | 22 / 95 |

      Labelled bHLH: 88 (31 binders) in db+lit+inf, 78–79 (24–25 binders) in db.
    - The literature fact-check is not blind to the rankings (the curated proteins came from top-20
      lists); the db-only curve (dashed) gives the same shape.
    - *Superseded, archived:* the C = 1 early-stopped bHLH/bZIP/other count version
      (`pca_sweep_supp/export_pca_sweep_supp_data.py` → `pca_sweep_purity.csv`, now in
      `figure_data/archive_pca_supp_C1_2026-10-05/` with its renders). Step 11 describes it.
11. **Results of the superseded C = 1 early-stopped sweep (bHLH family counts):**

    | | k ≤ 16 | k = 32–64 | k = 128–256 | k ≥ 512 (incl. no PCA) |
    |---|---|---|---|---|
    | MYC bHLH, top 20 | 0–1 | 19–20 | 20 | 17 (3 = JUN/JUND/JUNB) |
    | MYC bHLH / bZIP, top 100 | 4–20 / 0–3 | 77–91 / 0–3 | 92–97 / 3–8 | 83–86 / 7–10 |
    | MAD1 bHLH, top 20 | 0–2 | 16–20 | 20 | 20 |
    | MAD1 bHLH / bZIP, top 100 | 2–15 / 0–3 | 77–83 / 0 | 97 / 3 | 97 / 3 |

    - At k ≤ 16 the ranking fails for both concepts. Top hits are homeodomain, histone H1,
      keratin and VN1R proteins.
    - At k = 32–64 the family is found, but the tail is noisy (HDACs and miscellaneous proteins).
    - From k = 128 up, both concepts are bHLH-dominated, and their known partners stay near the
      top. For the MYC CAV, MAX, MNT, MAD1 and MXI1 rank about 5–25 at every k ≥ 32.
    - Only MYC admits bZIPs at k ≥ 512. With no PCA, JUN / JUND / JUNB rank 9 / 12 / 19;
      at k = 128 they rank 89–96.
    - Heads of adjacent k overlap by 88–96 of the top 100 for k ≥ 64. The concept vectors
      themselves are nearly unrelated (cosine to the no-PCA vector 0.05–0.33 for 32 ≤ k ≤ 512), so
      direction cosine is a poor guide to what a scan returns.
12. **How to frame it (C = 10⁻⁴ sweep, current):** PCA dimension is a resolution setting, shown
    as a curve. Below ~64 components a single-positive CAV loses the concept. Without PCA (or at
    k ≥ 256) the ranking stays bHLH-pure (95–100 of the top 100) but no longer separates CACGTG
    binders from other bHLH proteins (AUROC 0.53–0.63); k = 96–160 does (0.80–0.87). n = 2 related
    concepts, and k was not chosen on held-out data, so call it a plateau, not an optimum.
    *Earlier framing, for the superseded C = 1 sweep (step 11):* present PCA as an option with a lower failure threshold, not as a
    tuned setting.
    - Below roughly 32–64 components a single-positive CAV loses the concept.
    - Above that, rankings are robust. Retaining the low-variance tail (large k, or no PCA) can
      admit neighbouring families, depending on the concept.
    - That happened for MYC but not for MAD1. Why JUN enters for MYC is unexplained; it is not
      a shared zipper, because the MYC span contains none.
    - n = 2 related concepts, so do not claim a general specificity knob.
    - The default (no PCA) is what the main figures use, except Fig 5 B/C. Do not describe
      k = 128 as "optimal": it was not selected on held-out data.
    - **Fig 5 B/C — PROMOTED 2026-10-05: C = 10⁻⁴ converged, PCA-128.**
      - CAV: `converged_C1e-4/build_sweep_converged.py` (same recipe as step 7 but balanced LR
        C = 10⁻⁴, lbfgs `tol = 1e-10`, `max_iter = 1e5`; converged in 23 iterations). Vector =
        name `C1e-4_k0128` in `converged_C1e-4/myc_conv_cavs.npz`, also saved as
        `MYC_HUMAN.fasta_concept_C1e-4_pcaK128/concept_v1.npy` (shared scaler frame).
      - Scores: `converged_C1e-4/myc_conv_scores.tsv`, column `C1e-4_k0128` (step 8 scorer).
      - Export: `pca_sweep_supp/export_fig5BC_myc_data.py` → `figures/figure_data/fig5BC_myc_top20.csv`
        and `MYC_HUMAN.fasta_top20_C1e-4_pcaK128.tsv`; rendered by `figures.R` section
        "Figure 5 B/C" → `cav_fig5BC_myc_top20.{pdf,png}` (identical to the approved preview
        `previews/fig5B_myc_top20_C1e-4.png`).
      - The query MYC (P01106) is excluded from the ranking; included, it ranks 10th.
      - Colour = CACGTG call (`converged_C1e-4/ebox_label/cacgtg_labels_augmented.tsv`, label set
        db+lit+inf: HOCOMOCO v12 / JASPAR 2024, plus the literature fact-check and inference for
        proteins with no database motif). The fact-check is not blind to the ranking: the curated
        proteins were the ones appearing in top-20 lists.
      - Top 20: MXD4, MXD1, MXD3, MXI1, MYCL, MYCN, USF2, USF1, MYCLP1, TFAP4, MAX, HES6, HEYL,
        HEY1, HELT, HEY2, SIM2, HES4, SIM1, HES2. All bHLH; **16/20 CACGTG binders**; non-binders
        TFAP4 (CAGCTG), SIM1, SIM2; MYCLP1 unlabelled (no motif data).
      - Why k = 128 at C = 10⁻⁴: within-bHLH CACGTG AUROC is ≥ 0.80 for k = 96–160 for both MYC
        and MAD1 (0.85 / 0.84 at k = 128), vs 0.53 (MYC) / 0.63 (MAD1) without PCA (db-only labels:
        0.84 / 0.83 vs 0.55 / 0.65); the top 20 shares
        16–20 proteins with k = 128 across that window (`converged_C1e-4/COMPARISON.md`,
        fine-grid section). State it as a plateau, not a tuned optimum.
      - The supplemental sweep (step 10) uses the same C = 10⁻⁴ estimator (promoted 2026-10-05).
    - *Superseded (C = 1, early-stopped):* the MYC top 20 from the k = 128 sweep CAV (column `pca_k128` of
      `pca_sweep/myc_pca_sweep_scores.tsv`; the vector is column 5 of `myc_pca_sweep_cavs.npz`).
      It is the exact-PCA rebuild, not the May randomized-SVD PCA-128 CAV. The two share 97 of
      their top 100 (Spearman 0.935), but they are not identical, so cite this one.
      - All 20 hits are bHLH: MYCN, MYC, MYCL, MYCP1, USF2, USF1, TFAP4, MAX, MAD4, MNT, SRBP1,
        SRBP2, MAD1, MXI1, MLXPL, MAD3, MLXIP, MLX, TFE2, ITF2.
      - MYC itself (rank 2) is the training positive and should be marked as such in the panel.
      - The methods must say this CAV used PCA, cite this supplemental, and give the no-PCA
        comparison: 17/20 bHLH, with JUN/JUND/JUNB at ranks 9/12/19.
13. **Superseded:** the preview `pca_sweep/myc_pca_sweep_preview.png` used the PF00010-only
    bHLH list and undercounts purity. Notes are in `pca_sweep/PCA_KNOB_NOTES.md`.

### 7.6 Hierarchical circle pack of the Pfam motif CAV library
**DECIDED 2026-10-05.** Working directory:
`motif_clustering/results/converged/C0.0001_sharedscaler/threshold_cuts/diffusion_match4_nopromote/`.
- Figure: `pfam_circlepack_clanid_diff_match4_nopromote.{png,pdf}`, a PNG at 300 dpi plus a
  vector `cairo_pdf`, 12.8 × 13.2 in. Placement in the numbered set is `[TK]`.
- The superseded version (old early-stopped library, cut counts 314 – 18,984, top 3 promoted)
  and the cosine/diffusion variants explored on 2026-10-05 are archived in
  `methods_archive/section7.6_superseded_2026-10-05.md`.

**Methods prose:**
> The 27,472 unit-normalised Pfam CAVs (C = 10⁻⁴ library, §1.6) were connected in a mutual
> 15-nearest-neighbour graph by cosine similarity (edges below 0.10 removed), and clustered
> hierarchically by Ward linkage on truncated diffusion-map coordinates approximating the
> random-walk (walktrap) distance at walk length t = 16. To set the cut levels, an
> average-linkage tree on cosine distance (1 − cosine similarity) was cut at cosine distances of
> 0.85, 0.80, 0.70 and 0.50, giving 466, 825, 1,955 and 6,032 clusters. The diffusion tree was
> then cut at these cluster counts, and the resulting four-level hierarchy was drawn as a circle
> pack with ggraph, colouring each family by its Pfam 38.2 clan.

14. **Library:** the converged C = 10⁻⁴ shared-scaler Pfam library,
    **`/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_converged_C0.0001_sharedscaler/cavs/*/L25_concept_v1.npy`**
    (§1.6).
    - All 27,472 motifs are present. Positives whose window starts past the 2,000-token
      embedding cap are dropped, so the 13 concepts that failed in the original library train.
    - Each weight vector is L2-normalised (`cav_cluster.py::load_directions`), giving 27,472
      directions in 1,152 dimensions.
15. **Mutual-kNN cosine graph** (`cav_cluster.py::knn_graph`).
    - Each concept keeps its k = 15 nearest neighbours by cosine similarity, and neighbours
      below 0.10 are dropped.
    - An edge is kept only if it is mutual; its weight is the cosine.
    - Concepts left with no edge are reattached to their single nearest neighbour (1,555).
    - Result: 98,190 undirected edges, mean degree 7.1, 243 connected components. The largest
      component holds 26,486 concepts.
16. **Truncated diffusion-map coordinates, t = 16** (`cav_cluster.py::diffusion_coords`).
    - The target is walktrap's random-walk distance
      r_ij = √(Σ_k (P^t_ik − P^t_jk)² / d_k), with P = D⁻¹W.
    - It is computed from the top eigenpairs of S = D^−½ W D^−½: drop the trivial eigenvector,
      scale each remaining one by λ^t, then by D^−½.
    - Only **343 of 27,472** eigenpairs are used (`--n-eigs 100` widened by the 243 components),
      so Euclidean distance approximates r_ij; it equals it only with the full spectrum. (The
      `cav_cluster.py` docstring says "equal"; do not copy that.)
    - 110 eigenvalues exceed 0.999 (component indicators, kept). The smallest kept eigenvalue
      is 0.979.
    - Output: 342-dimensional coordinates.
17. **Ward linkage** on those coordinates (`scipy.cluster.hierarchy.linkage`, `method="ward"`).
    Merge heights span 1.8e-17 to 2.16.
    - Command, logged in `C0.0001_sharedscaler/cluster.log`, run by
      `results/converged/cluster_Csweep.sh`:
      `cav_cluster.py --distance diffusion --linkage ward --knn 15 --knn-floor 0.10 --n-eigs 100 --walk-steps 16 --cut-mode maxclust --no-shuffle-ids`
    - Tree: `C0.0001_sharedscaler/clusters_diffusion.linkage.npy` (+ `.names.txt`).
18. **Cut levels come from a cosine reference tree of the same library**
    (`threshold_cuts/run_threshold_cuts.sh`).
    - Average linkage on cosine distance (`--distance cosine --linkage average --cut-mode distance`),
      cut at the original fixed cosine-distance thresholds. Output: `threshold_cuts/cosine_average.csv`.
    - The four thresholds used, and the cluster counts they give:

      | Cosine distance | 0.85 | 0.80 | 0.70 | 0.50 |
      |---|---|---|---|---|
      | Clusters | 466 | 825 | 1,955 | 6,032 |

    - The diffusion tree is cut with `fcluster(criterion="maxclust")` at those counts
      (`--linkage-in`). Ward heights: 0.270, 0.0617, 0.0171, 0.00403.
      Output: `threshold_cuts/diffusion_threshcounts_all.csv`, which also has the 211-, 1,327-,
      3,671-, 9,217- and 13,484-cluster cuts; scrambled copy `diffusion_threshcounts_all.scrambled.csv`.
    - **Why borrow counts from a cosine tree:** Ward heights on diffusion coordinates have no
      interpretable scale, while cosine distance does. Earlier versions inherited counts from the
      OLD library's cosine tree (314 – 18,984); these are re-derived on this library.
    - **Why not 0.95:** on this library 0.95 gives 60 clusters, and on the diffusion tree one of
      those holds 27,308 of 27,472 families, so that level is uninformative and is omitted.
19. **No promotion.** Unlike earlier versions, no cluster is promoted (`--top-n 0`), so the
    drawn hierarchy is exactly the four cuts.
    - The largest top-level circle holds 4,264 families; the next largest hold 2,471 and 2,068.
    - At the 825 level the largest cluster is 547.
    - Input to the render: `diffusion_match4_nopromote/clusters_promoted.csv`, identical to the
      four cut columns.
20. **Clan labels:** Pfam 38.2 (`results/figures/Pfam-A.clans.tsv`, 30,134 mappings, 892 clans;
    890 appear in the figure).
21. **Circle-pack layout and rendering** (`diffusion_match4_nopromote/pfam_circlepack_promoted.R`,
    R `ggraph`; generated by `threshold_cuts/diffusion_match4.sh` from the master template).
    - Tree: origin → four cut levels → 27,472 families (36,751 nodes, depth 5), packed with
      `create_layout(layout = "circlepack")`.
    - Leaves have equal weight, so a circle's area scales roughly with its family count.
    - Cluster circles are outlines only. Each family is a dot coloured by one fixed random colour
      per clan (`clan_colours_R.csv`); families without a clan are grey (`#c9c8c0`). No legend.
22. **Clan agreement at the drawn levels** (`threshold_cuts/diffusion_threshold_cuts.scores.csv`,
    Pfam 38.2):

    | Clusters | 466 | 825 | 1,955 | 6,032 |
    |---|---|---|---|---|
    | ARI | 0.132 | 0.212 | 0.157 | 0.090 |
    | AMI | 0.596 | 0.657 | 0.639 | 0.542 |

    - Leaf-order clan contiguity is 64.4, with median clan span 39 leaves. For comparison the old
      library gave 10.7 and 258.
    - **A plain cosine/average tree scores higher on flat-cut ARI** on this library (0.36 at
      distance 0.75), but chains a large core (8,627 families at 0.85, still 4,571 at 0.70). The
      diffusion tree was chosen for its even hierarchy, its higher contiguity, and the case
      studies below. State this honestly if the choice is questioned.
23. **Case studies at the drawn levels** (`threshold_cuts/case_membership.tsv`; CAV-cosine ranks
    in `results/converged/case_studies.tsv`):
    - **DUF2797 (PF10977) → GIY-YIG clan (11 families):** all 11 GIY-YIG families share
      DUF2797's 13-family cluster at the 466, 825 and 1,955 levels, and 9 of 11 at 6,032.
      Median CAV-cosine rank of GIY-YIG families from DUF2797: 6.
    - **CILP_5th (PF23708) → Cystine-knot clan (14 families):** 11 of 14 in a 66-family cluster
      at 825; 8 in a 21-family cluster at 1,955; 8 of 12 at 6,032. Median rank 7.
    - **CLIP (PF12032, PF18398, PF18322, PF18399) → Defensin clan (16 families):** 12 in a
      17-family cluster at 1,955. Median rank 26.

---

## 8. Rendering and software

1. Every panel is rendered by a single R script, `figures/figures.R`, run via
   `figures/00_make_figures.sh`. It reads only flat files from `figures/figure_data/`;
   no model, GPU or large dataset is needed to redraw the figures.
2. R 4.5.3, with `tidyverse` (ggplot2 4.0.3, dplyr 1.2.1), cowplot 1.2.0, ggridges 0.5.7,
   ggrepel 0.9.8, ggrastr 1.0.2, RColorBrewer, `png`, `grid`. The Pfam circle pack (§7.6)
   additionally uses ggraph 2.2.2, tidygraph 1.3.1 and igraph 2.3.3 (from `~/R/4.5`).
   Python (conda environment `core_pkgs4`): Python 3, torch 2.7.1+cu118,
   transformers 4.57.6, scikit-learn 1.8.0, scipy 1.17.1, numpy 2.4.3, pandas 2.3.3,
   statsmodels 0.14.6, scanpy 1.11.5, anndata 0.12.10, umap-learn 0.5.11.
3. Typeface: the house font is set to Arial. Arial is not installed on the analysis
   system, so **Liberation Sans**, a metrically compatible substitute with identical
   glyph widths, is what actually renders. State the substitute.
4. Colour: the Okabe–Ito colourblind-safe palette throughout, with a single-hue
   sequential green ramp for the ordered rank tiers in Fig 2B (the two "barely found"
   tiers shade into purple so they read as a distinct family), and viridis for the
   Fig 5D bitscore heatmap.
5. Dense scatters are rasterised in place with `ggrastr` at 150 dpi while axes and text
   stay vector; PDFs are written with `cairo_pdf`, PNGs at 300 dpi (600 dpi for Fig 4).
6. Figure widths follow journal columns: 3.6 in single column (Fig 2), 5.5 in ≈ 140 mm
   1.5-column (Fig 5), 7.2 in ≈ 183 mm and 7.5 in ≈ 190 mm double column (Figs 3, 6).
7. Panels whose source CSV is absent become void placeholders rather than being dropped,
   so a partial rebuild keeps the intended layout.
7b. **Path portability.** Driver shell scripts resolve every path through
   `config/paths.sh`. The Python scripts that precompute Figure 4's inputs
   (`figures/figure_data/vav_motif_repro/score_*.py`,
   `figures/figure_data/pfam_clan_tsne/render_pfam_clan_umap_3d.py`) resolve
   repository-internal paths relative to their own file, and take the three external
   resources from the environment with the original cluster paths as defaults:
   `ESM_MODEL_DIR` (the ESM-C weights), `PFAM_CAV_DIR` (the 27,472-motif CAV library)
   and `ESM_C_REPO` (the windowing/pooling helpers imported from
   `ahmad_workspace/esm_c`). Behaviour is unchanged when the variables are unset.
   Two exceptions are deliberate: `go_eval/COMMANDS_files/*.sh` and
   `ec_eval/COMMANDS_files/*.sh` are generated artifacts kept as a record of what ran
   (regenerate with `shared_protein/02_make_cav_commands.sh`), and the
   `interface_search/fastas/pairs/bzip/make_bzip_*.py` exploratory scripts feed no
   current figure.
8. Compute: CAV training and all embedding ran on a SLURM cluster, GPU partition, one
   GPU per task, as job arrays. `[TK]` — GPU model and total compute, if the journal
   wants it.
8b. **Pfam CAV training cost (C = 10⁻⁴ library, §1.6; job 24047926, 2026-09-29).**
   - **Scale:** all 27,472 Pfam CAVs, each on up to 1,000 positives plus the 10,000 shared
     negatives, in 1,152 dimensions.
   - **Hardware:** CPU only, on 2 × 48-core AMD EPYC 7642 nodes (UArizona Puma `standard`
     partition). The run was a 40-task array, 4 cores per task, spread over 12 nodes.
   - **Wall-clock:** 21 min for the whole library (12:15–12:36). Tasks took 15–21 min each.
   - **Total compute:** 11.3 task-hours = 45 core-hours, including loading the precomputed
     embeddings from shared storage.
   - **Per CAV:**
     - Median logistic-regression fit 1.3 s on 4 cores (95th percentile 1.8 s, max 2.4 s),
       from a 400-motif sample of `L25_converged_report.json`.
     - Median 23 L-BFGS iterations (max 29) to tol 10⁻¹⁰.
     - Averaged over each task, including I/O, that is about 1.5 s per CAV.
   - **What this excludes:** embedding the positives (GPU, done once, reused across C values).
   - **Other C values are slower:** C = 1 (median 77–82 iterations, about 5 s per motif) took
     1 h 18 min wall-clock. Smaller C converges in fewer iterations.
   - Suggested sentence: "Training all 27,472 Pfam CAVs took 21 min wall-clock (45 CPU
     core-hours; AMD EPYC 7642) from precomputed embeddings, a median of 1.3 s per concept."
9. **Upstream data freshness.** The figure PDFs were rendered 2026-09-15 03:16. The
   `figure_data/` CSVs behind them date from 2026-07-13 (GO/EC evaluation, combined UMAP,
   EC hierarchy) through 2026-09-06 (SERK/CIK) to 2026-09-13 (single cell). Nothing was
   regenerated after the render, so the figures are internally consistent with their
   inputs — the only staleness is the GO CAV retrain in §10.1, which postdates all of it.

---

## 9. Data and code availability (in progress)

Four Zenodo records are planned — background embeddings, and the GO, EC and 27,472-motif
Pfam CAV libraries — each self-sufficient for scoring. Plan and open items in
`distribution/DISTRIBUTION_PLAN.md`. Per-motif score thresholds for the Pfam library are
built by `distribution/pfam_thresholds/build_motif_thresholds.py`.

---

## 10. Open items before submission

1. **RESOLVED 2026-10-01 for Figures 2 and 3.** They are now rendered on the C = 10⁻⁴
   library, which is trained and scored entirely in the shared frame (§1.6), so the
   frame mismatch below no longer affects any Fig 2/3 panel. The record is kept for
   provenance. Remaining follow-ups from the swap:
   - ~~re-pick the GO callout examples~~ done 2026-10-01 (§2.4 C, `callouts_C1e-4/`);
   - regenerate the Pfam figures (Fig 4, §7.6) if the C = 10⁻⁴ estimator is to be
     described paper-wide (§1.6).

   **The GO scaler-frame fix postdates the rendered figures.** 3,930 of 6,677 GO CAVs
   were trained in their own per-concept scaler frame (because
   `shared_protein/02_make_cav_commands.sh` never passed `--scaler-pkl`) while being
   *scored* in the shared frame. Magnitude of the mismatch: `scale_` differs from shared
   by a median of 2.1e-2, max 1.3e-1, relative to the median scale. The retrain
   (`distribution/retrain_go/01_retrain_go_cavs.sbatch`, 40-task array) ran 2026-09-16 to
   2026-09-17; the Figure 2/3 PDFs were rendered 2026-09-15 03:16 and
   `go_scaler_frames.tsv` still reports 3,930 own / 2,747 shared, i.e.
   `02_verify_retrain.sh` has not been re-run to confirm own = 0. **Either re-render
   Figures 2 and 3 on the retrained library, or the methods cannot claim a single shared
   coordinate system for the GO results.** Classify frames via `n_samples_seen_` (an
   own-frame scaler saw 10,000 + n_positive); md5 of the joblib file is not a valid test.
   EC and Pfam are each internally consistent and unaffected.
   **The defect bites Fig 3C/D/E differently** (letters as of 2026-10-06: C UMAP, D cosine, E rank). Figures 2 and 3A have a genuine
   train/score frame mismatch. Figures 3D and 3E compute cosine similarity directly
   between GO and EC *concept vectors*, and 3C UMAPs 6,555 GO vectors alongside 952 EC
   vectors — there is no scoring step at all. There, the problem is that ~59% of the GO
   weight vectors are expressed in a different basis from the shared-frame GO CAVs and
   from every EC CAV, so the geometry being compared is inhomogeneous. The discriminating
   check is whether the retrained library moves 3E's rank-1 percentage or 3D's median
   separation; if it does, those panels need re-rendering regardless of what is decided
   for Figure 2.
2. **Figure 4A's library and n.** Reconcile the `figures.R` header (7,688 of 7,692) with
   the 27,472-motif library used by every other Fig-4 component (§4.1 item 1).
3. **Fig 6 DE follow-ups (after the 2026-09-29 switch to pseudobulk DESeq2, §6.3).**
   (a) The two single-cell supplements (§7.1, §7.2) still read the mixed-model DE. Switching
   either is a one-line `DE_DIR` change in its exporter, but the x axis and gene universe
   change with it. (b) `FIG6_DE_TEXT.md`'s Results/Methods prose describes the mixed
   model (r = 0.60/0.73/0.27, "pseudobulk failed"). Rewrite it from §6.3 before drawing
   paper text from it. (c) `de_methods.py::run_deseq2_pair` still has the per-donor
   aggregation bug (§6.3 item 18). (d) The earlier instruction here to "correct" panel D
   to 514/486 was wrong. 1,000/1,000 is the analysed count and 514/486 the plotting
   subsample (§6.2 item 9).
4. **Unbuilt panels:** Fig 2A is a blank reserved slot in the assembled grid and needs
   content or removal. Fig 6A is a blank bar awaiting an externally supplied schematic.
   Fig 5 A/B/C are placeholder boxes (A and B are Illustrator; C is TBD and was intended
   to carry the MYC case). The MYC top-20 panel is now rendered standalone as
   `cav_fig5BC_myc_top20.{pdf,png}` (2.6 × 3.6 in) for assembly; it is not yet in the `cav_fig5` grid. Fig 4A's clan labels and Fig 5D's interface-class axis label are
   added by hand after rendering — worth one sentence so a reader knows they are not
   script output.
5. **Remaining `[TK]` values:** QuickGO download date, GO OBO release, DeepGo-SE
   version/checkpoint, ec2go snapshot date, EC training-overlap count, count of EC concepts carrying a
   proFAB `target_` split (§1.8), compute totals.
   (Software versions and the GO annotation cutoff are now filled in.)
6. **A motif-clustering methods section may be needed.** `motif_clustering/` contains a
   full four-way dendrogram comparison on the 27k library (cosine/average,
   cosine/ward, diffusion/ward, snn/average, all cut with `fcluster` at matched cluster
   counts and scored by ARI / AMI / per-clan F1 against Pfam clans), plus the circle-pack
   visualisations and the context-leakage analysis. The circle pack itself (C = 10⁻⁴ library,
   diffusion t = 16, cut at cosine-derived counts 466 – 6,032) is documented in §7.6. Its figure number is still
   `[TK]`, and the dendrogram comparison and context-leakage analysis remain undocumented. One finding
   worth carrying into the text regardless: the preATP/ATP-grasp co-clustering is **not**
   context leakage but a four-hop seam of clan-misannotated and merged-domain Pfam
   families, with eukaryotic glutathione synthetase as the key case.
7. **Which MYC CAV Figure 5 B/C uses — DECIDED 2026-10-05: C = 10⁻⁴ converged PCA-128, top 20**
   (§7.5 step 12; replaces the 2026-09-29 decision for the C = 1 early-stopped k = 128 sweep CAV).
   The PCA supplemental was swapped to the matching C = 10⁻⁴ CACGTG sweep the same day (§7.5 step 10).
   Remaining: state in the methods that this CAV uses PCA, which the rest of the paper does not.
   Do not use the Jul-8 `MAD1_HUMAN.fasta_query_results.tsv` in a panel: it matches the k = 128
   sweep only at Spearman 0.87, so its provenance is unconfirmed.
