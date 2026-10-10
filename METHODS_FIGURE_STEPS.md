# Extended methods — how each figure in the paper was produced

Companion to the Online Methods in `text.txt`. For each panel it gives the scripts, inputs,
parameters and the numbers the text cites, so a reader can be pointed at the code.
Panel letters are the **final letters in `text.txt`**. Several file and object names in the
code still carry older letters or numbers (table in §0).

The full pre-cleanup notes (dated decision logs, superseded libraries and analyses, demoted
supplements) are archived at `methods_archive/METHODS_FIGURE_STEPS_full_2026-10-09.md`.
Open discrepancies between the code and `text.txt` are listed in §9.

Conventions:
- `$BIOCAV_REPO` is the method repository (github.com/mcwhitelab/bioCAV).
- Paths are relative to this repository unless absolute.
- Every driver sources `config/paths.sh`.
- Large inputs live on `/xdisk` and are not in git.

---

## 0. Figure inventory

| Paper | Produced by | Rendered file(s) | Names still in the code |
|---|---|---|---|
| Fig 1 | Illustrator schematic | — | — |
| Fig 2 | `figures/figures.R` | `figures/figures/cav_fig2.{pdf,png}` | R grid has A blank, B ranks, C/D violins. Paper: A, B schematics; C ranks; D AUROC; E AUPR |
| Fig 3 | `figures.R` | `cav_fig3.{pdf,png}` | panel objects `p_2d`, `p_2e_solo`, … |
| Fig 4 A | `motif_clustering/…/diffusion_linkages/` (§4.1) | `pfam_circlepack_clanid_average_mainfrac_045-005_promote2.{pdf,png}` | not drawn by `figures.R`; `cav_fig4` panel A is an obsolete 3-D UMAP |
| Fig 4 B, C | `figures.R` | `cav_fig4.{pdf,png}` (panels B, C) | — |
| Fig 4 D–J | structural analyses outside this repository (§4.4) | Illustrator / PyMOL | — |
| Fig 5 A, C | Illustrator schematics | — | — |
| Fig 5 B | `figures.R`, "Figure 5 B/C" block | `cav_fig5BC_myc_top20.{pdf,png}` | data `fig5BC_myc_top20.csv` |
| Fig 5 D, E | `figures.R`, Figure 5 block | `cav_fig5.{pdf,png}` (panels D, E; A–C are placeholders) | data `fig6_heatmap.csv`, `fig6_densities.csv`, `fig6_auroc.csv` |
| Fig 6 A, B | schematics | — | — |
| Fig 6 C, D | `figures.R`, Figure 6 block | `cav_fig6.{pdf,png}` | exporter `export_fig5_case_studies.py`; renders three cases (§9) |
| Supp File 1 | §4.1 | `clusters_average_mainfrac_steps005{,.scrambled}.csv` | — |
| Supp File 2 | §4.4 | not in this repository (§9) | — |
| Supp (not cited in text) | `figures.R` | `fig_supp_pca_sweep`, `fig_supp_projection_comparison`, `fig_ec_eval` | — |

**Redrawing.** `figures/00_make_figures.sh` runs `figures/figures.R` from `figures/`.
- The script reads only flat files in `figures/figure_data/`, plus `figures/figures/pfam_clan_umap_3d_mosaic.png`. Without that PNG the whole Fig 4 block is skipped.
- Each block is wrapped in `file.exists()`, so a missing input skips that figure silently. Check the log for "Skipping".
- All computation happens upstream. The R script itself computes only the per-panel Spearman correlations.

---

## 1. Training concept activation vectors (shared protein preamble)

### 1.1 Model and embedding
- **Model.** ESMplusplus_large, i.e. ESM-C 600M: 36 layers, hidden size 1152, loaded with `trust_remote_code=True`. ESM-C 600M is MIT-licensed. The `NOTICE.txt` beside the local weights is stale.
- **Embedding script.** `mcwlab_utils/hf_embed_new.py`, the copy at `/groups/clairemcwhite/claire_workspace/github/mcwlab_utils/` (commit `f733f2cc3cb1`).
  - Layer `-l -11` = `hidden_states[25]`, which is **"layer 26"** in the paper. The paper counts the embedding layer as layer 1, so a file named `L{k}` is displayed as layer k + 1.
  - Mean pooling, batch size 1, `--max_length 2048`.
  - `--get_aa_embeddings` for anything that is span-pooled.

### 1.2 Background and shared frame
- **Background.** 10,000 random UniRef50 proteins: `shared_protein/00_get_negatives.sh` → `shared_protein/reference_population/neg_10000.pkl`.
- **Shared scaler.** `$BIOCAV_REPO/scripts/fit_reference_scaler.py` fits one `StandardScaler` on the background → `scaler_v1.pkl`. Every library CAV is trained and scored in this frame.

### 1.3 Estimator
- **Script.** `$BIOCAV_REPO/scripts/train_cav_from_embeddings.py` → `core/src/train_cav.py`. Spans are pooled by `prepare_embeddings.py`.
- **Model.** scikit-learn `LogisticRegression`: L2, `solver='lbfgs'`, `class_weight='balanced'`, `random_state=42`, fit to convergence (`tol=1e-10`, `max_iter=100000`).
- **Libraries** (GO, EC, Pfam): C = 10⁻⁴, which is ≈ 1/n for n ≈ 10⁴ training embeddings. sklearn sums the loss over samples, so the effective regularization scales with 1/(C·n).
- **Concept vector.** `coef_[0] / ‖coef_[0]‖`. The intercept and the classifier are discarded.
- **Scoring.** Standardize with the training scaler, then take the dot product with the unit CAV. The result is a signed projection, not a probability.
- **Why convergence matters.** With sklearn's default `tol=1e-4`, lbfgs stops after 12–50 iterations, and the direction then depends on the stopping rule. All libraries in the paper are converged fits.

### 1.4 Libraries used

| Library | n | Positives | Location |
|---|---|---|---|
| GO | 6,678 | proFAB `random_positive_train_max1000.span` | `/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_sharedscaler` |
| EC | 953 | same | same |
| Pfam | 27,472 | ≤1,000 InterPro-annotated domain spans | `/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_converged_C0.0001_sharedscaler` |

- **GO/EC build.** `distribution/retrain_go_ec_C1e-4/`. The pipelines read the GO/EC CAVs through a symlink overlay of proFAB, `/xdisk/clairemcwhite/clairemcwhite/go_ec_cavs_C1e-4_overlay`, built by `04_build_overlay.py`.
- **EC 3.4.22.28.** All 76 positives are 2,052–2,432 aa polyproteins. Each was truncated to its first 2,046 aa, which keeps the 3C protease region (`03_ec_3-4-22-28_truncated.sh`).
- **Pfam build.** `distribution/retrain_pfam_converged/`, which imports constants from `distribution/retrain_pfam_pca128/train_pfam_pca128.py`.
- **Training cost** (job 24047926), as quoted in the text:
  - all 27,472 Pfam CAVs in 21 min wall-clock;
  - a 40-task × 4-core array on AMD EPYC 7642, 45 core-hours in total;
  - median fit 1.3 s per CAV (23 lbfgs iterations);
  - embedding time excluded.

### 1.5 proFAB data and splits
1. `shared_protein/01_download_profab_datasets.sh` downloads the GO (10 parts) and EC (4 parts) archives.
2. `shared_protein/format_proFAB.sh` → `format_proFAB.py --max-n 1000 --seed 42` writes `{group}_{polarity}_{split}_max1000.span` files.
3. `shared_protein/02_make_cav_commands.sh` emits the per-concept training commands. The original library was run as SLURM arrays (`go_eval/SBATCH_files/array_GO.sbatch`, `ec_eval/SBATCH_files/array_EC.sbatch`). That run also cached the proFAB test-set embeddings that the LLR uses.

Which proFAB files are used:
- **Training positives:** `random_positive_train_max1000.span` only.
- **Training negatives:** the shared 10k background. proFAB's own negatives are not used for training.
- **LLR Gaussians and AUROC/AUPR negatives:** `random_{positive,negative}_test_max1000.span`.
- **`temporal_` and `target_` splits:** generated by `format_proFAB.py` but **never read**.

The paper's GO "temporal" hold-out is a separate QuickGO download (§2.1) and has nothing to do with proFAB's `temporal_` split.

### 1.6 LLR
- **Gaussians.** For each concept, `eval_go_temporal.py` / `eval_ec.py` fit N(pos) to its CAV scores on `random_positive_test_max1000.span` and N(neg) on `random_negative_test_max1000.span`.
  - These are **test** positives. The docstring of `eval_go_temporal.py` says "training positives"; the code is right and the docstring is wrong.
- **Definition.** LLR(x) = log N(x; pos) − log N(x; neg), via `scipy.stats.norm.logpdf`.
- **Gate.** LLR ≤ 0 counts as "no prediction" in Fig 2C and Fig 3A.
- **Library dependence.** The Gaussians come from whichever CAV is evaluated, so the gate was recomputed for the C = 10⁻⁴ library.

---

## 2. Figure 2 — GO annotation recovery

All panels come from `figures/rerun_fig23_C1e-4.sh`, which runs each panel's pipeline against the C = 10⁻⁴ overlay and writes `figures/figure_data_staging/C1e-4/`. The output was promoted to `figures/figure_data/`.

### 2.1 Temporal validation set
- **Source.** QuickGO GPAD files in `go_eval/static/go_annotations/{mf,bp,cc}.gpad`.
  - Downloaded 2026-05-27, filter `evidenceCode=ECO:0000314&evidenceCodeUsage=descendants`.
  - GO release 2026-04-27; ontology from `go.obo` in the same folder.
- **Gold standard** (`$BIOCAV_REPO/specific_scripts/eval_go_temporal.py`):
  1. `DB == UniProtKB` only.
  2. Annotation date after 2021-12-31 (`--date-cutoff 20211231`).
  3. Unique (protein, term) pairs.
  4. Accessions normalized by stripping `-PRO_…` and isoform suffixes.
  5. Terms need ≥ 5 proteins (`--min-n 5`) and a trained CAV.
- **Known-positive exclusion** (`go_eval/06_eval_exclude_known_positives.sh`, `--exclude-known-positives`). This drops a pair when the protein is already in that term's proFAB random train or test positives, then re-applies `--min-n 5`.
  - The term set is pinned by `--go-terms-file static/fig2_go_terms_{ont}.txt` to 136 MF, 256 BP and 162 CC terms.
  - Outputs go to `go_eval/results/temporal_eval_{ont}_noknown/` and `figures/figure_data_staging/noknown/`. `rerank_go.py` reads the latter.
- **Caveat worth a sentence.** The IDA-descendant filter also pulls in IPI, so "protein binding" (GO:0005515) is 90% of post-2021 MF rows. It drops out only because proFAB has no CAV for it.

### 2.2 Comparators
- **DeepGO-SE.**
  - Run locally from `/xdisk/clairemcwhite/deepgose_latest.sif` with `apptainer exec --nv`, model data in `/xdisk/clairemcwhite/data`.
  - Command: `python predict.py -if <fasta> --device cuda --batch-size 1`. Every other setting is the default, including the **0.1 output threshold**.
  - It was run on both the validation proteins and the proFAB negative test proteins (`go_eval/SBATCH_files/deepgose_{val_proteins,negatives}.sbatch`). `go_eval/02_format_outside_tool_preds.sh` collects the output into `go_eval/outside_tools/deepgose_preds_{bp,cc,mf}.tsv`.
  - Any (protein, term) pair it did not emit is scored 0. Positives and negatives share the threshold, so the zero-fill is symmetric.
- **DIAMOND baseline** (since 2026-10-09). DIAMOND 2.0.9 `blastp --ultra-sensitive --max-target-seqs 0 --evalue 1`, self-hits removed.
  - Search: `go_eval/rank_C1e-4/run_diamond_ultra.sbatch`, all validation proteins against the union of all CAV training positives (`train_positives_all.dmnd`). Hits are in `go_eval/blast_baseline/C1e-4_extra_ultra/`, which `rerank_go.py` reads.
  - Each candidate term is scored by the best bitscore to its CAV training positives. Ties get the average rank (`go_eval/blast_baseline/blast_ranks.py`).
  - "No prediction" means no hit (e ≤ 1) to the true term's training positives.
  - GO:0032436's positives were all already in the database, so its member rows (from `extend_blast.sh`) are reused as they are.
  - **Previous setting** (`--more-sensitive`, e ≤ 1e-3; `C1e-4_extra/`, `extend_blast.sh`) and the comparison with BLASTP (BLAST+ 2.17.0), jackhmmer (3 iterations) and ultra-sensitive at e ≤ 1e-3: `go_eval/blast_baseline/homology_compare/` (README, `summary.md`). No Rank 1 or top-10 cell moves by more than about 3.3 points between any of them. The only exception is EC L4 Rank 1, 100 → 95.2% under BLASTP and jackhmmer (1 of 21 pairs). Ultra-sensitive at e ≤ 1 is the strongest search in most cells.
  - Backup of the pre-switch data, `figures.R` and this file: `figures/figure_data/_backup_pre_fig23_diamond_ultra_2026-10-09/`.

### 2.3 Panel C — rank of the true term, by family count
- **Driver.** `go_eval/rank_C1e-4/rerank_go.py`, which imports `go_eval/rank_all_terms/compute_ranks_all_terms.py` and `blast_ranks.py`. Output: `figure_data/go_rank_families_{mf,bp,cc}.csv`.
- **Candidates** (column `*_rank_ont_filt_pf`):
  - all non-obsolete, same-ontology GO CAVs: MF 1,114, BP 4,298, CC 660;
  - 484 obsolete terms are dropped via `parse_obo(drop_obsolete=True)`;
  - also removed: the true term's `is_a`/`part_of` ancestors, and the protein's other GPAD and proFAB-positive terms, closed over their ancestors.
- **Rank.** 1 + the number of candidates scoring strictly above the true term. CAVs are ranked by raw projection. DeepGO-SE is ranked over the same candidates.
- **No prediction.** CAV LLR ≤ 0; DeepGO-SE did not emit the true term; DIAMOND found no hit.
- **Family bins.** A separate DIAMOND run, not the ultra-sensitive baseline above and unchanged by it: `go_eval/diversity/term_diversity.py` runs DIAMOND all-vs-all on each term's training-positive FASTA. An edge needs ≥ 30% identity **over ≥ 50% of the shorter sequence**; families are the connected components.
  - The bin edges <20 / 20–59 / ≥60 come from `go_eval/blast_baseline/family_cutoff_screen.py` (`preview_family_screen.png`).
  - Pairs per bin: 818 / 4,158 / 3,038. Pairs per ontology: MF 1,223, BP 2,739, CC 4,052.
- **Result** (Rank 1 / top 10 %):

  | | CAV | DeepGO-SE | DIAMOND |
  |---|---|---|---|
  | MF | 17.8 / 41.2 | 10.1 / 22.7 | 12.0 / 36.7 |
  | BP | 5.5 / 20.9 | 1.5 / 6.2 | 3.0 / 13.4 |
  | CC | 8.0 / 26.5 | 3.5 / 15.7 | 8.8 / 30.5 |

  - DIAMOND here is ultra-sensitive at e ≤ 1. With the old setting (`--more-sensitive`, e ≤ 1e-3) it was MF 11.4 / 34.3, BP 2.9 / 12.5 and CC 8.3 / 29.2.
  - CAV leads at < 60 families in every ontology. At ≥ 60 families, and for CC overall, DIAMOND is level or ahead (BP ≥ 60 top 10: DIAMOND 10.6 vs CAV 10.0), so do not write that CAV leads DIAMOND throughout.

### 2.4 Panels D, E — per-term AUROC and AUPR
- **Sources.** CAV values come from `summarize_temporal_eval.py` (columns `auc_/aupr_val_vs_test_neg`). DeepGO-SE values come from `compare_tool_temporal.py` (`tool_auc/tool_aupr`). Both feed `figure_data/temporal_tool_comparison_{mf,bp,cc}.csv` Both scripts are run by `go_eval/06_eval_exclude_known_positives.sh`, and again on the C = 10⁻⁴ overlay by `figures/rerun_fig23_C1e-4.sh`.
- **Positives and negatives.**
  - Positives are the term's known-positive-excluded validation proteins: 3–409 per term.
  - Negatives are its proFAB `random_negative_test` proteins: 9–77 per term.
  - Both methods are scored on identical sets.
- **Metrics.** `roc_auc_score` and `average_precision_score`, one value per term (MF 136, BP 256, CC 162 terms). Chance AUPR is each term's positive fraction (median 0.21–0.28).
- **Medians:**

  | | CAV AUROC / AUPR | DeepGO-SE AUROC / AUPR |
  |---|---|---|
  | MF | 0.938 / 0.882 | 0.571 / 0.374 |
  | BP | 0.838 / 0.733 | 0.500 / 0.318 |
  | CC | 0.897 / 0.857 | 0.556 / 0.463 |

  - CAV AUROC exceeds DeepGO-SE on 87% (MF), 89% (BP) and 94% (CC) of terms.

---

## 3. Figure 3 — EC recovery and GO–EC relationships

### 3.1 Panel A — EC rank (124 pairs)
- **Gold standard.** EC-Bench `test_ec.csv`, copied verbatim into `ec_eval/static/ec_annotations/EC-Bench/`.
  - Source: github.com/dsaeedeh/EC-Bench, commit `5e98a2c`; Swiss-Prot 2023-01 entries absent in 2018-02.
  - Reformatted by `ec_eval/01_reformat_gold_standards.sh`. Embeddings: `test_ec.fasta.pkl`.
- **Exclusions.**
  - The protein is already a proFAB train **or test** positive of the true EC (209 → 138 pairs). Implemented in `ec_eval/blast_baseline/ec_ranks_noknown.py`.
  - UniProt "fragment" entries (→ 124: L3 103, L4 21). List from `EC-Bench/fetch_uniprot_fragments.py`, applied in `ec_eval/rank_C1e-4/rerank_ec.py`.
- **Not a temporal hold-out.** EC-Bench is a 2018 → 2023 split, but relative to EC-Bench's own training set. Our CAVs are trained on proFAB (a 2021 snapshot), so separation comes from the overlap removal above, not from dates.
- **Candidates.** Same-level EC CAVs (153 at level 3, 733 at level 4), minus the protein's other EC-Bench and proFAB-positive ECs at that level.
  - Projection scores are not comparable across EC levels, so pooling levels would let deeper CAVs win on scale alone.
- **DIAMOND.** Same settings as §2.2 (ultra-sensitive, e ≤ 1), against each EC CAV's training positives (`ec_eval/blast_baseline/C1e-4_extra_ultra/{ec_term_members,test_vs_ectrain}.tsv`, read by `rerank_ec.py`).
- **CLEAN-SupCon.** Not run by us. We used the `clean-supcon` column of EC-Bench's `results_test_100_go_3.csv`, scored as exact match.
  - The pairs file `figure_data/ec_specificity_ranks.csv` (with `clean_supcon_correct`) is written by `ec_eval/04_compare_ec_tools.sh` and read by `rerank_ec.py`.
- **Output.** `figure_data/ec_ranks_samelevel_noknown.csv`.
- **Result:**

  | | CAV rank 1 / top 3 | CAV no prediction | DIAMOND rank 1 | CLEAN correct |
  |---|---|---|---|---|
  | L3 (n = 103) | 55.3 / 75.7% | 13.6% | 77.7% (old setting 75.7%) | 71.8% |
  | L4 (n = 21) | 90.5 / 100% | 0% | 100% | 90.5% |

  - The two level-4 CAV misses are fungal PALs (4.3.1.24), each ranked 2nd behind the HAL CAV (4.3.1.3). The PAL CAV's training set is mostly plant enzymes.

### 3.2 Panel B — score along the EC hierarchy
- **Driver.** `hierarchy_decay/00_run.sh` → `ec_hierarchy_decay.py`, written to `hierarchy_decay/results_C1e-4/` and then `figure_data/ec_hierarchy_decay.csv`.
- **Protein selection.** Proteins with a level-4 EC whose whole lineage has CAVs, and that are leak-free at every level (in no train or test span of any lineage CAV).
  - Rule: one protein per level-2 class, the first accession alphabetically.
  - Shown: C0HLQ6 (1.2.1.12), P0DQP1 (2.1.1.72), P0DUB9 (3.1.1.3), P0DUE1 (3.2.2.6), A0A345BJN1 (4.3.1.24).
  - C0HLQ6 is a 16-aa UniProt fragment, but it is the only leak-free class-1 protein.
  - The 3.1 pick does not follow the alphabetical rule (A0A345BJN2 would). See §9.
- **Plotted values.** Raw projection scores at levels 1–4.

### 3.3 Panel C — combined GO/EC CAV UMAP
- **Pipeline.** `go_eval/05_visualize_umap.sh` → `motif_clustering/combined_cav_umap_3d.py --dims 2 --coords-csv-out`, written to `figure_data/go_ec_combined_umap.csv`.
  - `combined_cav_umap_3d.py` imports `motif_clustering/cav_viz_pattern.py` and `ec_eval/scripts/ec_names.py`.
- **Namespace lists.** Built by `go_eval/scripts/build_go_namespace_lists.py`; the C = 10⁻⁴ copy is in `go_eval/results/go_namespace_lists_C1e-4/`.
- **Size.** n = 7,509 vectors: EC 953, MF 1,161, BP 4,688, CC 707.
- **Reproducibility.** UMAP uses `random_state=42` but is not bit-reproducible across runs.

### 3.4 Panels D, E — GO–EC cosine agreement
- **Pipeline.** `go_ec_crossval/00_run_validation.sh` → `go_ec_cosine_validation.py`, written to `go_ec_crossval/results_C1e-4/` and `figure_data/go_ec_cosine_{pairs,background_sample}.csv`.
  - Curated mapping: `go_ec_crossval/static/ec2go.txt`.
  - The background is a 20,000-value sample drawn with `rng(0)`.
- **What is computed.** For each of the 239 GO terms that have a curated EC partner and CAVs on both sides:
  - the cosine to that partner;
  - the partner's rank among all 953 EC CAVs.
- **Results.** Matched median cosine 0.831 vs background 0.021. The partner ranks first for 88.3% of terms and in the top 5 for 98.7%.
- **Control.** GO and EC CAVs share no supervision: different label sources and different positive sets.

---

## 4. Figure 4 — Pfam CAVs

### 4.1 Panel A and Supplemental File 1 — circle pack of 27,472 Pfam CAVs
Working directory: `motif_clustering/results/converged/C0.0001_sharedscaler/diffusion_linkages/`. Results are gitignored.

1. **Load.** Unit-normalize the C = 10⁻⁴ Pfam CAVs (`L25_concept_v1.npy`), with `motif_clustering/cav_cluster.py::load_directions`.
2. **Mutual k-NN cosine graph** (`knn_graph`).
   - k = 15, neighbours below 0.10 dropped, mutual edges only.
   - Isolated CAVs are reattached to their single nearest neighbour.
   - Result: 98,190 edges and 243 components; the largest has 26,486 CAVs.
3. **Diffusion-map coordinates** (`diffusion_coords`).
   - t = 16, from the top eigenpairs of D^−½WD^−½, giving 342 dimensions.
   - These **approximate** the random-walk distance; the `cav_cluster.py` docstring's "equal" is wrong.
   - Saved as `phi.npy`. ARPACK is not bit-reproducible, so keep this file.
4. **Average linkage.** Euclidean distance on `phi`, saved as `linkage_average.npy` (`linkage_compare.py`).
5. **Heights.** Expressed as fractions of the height at which the largest component is complete (0.0935); above that, merges only join disconnected components.
   - Cuts at 0.45 / 0.20 / 0.10 / 0.05 give 455 / 833 / 1,722 / 3,682 clusters (`average_mainfrac_045-005_promote2/make_cuts.py`).
6. **Promotion.** The two nested coarse blobs (20,446 and 12,069 CAVs) are replaced by their children before layout (`motif_clustering/promote_largest_clusters.py --top-n 1`, run twice).
   - As a result, the drawn rings mix levels.
7. **Render.** R ggraph circlepack, coloured by Pfam 38.2 clan (`Pfam-A.clans.tsv`; `motif_clustering/update_pfam_clans_38_2.py`).
   - Template: `diffusion_linkages/inputs/pfam_circlepack_template.R` (vendored copy of the old `retrained_shared_scaler/promoted_top3/old_diffusion` template, with clan colours `inputs/clan_colours_R.csv`). Family metadata comes from `inputs/pfam_family_metadata.csv`.
   - Script: `render_average.sh` (points at this variant; header fixed 2026-10-09).
   - Output: `average_mainfrac_045-005_promote2/pfam_circlepack_clanid_average_mainfrac_045-005_promote2.{pdf,png}`.
   - **Layout (chosen 2026-10-09):** the ggraph default circlepack layout on the input sorted by `dendrogram_order`, with no sibling shuffling.
     - It is the `_default` render in `average_mainfrac_045-005_promote2/jostle/` (`pfam_circlepack_jostle.R` with `pack_seed <- NA`), copied over the figure, and its layout is copied into `circlepack_layout_average_mainfrac_045-005_promote2.rds`.
     - Circle packing depends on input row order. The earlier layout (unsorted, Ward row order) put a large circle in the centre; it is archived in `archive/2026-10-09_circlepack_cleanup/`.
     - Random sibling orders (`pack_seed` 1, 2; `sort.by` = runif) were tried and not used; archived in `archive/2026-10-09_circlepack_cleanup/`.
     - The R script reuses the cached `.rds`; delete it to re-pack.
   - **Interactive web version** (for mcwhitelab.org/pfam-circlepack/): `average_mainfrac_045-005_promote2/web/index.html`, a single self-contained file of about 3 MB.
     - Build: `export_layout.R` writes the cached layout `.rds` to CSV, then `build_page.py` produces the page.
     - Page code is the existing lab page (`template_mcwhitelab_2026-10-09.html`) with new data. Colours come from `clan_colours_R.csv`, as in the PNG.
     - Hover shows the cluster IDs at h0.45 / 0.20 / 0.10 / 0.05, taken from the Supplemental File 1 scrambled CSV, so they match the supplemental table.
8. **Supplemental File 1.** `clusters_average_mainfrac_steps005.csv` (`make_allsteps_csv.py`).
   - 24 cuts, unpromoted, in columns `cut_hXXX_kNNNNN` (XXX = fraction × 100):
     - every 0.05 from 1.00 down to 0.05 (257 – 3,682 clusters);
     - plus 0.04 / 0.03 / 0.02 / 0.01, giving 4,566 / 5,911 / 8,214 / 12,431 clusters (added 2026-10-09).
     - The finest levels are mostly small clusters: singletons number 1,827 / 2,639 / 4,338 / 7,911 at the four added levels.
     - Adding the finer levels left the 20 original columns and their scrambled IDs unchanged. The filename keeps `steps005`.
   - Rows are sorted by `dendrogram_order`, which is the average-tree leaf order (`leaves_list(linkage_average.npy)`).
     - Before 2026-10-09 the rows kept the Ward tree's row order, and on 2026-10-08 the column briefly held Ward values.
     - The figure layout was re-packed from the sorted files on 2026-10-09 (step 7).
   - The distributed copy is `.scrambled.csv`, which randomizes cluster IDs within each column (`motif_clustering/scramble_cluster_ids.py`, seed 42). Membership is identical.
9. **Clan agreement** as drawn: AMI 0.66 / 0.67 / 0.63 / 0.62 across the four rings.

### 4.2 Panel B — VAV_DROME across layers
- **CAVs.** For PF00621, PF00130, PF00017 and PF00018, trained at all 36 layers into `figures/figure_data/vav_motif_repro/tcav_outputs_all_layers/` (log: `train_all_layers.log`).
  - Negatives: the shared 10k background.
  - These are a separate small library, not the converged Pfam library.
- **Scoring.** `figures/figure_data/vav_motif_repro/score_vav_all_layers.py`, applied to Q9NHV9:
  - window length = the motif's median training-span length, stride 1;
  - each window: standardize, then dot product;
  - each residue gets the mean score of the windows covering it.
  - Output: `vav_layerwise_scores.json`.
- **Layer choice.** Signal is strongest at layers 21–31; layer 26 (L25) is used downstream.

### 4.3 Panel C — human Plexin-A1
- **Scoring.** `figures/figure_data/vav_motif_repro/score_bigreceptors.py`, with the same window procedure at layer 26.
  - Output: `extra_proteins/candidates_bigreceptors_scores.json`.
  - The rendering helper `figures/draft_fig_extra_proteins_common.R` is sourced by `figures.R`.
- **Library.** The script's default `PFAM_CAV_DIR` is the original 27,472-motif library at `/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs` (early-stopped, per-motif scaler), not the converged library of §4.1 (§9).
- **Score scale.** Raw scores are shown; the domains stand out above ≈ 5.

### 4.4 Panels D–J and Supplemental File 2 — structural follow-up
- **Candidate selection.** Candidates were picked from Pfam CAV neighbourhoods. `motif_clustering/results/converged/case_studies.py` records, for example:
  - DUF2797 → GIY-YIG clan, median neighbour rank 6, sharing an 11/13 cluster;
  - CILP_5th → Cystine-knot clan, median rank 7.
- **PF18419 (E).** The ATP-grasp_6 reading is documented in `motif_clustering/context_leakage/FINDINGS.md`: a clan-misannotated, merged-domain seam involving eukaryotic glutathione synthetase.
- **Fold comparisons (D–F, J)** were run with the `pfam-fold-overlay` procedure:
  - one representative domain per family (best experimental structure, otherwise AlphaFold);
  - all-vs-all TM-align, then superposition.
- **Zn-site screen (G–I, Supplemental File 2)** used the `pfam-zinc-site-screen` procedure:
  - AlphaFold DB models restricted to InterPro Pfam bounds, ≤ 3 members per family;
  - sites of 4 Cys/His ligand atoms, centroid distance ≤ 3.0 Å, pairwise 2.6–5.0 Å;
  - a site must recur in ≥ 2 members;
  - controls: 58/60 known zinc-module families detected, 0/148 random families.
- **Location.** The run outputs for D–J and Supplemental File 2 are **not in this repository** (§9).

---

## 5. Figure 5 — interface specificity

Working directories are under `interface_search/fastas/pairs/`, which is gitignored (§9).

### 5.1 Panels A, B — single-example MYC CAV
- **Concept.** MYC_HUMAN (P01106) residues 366–420, span-pooled from the full-length embedding. The span is the bHLH, excluding the zipper.
  - Embedding, span and proteome search: `interface_search/run_spanfile_myc_contacts.sh`.
- **CAV** (`myc_mad/converged_C1e-4/build_sweep_converged.py`):
  1. standardize with the shared scaler;
  2. exact PCA (`svd_solver='full'`) fit on the 10k standardized background, keeping k = 128;
  3. balanced L2 logistic regression with C = 10⁻⁴, converged;
  4. fold back to 1152 dimensions (`components_.T @ coef`) and unit-normalize.
  - Saved as name `C1e-4_k0128` in `myc_conv_cavs.npz`.
- **Why PCA.** Without it, a C = 10⁻⁴ single-positive CAV is the standardized MYC window itself (cos 0.9998), i.e. a similarity search.
- **Scan.** 19,997 reviewed human proteins.
  - Window = span length (55), stride 1, max over windows. This matches `query_proteins_by_cav.py --sliding-window`.
  - Output: `converged_C1e-4/myc_conv_scores.tsv`.
- **Labels.** `converged_C1e-4/ebox_label/`: HOCOMOCO v12 core plus JASPAR 2024, a literature check, and inference for proteins with no database motif (label set db+lit+inf).
  - The literature check is not blind: the curated proteins came from top-20 lists.
- **Export.** `myc_mad/pca_sweep_supp/export_fig5BC_myc_data.py` → `figure_data/fig5BC_myc_top20.csv`. The query (MYC) is excluded from the ranking.
- **Top 20.** All bHLH; 16 are CACGTG binders.
  - The non-binders are TFAP4 (CAGCTG), SIM1 and SIM2.
  - MYCLP1 (a pseudogene) has no motif data.
- **Choice of k.** Within-bHLH CACGTG AUROC is ≥ 0.80 for k = 96–160 and 0.53 without PCA. Describe this as a plateau, not a tuned optimum.
  - Supplement `fig_supp_pca_sweep`: `pca_sweep_supp/export_pca_sweep_cacgtg_data.py` → `pca_sweep_cacgtg.csv`.

### 5.2 Panels C–E — SERK vs non-SERK LRR-RK co-receptor interfaces
Working directory: `lrr_serkcik/`.

- **Interface region.** For each receptor, the 150 residues ending 20 before the transmembrane helix (`serkcik.span`, `homolog_usable.json`).
  - Embed the full-length protein, then mean-pool the window (in-context). The CAV and the cosine baseline use the same embeddings.
- **Training set.** 32 Arabidopsis receptors.
  - 14 SERK-binding: BRI1, FLS2, EFR, PSKR1, RLK5, HSL1, HSL2, TDR, RLP23, PEPR1, PEPR2, ERECTA, GSO1, GSO2.
  - 18 negatives: CLV1, RPK2, BAM1–3, CLV2, PRK1–6, SRF2, SRF4–8.
  - BAK1, SERK4 and SERK5 are dropped from the original 35 (asserted in `export_fig5E_data.py`).
- **Test set.** 36 non-Arabidopsis receptors (24 SERK, 12 non-SERK): 8 curated (rice, tomato, maize, canola) and 28 found by homology to the training genes.
- **CAV.** Balanced logistic regression, C = 1, with a `StandardScaler` fit on the 32 training interfaces.
  - At C = 10⁻⁴ (C·n ≈ 0.003) the CAV shrinks to the class-mean difference and AUROC falls to 0.53.
- **Panel D.** All-vs-all BLASTp over the 150-aa windows (`make_hardneg_blast_heatmap.py`), exported to `figure_data/fig6_heatmap.csv`. It shows selected hard cases, not a failure rate.
- **Panel E.** `export_fig5E_data.py` → `figure_data/fig6_{densities,auroc}.csv`.
  - Signed cosine 1-NN margin = cosine to the nearest SERK training interface minus cosine to the nearest non-SERK one.
  - AUROC: cosine margin 0.986, CAV 1.000. All the discordant pairs involve one gene, P0DL10 (maize TD1).
  - BLAST 1-NN also reaches 1.000 on this split, so the claim is that the CAV sharpens the separation, not that homology fails (`BASELINE_FRAMING.md`).

---

## 6. Figure 6 — single-cell transition axes

All analysis code is in `single_cell/official_tokenizer/` (see its README). Embeddings and runs are in `/xdisk/clairemcwhite/clairemcwhite/sc_geneformer_official/`.

### 6.1 Data and embedding
- **Atlas.** CELLxGENE dataset `b617ee1b-f8c8-4de9-b82b-e803ab93550d` (Guimarães et al. 2024).
  - Downloaded and split into concepts by `single_cell/01_download_data.sh` → `build_cav_library_agent.py`. This step is LLM-assisted metadata inspection; the prompt is the `--prompt` argument in the driver.
  - Concepts are cell type × tissue × disease state, with ≤ 500 cells, groups under 50 cells skipped, and blood/PBMC/lymph node excluded. 93 CAVs in total.
- **Background.** 10,000 cells drawn uniformly (seed 42) from the full human Census (`census_version="stable"`), with no filtering (`single_cell/00_get_negatives.sh`; join IDs in `single_cell/reference_population/`).
- **Embedding.** Geneformer-V2-104M, official V2 tokenization reimplemented in `01_embed_official.py`:
  - Ensembl IDs → `ensembl_mapping_dict_gc104M`, with duplicates summed;
  - counts / n_counts × 10⁴ / gene non-zero median;
  - top 4,094 genes, framed by `<cls>` … `<eos>`;
  - readout `cls_penult`, the penultimate-layer `<cls>` state, 768-d.
  - Ran on CPU (`10_embed_cpu.sbatch`).

### 6.2 CAVs and the transition axis
- **CAVs.** `02_train_and_coords.py`.
  - Shared `StandardScaler` and PCA(128, seed 42), both fit on the 10k reference cells.
  - Balanced logistic regression, C = 1, lbfgs, `tol=1e-10`, retrained on the stored training cells (`00_recover_train_ids.py`).
  - The recipe is the conv1 arm of `single_cell/convergence_check/01_refit_cavs.py`. On the original embeddings it reproduces `convergence_check/coords_conv1.tsv` (r ≥ 0.9999999).
- **Lung 10x restriction.** The lung pairs use 10x 3′ v2 cells only (`single_cell/scripts/train_10x_lung_cavs.py`, `extend_hierarchy_lung10x.py`). "lung10x" is a naming convention, not an atlas tissue value.
- **Axis.** `16a_build_diff_coords.py` computes axis = unit(unit(tumour CAV) − unit(matched normal CAV)). The cell score is the projection onto it (`runs/cls_penult/coords_diff.tsv`; the column keeps its old `L2__…` name).

### 6.3 Gene correlation and pseudobulk DE (panels C, D)
- **Paired donors.** Both analyses use only donors that contribute both conditions. Donor-blind pooling is open to a Simpson's-paradox confound.
- **Gene correlation.** Pearson r between log1p(CP10K) expression and the cell score, over paired-donor cells (`single_cell/scripts/paired_donor_gene_corr.py`; outputs `runs/cls_penult/gene_corr_allpairs_diff/`).
- **Pseudobulk DE** (`single_cell/scripts/paired_donor_pseudobulk.py` + `paired_donor_pseudobulk_de.R`; all-pairs copy `single_cell/convergence_check/13a_pseudobulk_allpairs.py`, `13b_pseudobulk_de_allpairs.R` → `convergence_check/de_pseudobulk_allpairs/min10/`):
  - raw UMIs summed per donor × condition;
  - samples under 10 cells dropped, together with the donor's paired sample;
  - `edgeR::filterByExpr(group = condition)`;
  - DESeq2 1.46.0 with `~ donor + condition` and normal as reference, Wald test, BH correction;
  - the plotted value is the **unshrunk MLE** log2FC (column `log2fc_mle`);
  - environment: `/xdisk/clairemcwhite/clairemcwhite/envs/pseudobulk_de` (R 4.4.3, Bioconductor 3.20).
- **Comparison statistic.** Spearman between log2FC and gene r over genes tested in both:

  | Pair | Genes | Donors | ρ |
  |---|---|---|---|
  | neutrophil / breast | 1,326 | 2 | 0.77 |
  | epithelial / lung (10x) | 1,446 | 4 | 0.80 |
  | neutrophil / lung (rendered; not in text) | 1,846 | 5 | 0.75 |

  - Continuum genes removed by `filterByExpr` (AKR1C2, FGB, S100P) appear in the strips but not in the scatter.
  - These tumour-heterogeneous genes are expressed in only 2–3 of the 4 donors.
- **Export and render.** `single_cell/08_export_fig5_data.sh` → `scripts/export_fig5_case_studies.py` (which imports `cav_continuum_viz.py`) → `figure_data/sc_{pair_meta,de_vs_cav,continuum_cells,continuum_genes}.csv`, then the Fig 6 block of `figures.R`.
  - Strips show a 1,000-cell plotting subsample. All statistics use every cell.

### 6.4 Regularization robustness (text, no panel)
- **Gene-ranking Spearman, C = 1 vs C = 10⁻⁴.** `single_cell/official_tokenizer/20_compare_C.py` → `compare_C.txt` (column `gene_r_spearman`). Breast neutrophils 0.97, lung epithelial 0.91.
- **Top-40 overlap** (top 20 + bottom 20 by r). `25_C_modules_breastneut_lungepi.py` → `figures/previews/sc_C_modules_breastneut_lungepi_{genes,metrics}.tsv`. Shared genes: 33 (breast) and 27 (lung).
- **Arms compared.** `runs/cls_penult` vs `runs/cls_penult_C1e-4`.

---

## 7. Supplementary figures rendered but not cited in `text.txt`
- **`fig_supp_pca_sweep`.** MYC and MAD1 (Q05195 61–94) single-positive CAVs at k = 8–512, plus no PCA. C = 10⁻⁴, judged by within-bHLH CACGTG AUROC. See §5.1. This supports the Discussion's sentence on PCA dimensionality.
- **`fig_supp_projection_comparison`.** DE-vs-CAV ρ for four axes: raw cancer CAV, minus L0, L2, and cancer − normal (`single_cell/10_export_projection_comparison_supp_data.sh` → `proj_de_vs_cav.csv`). Only the difference axis works for all three pairs (ρ 0.77 / 0.80 / 0.75).
- **`fig_ec_eval`.** Per-term EC recall vs EC-Bench tools (`ec_eval/02`–`04`). Superseded by Fig 3A.

---

## 8. Software and rendering
- **R.** R 4.5.3 with tidyverse (ggplot2 4.0.3), cowplot 1.2.0, ggridges 0.5.7, ggrepel 0.9.8 and ggrastr 1.0.2. The circle pack also uses ggraph 2.2.2, tidygraph 1.3.1 and igraph 2.3.3 (`~/R/4.5`).
- **Python.** Conda env `core_pkgs4`: torch 2.7.1, transformers 4.57.6, scikit-learn 1.8.0, scipy 1.17.1, numpy 2.4.3, pandas 2.3.3, scanpy 1.11.5, umap-learn 0.5.11.
- **Fonts.** Arial is requested, but **Liberation Sans** (metric-compatible) is what renders.
- **Graphics.**
  - Okabe–Ito palette.
  - Dense scatters are rasterized with ggrastr at 150 dpi.
  - PDFs are written with `cairo_pdf`; PNGs at 300 dpi (600 dpi for Fig 4).
- **Data release.** Four Zenodo records are planned (background embeddings, GO, EC and Pfam libraries); see `distribution/DISTRIBUTION_PLAN.md`. Per-motif Pfam thresholds come from `distribution/pfam_thresholds/build_motif_thresholds.py`.

---

## 9. Discrepancies between the code and `text.txt` (resolve before submission)
1. **Fig 6 cases.**
   - `figures.R` renders three cases: neutrophil/breast, epithelial/lung (10x) and neutrophil/lung. The text discusses only the first two.
   - The Fig 6 legend's case letters say "(B)/(C)" inside "(C–D)".
   - The text's Methods (Single Cell Analysis, last paragraph) still name "colorectal fibroblasts" as a highlighted case.
2. **Fig 4 libraries.**
   - Panel A uses the converged C = 10⁻⁴ Pfam library.
   - Panel C (Plexin-A1) was scored with the original early-stopped, per-motif-scaler 27k library (`shamail/…/cavs`).
   - Panel B uses its own 36-layer CAVs.
   - The text describes a single library.
   - `cav_fig4` still carries the obsolete 3-D UMAP as panel A.
3. **Zn screen count.** The text says 84 families. The only recorded run, the skill's reference result, has 82 reproducible families on the earlier 12-cut clustering. No Supplemental File 2 or D–J outputs are in this repository. Locate them, or re-run the screen on the §4.1 clustering.
4. **Interactive circle pack.** ~~mcwhitelab.org/pfam-circlepack is not built from the final clustering.~~ Resolved 2026-10-09: `average_mainfrac_045-005_promote2/web/index.html` is built from the final layout (§4.1). The site still needs the new file uploaded. The old R builder is archived.
5. **SERK negatives.** The text calls CLV2 "CIK-dependent" (six). The analysis notes count five CIK founders (CLV1, RPK2, BAM1–3), with CLV2 among the hard negatives.
6. **Small method details missing from the text:**
   - the DIAMOND family edge also needs ≥ 50% coverage;
   - the EC exclusion also removes proFAB **test** positives and UniProt fragments;
   - EC 3.4.22.28 positives were truncated to 2,046 aa.
7. **Fig 3B rule.** The class-3.1 protein (P0DUB9) is not the first alphabetical leak-free accession (A0A345BJN2 is). Change the pick or the stated rule.
8. **Homology baseline settings (2026-10-09).**
   - The Fig 2B / 3A baseline is now DIAMOND `--ultra-sensitive`, e ≤ 1 (§2.2). The text's Methods (line 151) still says "more-sensitive mode with an E-value threshold of 10⁻³".
   - Suggested replacement: "We also included a sequence-similarity baseline using DIAMOND (v2.0.9) (Buchfink et al. 2021) in ultra-sensitive mode with an E-value threshold of 1, as recommended for annotation transfer by [CITATION], and self-hits excluded."
   - The citation for the ultra-sensitive / e ≤ 1 recommendation still has to be added.
   - The legends of Fig 2B and 3A label this baseline "BLAST" although it is a DIAMOND search; decide whether to relabel.
   - The family count sentence (line 152) refers to the separate all-vs-all DIAMOND run, which is unchanged.
9. **Untracked live code.**
   - `single_cell/official_tokenizer/` is untracked.
   - Live Fig 5 code under `interface_search/fastas/pairs/{myc_mad,lrr_serkcik}/` is gitignored (`.py`, `.sh`, `COMPARISON.md`, `BASELINE_FRAMING.md`).
   - So are the circle-pack scripts under `motif_clustering/results/…/diffusion_linkages/` and the render template.
   - Some numbers in the text come from gitignored outputs, e.g. `figures/previews/sc_C_modules_breastneut_lungepi_genes.tsv` (33/27) and `single_cell/official_tokenizer/compare_C.txt` (0.97/0.91).
   - All of these need `git add -f` (or narrower ignore rules) before the repository is public.
