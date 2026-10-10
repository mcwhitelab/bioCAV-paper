# bioCAV-paper

Analysis code and figure data for:

> **Representing biological concepts as directions in foundation model embeddings.**
> Ahmad Shamail and Claire D. McWhite, The University of Arizona.

A concept activation vector (CAV) is the unit weight vector of a linear classifier that
separates examples of a concept from a common background in a pretrained embedding space.
This repository reproduces the paper's analyses of:
- protein CAVs (ESM-C 600M) for GO, EC, Pfam and binding interfaces;
- single-cell CAVs (Geneformer V2).

The CAV method itself (training, scoring, localization) lives in a separate package:
**https://github.com/mcwhitelab/bioCAV**.

## Layout

| Directory | Contents |
|---|---|
| `config/paths.sh` | Every external path and environment setting. All drivers source it. Edit this first |
| `shared_protein/` | proFAB download and formatting, the 10k UniRef50 background, shared scaler/PCA, CAV training commands |
| `distribution/` | Final CAV libraries (GO/EC at C = 10⁻⁴; converged Pfam), Pfam score thresholds, Zenodo release plan |
| `go_eval/` | Fig 2: GO temporal validation, DeepGO-SE and DIAMOND baselines, family-count bins, ranks |
| `ec_eval/` | Fig 3A: EC-Bench evaluation, DIAMOND baseline, CLEAN-SupCon comparison |
| `hierarchy_decay/` | Fig 3B: CAV scores along the EC hierarchy |
| `go_ec_crossval/` | Fig 3D–E: GO–EC CAV cosine agreement against ec2go |
| `motif_clustering/` | Fig 3C combined UMAP; Fig 4A Pfam CAV clustering and circle pack (Supplemental File 1) |
| `interface_search/` | Fig 5: single-example MYC CAV and contrastive SERK/LRR-RK CAV |
| `single_cell/` | Fig 6: Geneformer embeddings, single-cell CAVs, transition axes, pseudobulk DE |
| `figures/` | `figures.R` renders every data panel from flat files in `figures/figure_data/` |
| `METHODS_FIGURE_STEPS.md` | Extended methods: the scripts, parameters and numbers behind each panel |

## Figure → code map

| Panel(s) | Driver / exporter | Data read by `figures.R` |
|---|---|---|
| 2C | `figures/rerun_fig23_C1e-4.sh` → `go_eval/rank_C1e-4/rerank_go.py` | `go_rank_families_{mf,bp,cc}.csv` |
| 2D–E | `go_eval/06_eval_exclude_known_positives.sh`; re-run on the C = 10⁻⁴ library by `figures/rerun_fig23_C1e-4.sh` | `temporal_tool_comparison_{mf,bp,cc}.csv` |
| 3A | `ec_eval/rank_C1e-4/rerank_ec.py` | `ec_ranks_samelevel_noknown.csv` |
| 3B | `hierarchy_decay/00_run.sh` | `ec_hierarchy_decay.csv` |
| 3C | `go_eval/05_visualize_umap.sh` → `motif_clustering/combined_cav_umap_3d.py` | `go_ec_combined_umap.csv` |
| 3D–E | `go_ec_crossval/00_run_validation.sh` | `go_ec_cosine_{pairs,background_sample}.csv` |
| 4A, Supp File 1 | `motif_clustering/cav_cluster.py`, `promote_largest_clusters.py`, `scramble_cluster_ids.py` | (rendered separately with ggraph) |
| 4B–C | `figures/figure_data/vav_motif_repro/score_{vav_all_layers,bigreceptors}.py` | `vav_motif_repro/*.json` |
| 5B | `interface_search/fastas/pairs/myc_mad/pca_sweep_supp/export_fig5BC_myc_data.py` | `fig5BC_myc_top20.csv` |
| 5D–E | `interface_search/fastas/pairs/lrr_serkcik/export_fig5E_data.py` | `fig6_{heatmap,densities,auroc}.csv` |
| 6C–D | `single_cell/official_tokenizer/` → `single_cell/08_export_fig5_data.sh` | `sc_*.csv` |

Panels 1, 2A–B, 4D–J, 5A, 5C and 6A–B are schematics or structural renderings made outside this pipeline.
Some file names carry older figure numbers (e.g. the `fig6_*` files feed Fig 5). The first
table of `METHODS_FIGURE_STEPS.md` maps old names to final panels.

## Redrawing the figures

```bash
# 1. point config/paths.sh at your copies of the external resources
# 2. render every panel from the shipped flat files (no GPU or model needed)
bash figures/00_make_figures.sh        # writes figures/figures/*.pdf|png
```

- `figures.R` needs R ≥ 4.4 with tidyverse, cowplot, ggridges, ggrepel and ggrastr.
- Each figure block runs only if its input files exist. Check the console for `Skipping …` lines.

Regenerating `figure_data/` from scratch needs:
- the ESM-C 600M weights (ESMplusplus_large);
- Geneformer-V2-104M;
- the proFAB, QuickGO, EC-Bench and CELLxGENE downloads;
- the CAV libraries;
- a SLURM cluster with GPUs for embedding.

CAV training itself runs on CPU. The numbered `NN_*.sh` scripts in each directory run in order.

## Data

The CAV libraries and the shared background embeddings will be released as four Zenodo
records: background, GO (6,678), EC (953) and Pfam (27,472). See
`distribution/DISTRIBUTION_PLAN.md`.

Large intermediates (embeddings, `results/`, `static/` downloads) are not tracked in git.

## Archived material

Superseded analyses, drafts and old figure data have been moved under `archive/` folders,
which are gitignored and kept only on the authors' system for provenance. See
`archive/2026-10-09_cleanup/MANIFEST.md` locally.

## License

MIT (see `LICENSE`).
