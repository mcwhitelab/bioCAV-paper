# Archived single-cell pipeline code (2026-09-13)

Nothing here is read by the paper figures any more. Kept for provenance; all
moved with `git mv` (except `mixedlm_scaleup/`, which was untracked), so history
is intact. Nothing was deleted.

| Archived | Why | Live replacement |
|---|---|---|
| `alt_00_diffe.sh` | Pseudobulk DESeq2 run (`--method deseq2`) over 25 pairs → `results/de_vs_cav/`. No figure ever read that directory; it failed outright on `neutrophil__breast` (3 donors). Its only surviving role is as the evidence behind the Methods sentence explaining why pseudobulk was rejected. | `scripts/paired_donor_de_mixedlm.py` |
| `mixedlm_scaleup/` | 29-task SLURM array that ran the mixed model atlas-wide → `results/de_mixedlm/`. Superseded for the figure pairs by the paired-donor re-fit; the remaining 22 pairs are not shown in any figure. | `scripts/paired_donor_de_mixedlm.py` |
| `07_paper_cases.sh` | Exploratory matplotlib case-study PNGs. Built around the skin/melanoma pairs, which were dropped as unfixably assay-confounded. | `08_export_fig5_data.sh` → native ggplot panels |
| `scripts/export_fig5_data.py` | First figure-data exporter; donor-blind (pooled across all donors), so vulnerable to the Simpson's-paradox confound. | `scripts/export_fig5_case_studies.py` |
| `scripts/export_lung10x_fig5_data.py` | Separate lung10x exporter; folded into the single consistent path. | `scripts/export_fig5_case_studies.py` |
| `scripts/draft_adamdec1_module_check.py` | One-off investigation of the ADAMDEC1 / quiescent-fibroblast module. Its conclusion is now a supplemental figure. | `scripts/export_subpop_supp_data.py` |

Data directories left in place but no longer read by any figure:
`results/de_vs_cav/` (DESeq2) and the non-case-study pairs of `results/de_mixedlm/`.
Figure 6 and its supplementals read `results/de_mixedlm_paired_donors/`.

Note: `scripts/lung10x_de_and_corr.py` stays live — its DE half is superseded,
but it is still the only producer of the lung10x row in
`results/gene_correlation_celltype_tissue/`, which `06_viz_continuum.sh` reads.
See its docstring.
