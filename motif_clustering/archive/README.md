# Archived motif_clustering drivers

`01_visualize_umap.sh`, `02_dendrogram.sh` — archived 2026-09-16.

Both resolve `$AHMAD_CAV_OUTPUTS_DIR`
(`/groups/clairemcwhite/ahmad_workspace/esm_c/tcav_outputs_esmplusplus_all`), the
**7,692-motif** Pfam CAV library. The paper now uses the **27,472-motif** library at
`/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all`, which is what
`03_dendrogram_comparison.sh`, `cav_cluster.py`, `context_leakage/`, and
`figures/figure_data/vav_motif_repro/` all read, and it is the library being released
(see `distribution/DISTRIBUTION_PLAN.md`, B2).

Their outputs — `results/figures/pfam_umap*.html`, `pfam_umap*_3d.gif`,
`results/figures/pfam_clusters.csv` — are superseded by `pfam_clusters_27k.csv` and are
not consumed by any current figure script. They are left in `results/` untouched; delete
them only if you are sure nothing references them.

To revive either script against the 27k library, replace `$AHMAD_CAV_OUTPUTS_DIR/cavs`
with `/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs`.
