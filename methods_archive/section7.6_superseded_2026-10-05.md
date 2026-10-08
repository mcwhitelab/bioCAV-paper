# Superseded METHODS_FIGURE_STEPS.md §7.6 (archived 2026-10-05)

Replaced by the C = 10⁻⁴ diffusion_match4_nopromote circle pack. Kept for provenance: the old
early-stopped library's circle pack (steps 14–25), the t = 4 sensitivity table, and the
2026-10-05 cosine/diffusion candidate notes.

### 7.6 Hierarchical circle pack of the Pfam motif CAV library (diffusion distance, t = 16)

**Methods prose (C = 10⁻⁴ shared-scaler library, decided 2026-09-29; drafted 2026-10-05):**
> The 27,472 unit-normalised Pfam CAVs were connected in a mutual 15-nearest-neighbour graph by
> cosine similarity (edges below 0.10 removed), and clustered hierarchically by Ward linkage on
> truncated diffusion-map coordinates approximating the random-walk (walktrap) distance at walk
> length t = 16. The tree was cut at 314, 2,769, 4,629, 9,172 and 18,984 clusters. The three
> largest top-level clusters, holding 59% of families, were replaced by their sub-clusters at the
> next cut, and the resulting five-level hierarchy was drawn as a circle pack with ggraph,
> colouring each family by its Pfam 38.2 clan.

**SUPERSEDING CANDIDATE (2026-10-05): cosine/average tree at fixed thresholds.** On the
C = 10⁻⁴ library the plain cosine tree scores better than diffusion against Pfam 38.2 clans:
ARI 0.362 at 0.75 (1,327 clusters) vs 0.212; best AMI 0.701 vs 0.657; best clan F1 0.762 vs
0.667. Its contiguity is slightly lower (51.6 vs 64.4).

Cutting the diffusion tree at the counts these thresholds give fails: at 60 clusters one cluster
holds 27,308 families. Files: `motif_clustering/results/converged/C0.0001_sharedscaler/threshold_cuts/`
(`cosine_average.csv` / `.scrambled.csv`, `cosine_circlepack/` with top 3 promoted,
`cosine_circlepack_top4/`).

Draft prose (fill in N = number of clusters promoted):
> The 27,472 unit-normalised Pfam CAVs were clustered hierarchically by average linkage on cosine
> distance (1 − cosine similarity), and the tree was cut at fixed cosine distances of 0.95, 0.85,
> 0.80, 0.70 and 0.50, giving 60, 466, 825, 1,955 and 6,032 clusters. The N largest top-level
> clusters were replaced by their sub-clusters at the next cut before layout, so that no single
> cluster dominated the figure. The resulting five-level hierarchy was drawn as a circle pack
> with ggraph, colouring each family by its Pfam 38.2 clan.

**User preference (2026-10-05):** the circle pack the user prefers is the DIFFUSION tree cut at
the cosine-threshold counts: `threshold_cuts/pfam_circlepack_clanid_diffusion_threshcounts_k00060.{png,pdf}`,
from `threshold_cuts/clusters_promoted_top3.csv`.
- The tree is the t = 16 Ward tree, cut at 60 / 466 / 825 / 1,955 / 6,032 clusters, the counts
  that cosine distances 0.95 / 0.85 / 0.80 / 0.70 / 0.50 give on this library.
- Because the 60-cluster level holds one cluster of 27,308 families, promoting the top 3 moves
  nearly every family up one level. Families per level after promotion: 455 / 825 / 1,955 /
  6,032 / 6,032. So the drawn hierarchy is effectively four levels: the 466, 825, 1,955 and
  6,032 cuts. The largest top-level circle holds 4,264 families.
- The methods text must describe it this way, not as a five-level cut at 60 clusters.

Caveat for the legend (cosine version): average linkage chains a large core. Even with the top 3 promoted (78% of
families), one top-level circle holds 8,627 families, mostly without a clan.

Source (diffusion version above): `motif_clustering/results/converged/C0.0001_sharedscaler/`
- the clustering (`clusters_diffusion.csv` and `.scrambled.csv`);
- the promoted table (`clusters_promoted_top3.csv`: 16,258 families promoted out of clusters of
  11,613, 2,478 and 2,167; 1,859 top-level clusters after promotion);
- the render (`pfam_circlepack_clanid_k00314.{png,pdf}`).

Steps 14–25 below still describe the OLD early-stopped library's circle pack. Its pipeline
(steps 15–24) is unchanged; only the library (step 14), the counts and the file paths differ.

Working directory: `motif_clustering/`. The figure is
`results/retrained_shared_scaler/promoted_top3/old_diffusion/pfam_circlepack_clanid_k00314.{pdf,png}`.
It is a PNG at 300 dpi plus a vector `cairo_pdf`, 11.4 × 13.2 in. Placement in the
numbered set is `[TK]`.

14. **Library: the original independent-scaler Pfam CAVs, not the shared-scaler retrain.**
    - 27,459 layer-25 CAVs, `L25_concept_v1.npy`, from the 27,472-motif library at
      **`/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs`**. The 13 missing
      concepts are the training failures (§1.7).
    - Each CAV was trained in its own per-concept `StandardScaler` frame.
    - The shared-scaler retrain under `/xdisk/clairemcwhite/clairemcwhite/pfam_cavs_shared_scaler`
      has its own circle packs (`retrained_diffusion/`, `retrained_shared_scaler/*.png`). Those
      are **not** the figure.
    - The folder name `retrained_shared_scaler/` is only where the promoted renders were
      written. The library inside `old_diffusion/` is the old one.
15. **Direction matrix.** Each weight vector is L2-normalised to unit length, so the inputs
    are 27,459 directions in 1,152 dimensions (`cav_cluster.py::load_directions`).
16. **Mutual-kNN cosine graph** (`cav_cluster.py::knn_graph`).
    - Each concept keeps its k = 15 nearest neighbours by cosine similarity.
    - Neighbours below a cosine of **0.10** are dropped.
    - An edge is kept only if it is **mutual**, i.e. each concept is in the other's top 15.
      The edge weight is the cosine similarity.
    - Concepts left with no edge (1,304) are reattached to their single nearest neighbour.
    - Result: 88,946 undirected edges, mean degree 6.5, 186 connected components. The largest
      component holds 26,881 concepts.
17. **Truncated diffusion-map coordinates** (`cav_cluster.py::diffusion_coords`).
    - The target is walktrap's random-walk distance (Pons & Latapy 2006) at walk length
      **t = 16**:
      r_ij = √(Σ_k (P^t_ik − P^t_jk)² / d_k), with P = D⁻¹W.
    - This is computed spectrally, never forming P^t:
      - take the top eigenpairs of the symmetric S = D^−½ W D^−½;
      - drop the trivial eigenvector;
      - scale each remaining eigenvector by λ^t, then by D^−½ row-wise.
    - Only the top **286 of 27,459** eigenpairs are used (`scipy.sparse.linalg.eigsh`), so
      Euclidean distance between rows *approximates* r_ij. It equals r_ij only with the full
      spectrum. (The `cav_cluster.py` docstring says "equal"; do not copy that.)
    - The truncation is not negligible. The smallest kept eigenvalue is 0.972
      (0.972^16 ≈ 0.63), so the discarded eigenvalues just below it still carry substantial
      weight at t = 16.
    - The request of 286 is `--n-eigs 100` widened by the 186 connected components. 107
      eigenvalues exceed 0.999 (component indicators, kept so that separate components sit
      far apart). After removing the trivial vector, about 179 within-component vectors
      remain.
    - Output: 285-dimensional coordinates.
18. **Ward linkage** on the 285-d diffusion coordinates (`scipy.cluster.hierarchy.linkage`,
    `method="ward"`). Merge heights span 2.5e-17 to 3.77.
19. **Cuts.** The tree is cut with `fcluster(criterion="maxclust")` at **314, 2,769, 4,629,
    9,172 and 18,984** clusters. The corresponding Ward heights are 0.449, 0.0130, 0.00729,
    0.00299 and 0.000634.
    - These counts are not tuned on this tree. They are the cluster counts that the historical
      cosine/average-linkage tree gives at cosine-distance thresholds of 0.95, 0.85, 0.80,
      0.70 and 0.50 (`DEFAULT_MAXCLUST`).
    - Reusing them keeps every dendrogram construction comparable at matched resolution.
20. **Provenance: this run has no driver script. State the parameters explicitly.**
    - `03_dendrogram_comparison.sh` runs only `diffusion_ward` with the defaults (floor 0.15,
      200 eigenvectors, t = 4).
    - The figure's tree is the parameter-sweep run
      `results/dendrogram_comparison/diffusion_f0.10_e100_t16.*`, logged in `full_diff.log`.
      It is equivalent to:
      `cav_cluster.py --distance diffusion --linkage ward --knn 15 --knn-floor 0.10 --n-eigs 100 --walk-steps 16 --cut-mode maxclust --no-shuffle-ids`.
    - `results/figures/circlepack_input_diffusion.csv` re-cuts that tree at the five levels
      above. Its `.linkage.npy` is byte-identical (md5 `1fb921a1…`) to the sweep run's.
    - The command above is **reconstructed** from the run's filename and `full_diff.log`,
      which logs k, floor, mutual, the eigenvector request and t. No record of the literal
      invocation exists.
21. **Clan labels** were relabelled from Pfam 38.1 to **Pfam 38.2**
    (`update_pfam_clans_38_2.py`).
    - Only the clan columns change. The tree and cuts are untouched.
    - After relabelling, 13,164 of the 27,459 families belong to 891 clans.
22. **Promoting the three largest top-level clusters, in the data, before layout**
    (`promote_largest_clusters.py --top-n 3`).
    - At the 314-cluster cut, the three largest clusters hold 9,737, 5,613 and 1,635
      families: 16,985 families, 61.9% of the library.
    - For every family in those clusters, each cut column takes the value of the next finer
      cut, and the finest level keeps its own. Shifted IDs are offset by 10⁶ so they cannot
      collide with existing IDs.
    - The three clusters therefore never exist in the hierarchy the layout sees. Their
      sub-clusters become top-level circles.
    - Distinct clusters per level after promotion: 1,752 / 3,928 / 7,608 / 15,574 / 18,984.
    - Output: `old_diffusion/clusters_promoted_top3.csv`.
    - This replaces the earlier approach of hiding the largest outlines at draw time
      (`hide_top_n`). The rendering script runs with `hide_top_n <- 0`.
    - The legend must state that the outermost circles are not the 314-cluster cut
      throughout: for about 62% of families they are the 2,769-cluster cut.
23. **Circle-pack layout** (`old_diffusion/pfam_circlepack_promoted.R`, R `ggraph`).
    - The tree is origin → five cut levels → individual Pfam families. Cluster IDs are made
      unique per level by suffixing the column name.
    - The tree is built as a `tidygraph` edge list and packed with
      `create_layout(layout = "circlepack")`: 75,306 nodes, depth 6.
    - All leaves have equal weight, so a cluster circle's area scales roughly with its family
      count (packing gaps inflate it).
    - The layout is cached as `circlepack_layout_k00314.rds`.
24. **Rendering: clan-identity colouring only.**
    - Cluster circles are drawn as outlines only: fill white, stroke `#cbcac3`, line width
      thinning with depth from 0.35 to 0.04.
    - Each family is a filled dot with no stroke. Its colour is one random colour per Pfam
      clan, drawn from a fixed per-clan table, `clan_colours_R.csv`, so a clan keeps its
      colour across variants and Pfam releases.
    - The table was pre-built on 2026-09-05. Clans new to a release get hues sampled with
      `set.seed(42 + <rows already in table>)`, chroma 70–120 and luminance 35–65. State only
      that colours are random and fixed per clan; the original seed is not verified.
    - Families with no clan assignment are grey, `#c9c8c0`.
    - The colours are not meant to identify individual clans. They show which branches
      share a clan.
    - There is no legend.
25. **Sensitivity to walk length (t = 4).** The figure uses t = 16. No record shows t = 16 was
    deliberately selected over t = 4, so present this as a sensitivity check, not a choice.
    Both were scored against Pfam 38.1 clans (`score_dendrogram.py`;
    `results/dendrogram_comparison/*.scores.csv`). They are essentially tied:

    | Run | Best ARI (cut) | Contiguity ratio | Median clan span |
    |---|---|---|---|
    | diffusion, floor 0.10, t = 16 (**figure**) | 0.165 (1,326) | 11.7 | 233 |
    | diffusion, floor 0.10, t = 4 | 0.162 (1,326) | 11.0 | 276 |
    | diffusion, floor 0.15, t = 4 (`diffusion_ward`) | 0.158 (1,326) | 10.0 | 236 |
    | cosine / average | 0.207 (2,769) | 8.0 | 351 |
    | cosine / Ward | 0.166 (314) | 5.1 | 444 |

    - The contiguity ratio is the median clan span under a permuted leaf order divided by the
      observed median span. Higher means clans are more contiguous.
    - Diffusion's advantage is leaf-order contiguity, which is what the circle pack displays.
      It is not flat-cut ARI, where cosine/average is higher.
    - Do not describe t = 16 as optimised. The difference from t = 4 is too small to
      justify either.
