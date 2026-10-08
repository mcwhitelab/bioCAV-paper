# Distribution plan — background embeddings + GO / EC / Pfam CAV libraries

Draft 2026-09-15. Four Zenodo records, cross-linked, each self-sufficient for scoring.

---

## Blocking issues found while surveying (resolve before packaging)

### B1. GO CAVs were trained in one frame and are scored in another
`02_make_cav_commands.sh` emits `train_cav_from_span.sh $f` with **no `--scaler-pkl`**, so
`train_cav_from_embeddings.py` fit a *fresh* StandardScaler on (positives + 10k negatives)
for each concept. A later rerun (`cav_GO_COMMANDS_remaining.sh`) *does* pass
`--scaler-pkl reference_population/scaler_v1.pkl`, so the GO library is a mix of both.

Verified by comparing `mean_`/`scale_` arrays (md5 is not a valid test — joblib writes
differing bytes for identical content), on random samples:

| Library | trained with own scaler | trained with shared scaler | scored with |
|---|---|---|---|
| **GO** (6,677) | 49/80 (built 2026-05-21/22) | 31/80 (built 2026-05-25 → 06-10) | **shared** (`eval_go_temporal.py` default `--scaler-pkl`, `03_compare_to_outside_tools.sh`) |
| **EC** (952) | 0/60 | 60/60 | shared (`02_eval_ec.sh --scaler-pkl`) — **consistent** |
| **Pfam** (27,472) | 60/60 | 0/60 | **own** (`score_candidates_fig4b.py` loads `CAV_DIR/<motif>/L25_scaler_v1.pkl`) — **consistent** |

So the defect is **GO only**: ~61% of GO CAVs have weights fit in a per-concept frame but are
standardised at scoring time with the shared frame. EC and Pfam are each internally consistent.

Magnitude: for the own-frame GO CAVs, `scale_` differs from shared by median 2.1e-2
(max 1.3e-1) relative to the median scale. Small, but it is a systematic train/score mismatch
affecting roughly 4,000 of the 6,677 released GO CAVs.

Options:
- (a) **Retrain the affected GO CAVs with `--scaler-pkl`.** Training itself is CPU-cheap, but
  `train_cav_from_span.sh` step 5 deletes the pooled positives, so this means re-embedding
  ~4,000 `.span.fasta` files on GPU. Cleanest fix.
- (b) **Ship each GO CAV's own scaler and score with it**, as Pfam already does. Correct by
  construction; gives up cross-CAV GO score comparability, which several figures rely on.
- (c) **Quantify and document**: rescore a held-out set both ways, report rank correlation,
  ship the shared scaler with a stated caveat.

**RESOLVED 2026-09-15: option (a).** Claire chose to retrain. Exact classification of all
6,677 (`retrain_go/classify_scalers.py`, using `n_samples_seen_ > 10000` — exact, since an
own-frame scaler saw 10,000 + n_positive; zero invariant violations): **3,930 own / 2,747
shared**. Retrain submitted as slurm job **23877661** (`retrain_go/01_retrain_go_cavs.sbatch`,
array 1-40%10). Originals backed up to `retrain_go/backup_own_scaler_cavs_20260915.tar.zst`
(228 MB, 3,930 dirs). Verify with `retrain_go/02_verify_retrain.sh` — expect own=0.
Update the bioCAV README's one-coordinate-system claim once it passes.
**These CAVs are not being packaged yet** (Claire, 2026-09-15).

### B2. Pfam library choice — RESOLVED: release the 27,472 library
Claire confirmed: use `/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all`
(**27,472** motifs; `training_summary.json` reports 27,459 successful / 13 failed), not
`$AHMAD_CAV_OUTPUTS_DIR` (7,692). Already used by `figures/figure_data/vav_motif_repro/*.py`,
`motif_clustering/03_dendrogram_comparison.sh`, and `context_leakage/`.

Two consequences to handle:
- ~~`motif_clustering/01_visualize_umap.sh` and `02_dendrogram.sh` still resolve
  `$AHMAD_CAV_OUTPUTS_DIR` (7,692)~~ **RESOLVED 2026-09-16** — both are stale and were moved
  to `motif_clustering/archive/` (see its README). Their outputs are superseded by
  `pfam_clusters_27k.csv`. No script in the paper repo now references the 7,692 library.
- ~~The 27k library has no `motif_thresholds.json`~~ **RESOLVED — no retraining needed.**
  The file was never a separate computation: in `batch_train_cavs.py` it is a pure aggregation
  of fields `train_cav.py` already wrote into every `L{layer}_report_v1.json`. All 27k reports
  carry complete threshold blocks, so it re-aggregates directly.
  `distribution/pfam_thresholds/build_motif_thresholds.py` regenerates it: **27,459 motifs**
  (12 MB), exactly the `training_summary.json` success count; the 13 missing are the 13 failed
  trainings, listed in `motifs_without_thresholds.txt`.

### B3. `hf_embed_new.py` — RESOLVED: tagged `cav` @ `4b32101`
**Which copy the runs used:** `/groups/clairemcwhite/claire_workspace/github/mcwlab_utils/hf_embed_new.py`
(md5 `f733f2cc3cb1`). It was both `HF_EMBED_SCRIPT` in `config/paths.sh` *and* hardcoded as
`embed_script=` in `train_cav_from_span.sh`, so the background set (dataset 1) and every GO and
EC CAV (datasets 2, 3) went through it. Four other copies exist on the cluster — in particular
`/groups/clairemcwhite/mcwlab_utils/` (md5 `f37db72fe365`) is lab-shared and *newer* but is
**not** what produced the data. Cite only the tag below.

**Citation target (2026-09-16):** merged to `main` via PR #4 and released as tag **`cav`** at
commit **`4b32101`**, pushed and reachable on `origin/main`. Pin `mcwhitelab/mcwlab_utils@cav`
in all four records.

**Documented caveat — the tag is not byte-identical to the file that was run.** The merge
(`d27ec49`, main → test_fusion) combined the run version `f733f2cc3cb1` with main's
`f37db72fe365`, so the tagged blob is md5 **`3031537f69be`**. The differences are additive:
a new `maxpool` choice for `--strat`, a new `--aa_trim` option defaulting to `none`
(backward-compatible), and a refactor of the mean-pooling branch.

**Verified output-equivalent, not merely assumed.** Both versions were run on the same input
with the exact flag set the pipeline uses (`--strat mean --get_aa_embeddings
--get_sequence_embedding -l -11 -b 1 --max_length 2048`, no `--aa_trim`):
`sequence_embeddings` and `aa_embeddings` came back **bitwise identical** (`max|diff| = 0.0`).
State this equivalence in each record's `PROVENANCE.md` rather than claiming the tag is the
literal artifact.

**Dataset 4 is the exception: it never used `hf_embed_new.py`.** The Pfam library embeds through
its own `esm2_tcav/src/embed.py` (`process_dataset`), called from `batch_train_cavs.py`. Its
record must cite that code path, not `mcwlab_utils@cav`.

### B4. Redistribution permissions — RESOLVED for the internal sources
- `neg_10000.fasta` (dataset 1) and both Pfam CAV libraries (dataset 4) live in
  `ahmad_workspace` / `/xdisk/clairemcwhite/shamail`. **Ahmad Shamail is a co-author on the
  paper** (confirmed 2026-09-16), so these are co-authored lab outputs, not third-party
  material. List Shamail as a creator on datasets 1 and 4 at minimum; simplest is the same
  creator list across all four records, matching the paper's author list.
- Remaining obligations are **citation, not permission**:
  - **proFAB** — source of the GO/EC positive/negative splits (datasets 2 and 3). We
    redistribute only our derived CAVs and the `.span` files, not proFAB's own feature files.
  - **UniProt / UniRef50** — source of all sequences, including the 10,000 background proteins.
  Put both in each record's README and `CITATION.cff`.
- **CC-BY-4.0 for all four records — settled.** Both upstream dependencies are permissive and
  compatible:
  - **ESM-C / ESM++** — MIT (B5).
  - **proFAB** — MIT (*ProFab Copyright (C) 2022 CanSyL*), which permits redistribution,
    modification and sublicensing. Releasing derived work under CC-BY-4.0 is compatible.
  MIT's one condition is that the copyright and permission notice travel with copies of
  "the Software". What datasets 2 and 3 actually redistribute from proFAB is the `.span`
  accession lists — data, not code — so the condition is arguably not triggered, but
  **include the proFAB MIT notice in both records anyway**: it costs a `THIRD_PARTY.md` and
  removes the ambiguity entirely.
  - **UniProt / UniRef50** — source of all sequences (including the 10,000 background
    proteins); UniProt data is itself CC-BY-4.0, so it is consistent. Cite, don't relicense.

### B5. ESM-C licensing — RESOLVED, no restriction (was: non-commercial)
**Superseded 2026-09-16.** An earlier draft of this plan read
`models/ESMplusplus_large/NOTICE.txt` and concluded all four records were locked to
non-commercial by the EvolutionaryScale Cambrian Non-Commercial License. **That local NOTICE
file is stale.** ESM-C has since been relicensed:

- EvolutionaryScale was absorbed into **Biohub** (CZI Science + CZ Biohub + ES). Biohub now
  publishes **ESMC 300M / 600M / 6B as open weights under MIT**, alongside the 2026 preprint
  *"Language Modeling Materializes a World Model of Protein Biology"* (Candido et al.).
- The legacy `biohub/esmc-600m-2024-12` card carries `license: [mit, other]`. The `other` tag
  resolves to `THIRD_PARTY_NOTICE.md`, which covers only **software dependencies** (PyTorch,
  flash-attn, xformers, einops…) — it places **no restriction on the model weights**.
- Decisive for us: the checkpoint we actually used is Synthyra's ESM++, and
  **`Synthyra/ESMplusplus_large` states `License: mit` and "Checkpoint terms: MIT. The Hub
  model-card identifier is `mit`."** No Cambrian NC reference remains.

**Consequences:** embeddings and CAVs are **not** encumbered. CC-BY-4.0 *is* available for all
four records (see B4). No ESM NOTICE string is required, and the "title must begin with ESM"
constraint does not apply — the four records keep their names.

**One residual action:** `/groups/clairemcwhite/models/ESMplusplus_large/NOTICE.txt` still
carries the old Cambrian NC text, because the checkpoint was downloaded before relicensing.
Refresh the local checkpoint's license files (or drop a `LICENSE_UPDATE.md` beside it) so
on-disk provenance matches reality — otherwise the next person to check, like this plan did,
reaches the wrong conclusion.


---

## Cross-cutting format rules (apply to all four records)

1. **Primary format is `.npz` / TSV, not pickle.** `neg_10000.pkl` and every `scaler_v1.pkl`
   are joblib pickles: arbitrary code execution on load, and sklearn `StandardScaler`
   unpickles break across sklearn versions. Ship `.npz` as primary; keep the `.pkl`
   alongside in `compat/`, labeled "for byte-exact bioCAV pipeline compatibility."
2. **A CAV is not scoreable from `concept_v1.npy` alone.** `report_v1.json` carries the
   `intercept`. Every stacked coefficient matrix ships with an aligned intercept vector.
3. **Keep float64 for coefficients** — 6,677 × 1,152 float64 ≈ 62 MB, cheap, and preserves
   exact reproduction of the paper's numbers.
4. **Stack, don't ship 7,000 directories.** Per-concept dirs are ~97 KB each, ~90% of it a
   duplicated 28 KB scaler and a 59 KB report. Stacked matrix + metrics TSV + one tarball of
   raw reports cuts each library from ~650 MB to ~70 MB of primary artifacts.
5. **Two copies of every command.** `as_run/` = the verbatim script actually executed
   (provenance). `run/` = a parameterized version a downloader can execute. Note in the README
   that `00_get_negatives.sh` uses short flags (`-ss mean -s -l -11`) that differ cosmetically
   from the README's long-flag form, and that `train_cav_from_span.sh` has a hardcoded
   `biocav_dir` and a *relative* `ref_neg=reference_population/neg_10000.pkl` that only
   resolves from one working directory.
6. **Commands are meaningless without their inputs.** The `.span` files are tiny — ship them
   inside the GO and EC records.
7. Every record: `MANIFEST.tsv` (path, bytes, sha256), `LICENSE`, `CITATION.cff`, and a
   `PROVENANCE.md` naming the exact `mcwlab_utils` tag, ESM++ model dir + revision, conda env
   export, and the source paths on the cluster.

---

## Dataset 1 — Background reference population

**Title:** bioCAV reference background: 10,000 proteins and ESM++ large embeddings

```
background/
├── neg_10000.fasta                     4.5 MB, 10,000 seqs
├── neg_10000_embeddings.npz            (10000, 1152) float32 + seq_ids  ~46 MB
├── scaler_v1.npz                       mean_, scale_, n_samples_seen_, feature dim
├── compat/
│   ├── neg_10000.pkl                   46 MB  (joblib, exact pipeline input)
│   ├── neg_10000.pkl.seqnames
│   ├── neg_10000.pkl.description
│   └── scaler_v1.pkl                   28 KB
├── as_run/00_get_negatives.sh          verbatim, incl. the sbatch header
├── run/embed_background.sh             parameterized
├── run/fit_scaler.sh                   fit_reference_scaler.py --neg ... --out ...
├── PROVENANCE.md                       ESM++ large, layer -11 (= layer 25 of 36), mean pooling,
│                                       half precision, max_length 2048, batch 1
└── MANIFEST.tsv, LICENSE, CITATION.cff, README.md
```

Not shipped: `global_pca_v1*.pkl` — superseded by the scaler, retained only for pre-migration
analyses. Mention in README, don't distribute.

**Embedding command (as run):**
```bash
python hf_embed_new.py -f neg_10000.fasta -o reference_population/neg_10000.pkl \
    -ss mean -s -l -11 -m $ESM_MODEL_DIR -b 1 --max_length 2048
```
**Scaler command:**
```bash
python bioCAV/scripts/fit_reference_scaler.py \
    --neg reference_population/neg_10000.pkl --out reference_population/scaler_v1.pkl
```

Everything in datasets 2–4 is defined relative to this record. Zenodo
`isReferencedBy` / `isSupplementTo` links both ways.

---

## Dataset 2 — GO CAVs

**Title:** bioCAV Gene Ontology concept activation vectors (6,677 GO terms, ESM++ large L25)

```
go_cavs/
├── go_cav_coefficients.npz     coef (6677, 1152) float64 · intercept (6677,) ·
│                               go_ids (6677,) · scaler_source (6677,) {"shared"|"own"}  ← see B1
│                               classify all 6,677 by file mtime (own: ≤2026-05-22), validated
│                               against the 80 already labeled by array content
├── go_cav_metrics.tsv          go_id, namespace, n_positive, n_negative, cv_auroc_mean,
│                               cv_auroc_std, train_auroc, train_auprc, intercept, coef_norm,
│                               concept_cav_hash, scaler_hash, built_date, profab_part
├── go_cav_own_scalers.npz      mean_/scale_ for the ~4,000 concepts that fit their own (B1)
├── spans/                      random_positive_train_max1000.span, one per GO term (tiny)
├── fastas/                     random_positive_train_max1000.span.fasta  (optional, ~500 MB)
├── reports_raw.tar.zst         6,677 × report_v1.json incl. the 20 label-shuffle null CAVs
├── as_run/
│   ├── 02_make_cav_commands.sh
│   ├── train_cav_from_span.sh
│   └── cav_GO_COMMANDS.sh      the 6,678 emitted lines, verbatim
├── run/train_go_cav.sh         parameterized, always passes --scaler-pkl
└── PROVENANCE.md, MANIFEST.tsv, LICENSE, CITATION.cff, README.md
```

**Command trail to state in the README** (embedding → pooling → training, per concept):
```bash
# 1. positives FASTA from the proFAB span file
python bioCAV/specific_scripts/retrieve_fasta.py $SPAN $SPAN.fasta
# 2. embed
python hf_embed_new.py -f $SPAN.fasta -o $SPAN.fasta.pkl \
    --get_aa_embeddings --get_sequence_embedding --strat mean \
    -l -11 -m $ESM_MODEL_DIR -b 1 --max_length 2048
# 3. pool by span
python bioCAV/scripts/prepare_embeddings.py --pkl $SPAN.fasta.pkl \
    --info $SPAN.fasta.pkl.seqnames --spans $SPAN --out $SPAN.pos.npy
# 4. train  (C=1.0, cv_folds=5, seed=42)
python bioCAV/scripts/train_cav_from_embeddings.py --pos $SPAN.pos.npy \
    --neg reference_population/neg_10000.pkl --out ${SPAN%.span}_cav \
    --cv-folds 5 --scaler-pkl reference_population/scaler_v1.pkl
```
Attribution: positive/negative splits from proFAB; sequences from UniProt.

---

## Dataset 3 — EC CAVs

**Title:** bioCAV Enzyme Commission concept activation vectors (952 EC numbers, ESM++ large L25)

Identical layout to dataset 2, with `ec_cav_coefficients.npz` / `ec_cav_metrics.tsv`
(`ec_number` in `ecNo_1-3-1-14` form, plus a normalised `1.3.1.14` column) and
`as_run/cav_EC_COMMANDS.sh`. **No `own_scalers.npz` needed — all 952 use the shared scaler**
(verified on a 60-CAV sample). Same command trail; `--cv-folds 5`, `C=1.0`, seed 42.

Separate record rather than merged with GO because it has different provenance, a different
gold standard (EC-Bench / Price-149), and is clean of B1.

---

## Dataset 4 — Pfam motif CAVs

**Title:** bioCAV Pfam motif concept activation vectors (27,472 Pfam families, ESM++ large L25)

Source: `/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all` (per B2).

```
pfam_cavs/
├── pfam_cav_coefficients.npz     coef (27472, 1152) float64 ≈ 253 MB · intercept · pfam_ids
│                                 13 motifs failed training — exclude, list them in README
├── pfam_cav_metrics.tsv          pfam_id, n_positive, n_negative, cv_auroc_mean/std,
│                                 train_auroc, train_auprc, intercept, coef_norm, hashes
├── pfam_cav_scalers.npz          per-motif mean_/scale_ — REQUIRED here: no Pfam CAV uses the
│                                 shared scaler, and scoring loads the per-motif one (B1).
│                                 (27472, 1152) ×2 float32 ≈ 253 MB
├── motif_thresholds.json         27,459 motifs, regenerated by
│                                 distribution/pfam_thresholds/build_motif_thresholds.py (B2)
├── training_summary.json
├── reports_raw.tar.zst           incl. the 20 label-shuffle null CAVs per motif
├── as_run/                       the library build driver (from shamail)
└── PROVENANCE.md, MANIFEST.tsv, LICENSE, CITATION.cff, README.md
```

**Background provenance — verified directly.** Pfam caches its shared negatives at
`embeddings/ESMplusplus_large/_shared_negatives/L25_all.npy`, so this was checked by comparing
vectors rather than inferred: same 10,000 UniRef50 accessions **in the same order** as
`neg_10000.pkl`, and **median row cosine 0.99999995**, with only **1 row of 10,000** below
0.9999 (min 0.997583). Residual differences are float32-vs-float64 storage and half-precision
nondeterminism; no sequence exceeds the truncation limit, so truncation is not a factor.

This holds even though the two were produced by **different code**: dataset 1 via
`hf_embed_new.py --strategy mean`, Pfam via `esm2_tcav/src/embed.py` with
`shared_neg_pooling: "full_sequence"`. The two agree to numerical precision, so
`isDerivedFrom` → dataset 1 is justified and dataset 1 need not ship a second pooling — but
dataset 4's README must state its own embedding code path, not `hf_embed_new.py` (B3).

Differences from GO/EC to state explicitly in the README:
- `class_weight = "balanced"` (GO/EC use the default), `n_random_cavs = 20`,
  `random_mode = "label_shuffle"`, `use_pca = false`, `C = 1.0`, `cv_folds = 5`, seed 42.
- Positives capped at 1,000 per family.
- Files are `L25_`-prefixed; layer 25 of ESM++ large == the `-l -11` used everywhere else.
- Each motif is scored in **its own** scaler frame, not the shared one — by design here, and
  the reason `pfam_cav_scalers.npz` is mandatory rather than optional.


---

## Suggested rollout

1. **B2-B5 are closed, and licensing is fully settled: CC-BY-4.0 on all four records.**
   The only open blocker is B1 — the GO retrain (job 23885752, queued), which changes
   dataset 2's contents. Do not write the build scripts against dataset 2 until it lands.
2. Write `distribution/build_dataset_{1,2,3,4}.py` — each reads the cluster source paths and
   emits a record directory + MANIFEST; they must be re-runnable, since B1 may change dataset 2.
3. Verification gate: for 50 random concepts per library, score a held-out FASTA from the
   packaged `.npz` and from the original per-concept directory; require bitwise-identical
   rankings before upload.
4. Reserve four Zenodo DOIs, cross-link with `isDerivedFrom` → dataset 1, cite all four in the
   paper's Data Availability.
