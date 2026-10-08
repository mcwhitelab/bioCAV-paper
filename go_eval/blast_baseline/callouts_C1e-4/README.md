# GO callout examples on the C = 1e-4 library (2026-10-01)

Replaces `../callouts/` (MF) and `../bp_cav_only_rank1_cases.csv` (BP), which were picked
on the old CAV library and are kept unchanged for provenance.

## Selection rules
Both rules were recovered from the old files and reproduce them exactly on the data they
were picked from.
- **strict** (old `mf_callouts.csv`): CAV rank 1 among the Fig 2C candidates
  (`cav_rank_ont_filt_pf`), LLR > 0, DeepGo-SE emitted nothing for the term, and no DIAMOND
  hit to the term's training positives.
- **cav_only** (old `bp_cav_only_rank1_cases.csv`): CAV rank 1 with LLR > 0. DeepGo-SE is
  not in its top 10 (or not predicted), and BLAST's mid rank is not in its top 10 (or there
  is no hit).

  The old BP list (62) was picked on data that still included obsolete GO terms. On the
  obsolete-excluded data the old library gives 56.
- **Input:** live `figures/figure_data/go_rank_families_{ont}.csv` (C = 1e-4).

## Homology annotation (strict set)
The recipe was recovered from the transcript of the session that built the old callouts.
- **jackhmmer:** pyhmmer 0.12.3, 3 iterations, E = 10, incE = 1e-3. Each candidate is
  searched against `train_positives_all.fa`, which already contains every GO:0032436
  positive. Hits with E <= 1e-3 are ranked with `blast_ranks.py` over the Fig 2C
  candidates.
- **Pfam:** UniProt Pfam cross-references, shared with the true term's training positives.
  Computed only for candidates without a jackhmmer hit.
- **Tier:** jackhmmer hit > shared Pfam > none. `callout = True` means no detectable
  homology to the true term's positives.

**Calibration** (`calibration.txt`): the 60 old strict queries were re-run in the same job.
- Their jackhmmer hits are identical to the old run (12,599 hits).
- The jackhmmer hit flag and best E-value reproduce for all 11 old MF callouts.

## Files
- `strict_callouts.csv`: all ontologies, with an `ontology` column.
- `mf_callouts.csv`, `mf_callouts_cav_runnerups.csv`: same format as the old files.
- `{mf,bp,cc}_cav_only_rank1_cases.csv`: with CAV runner-ups, DeepGo-SE and BLAST top
  terms, the nearest DIAMOND neighbour and known GPAD terms, all from the C = 1e-4 library
  and extended BLAST.
- **Scripts:** `01_select.py` → `02_run_jack.sbatch` (`02_run_jack.py`) → `03_build_tables.py`.

## Counts
| | strict | strict, no detectable homology | cav_only | Old library (strict / no-homology / cav_only) |
|---|---|---|---|---|
| MF | 13 | 6 | 17 | 11 / 4 / 15 |
| BP | 28 | 11 | 58 | — / — / 56 (62 as picked) |
| CC | 50 | 34 | 72 | — |

Of the old MF "no detectable homology" callouts, only Q2RVS1 (EncFtn, ferroxidase activity)
survives. It is still rank 1 with no jackhmmer hit or shared Pfam.
