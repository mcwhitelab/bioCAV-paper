# CAV × BLAST integrated score + CAV-rank-1 / BLAST-low cases (2026-10-06, PREVIEW)

Not in figures.R or figure_data/. Pairs, candidates and filters are those of the live Fig 2
rank panel (Fig 2B since the 2026-10-06 relettering; was 2C): C = 1e-4 GO library,
known-positive-excluded pairs, non-obsolete same-ontology terms, strict _pf filter.
`01_integrated_ranks.py` rebuilds the full per-candidate CAV and BLAST score vectors and asserts
that CAV rank, BLAST mid rank and the BLAST hit flag reproduce figure_data/go_rank_families_{ont}.csv
exactly (1.0000 for all three ontologies) before computing anything new.

## Integrated scores (both products; computed per protein over its candidate terms)
- **rank product**: pct(CAV) × pct(BLAST), tie-averaged percentiles (1 = best). Parameter-free.
- **soft**: exp(z_CAV) × (ε + bits / max bits), z over the protein's candidates; ε = 0.1 fixed a
  priori (best BLAST term gains ≤ ln 11 = 2.4 CAV SDs). ε = 0.01 and 1 are reported as a
  sensitivity check only.
- A protein with no BLAST hit to any candidate gets the CAV ranking unchanged.
- Integrated "No prediction" = CAV LLR ≤ 0 AND no BLAST hit to the true term. The gain also holds
  with the CAV-only gate (integrated_summary.csv numbers drop ≤ 2.8 points top-10).
- Caveat: no parameter was fit to these pairs, but the variants were chosen after seeing
  that BLAST and CAV succeed on partly different pairs (README one level up), so this is an
  in-sample comparison, not a held-out one.

## Files
- `integrated_ranks_{ont}.csv` per-pair ranks (CAV, BLAST, rank product, soft ε 0.1/0.01/1, best-of)
- `integrated_summary.csv`, `integrated_summary_by_family.csv`, `integrated_tier_composition.csv`
- `preview_integrated_ranks.{png,pdf}` (`03_preview.R`), Fig 2 rank-panel style
- `cav_rank1_blast_low.csv` (`04_discordant.py`): CAV rank 1 with LLR > 0 and BLAST mid rank
  > 10 or no DIAMOND hit (e ≤ 1e-3) to the true term's training positives. UniProt metadata from the
  REST API; jackhmmer/Pfam status joined from ../callouts_C1e-4/strict_callouts.csv where present
  (the strict set additionally requires DeepGo-SE silence, so most rows here have no jackhmmer info).
  "No hit" is DIAMOND-only and not "no homology": short peptides (e.g. Bowman-Birk inhibitors)
  and paralog clusters appear here.
