# GO homology-transfer baseline (2026-09-29, exploratory)

DIAMOND nearest-neighbour baseline for Fig 2C on the known-positive-excluded pairs
(`figures/figure_data_staging/noknown`). Not in figures.R.

- `train_positives_all.fa` / `.dmnd` — union of the training positives of all 6,677 GO CAVs
  (128,187 proteins); `term_members.tsv` maps term → member. Two sequences had a stray `,`
  (a retrieval artifact), which was stripped.
- `val_queries.fa` — the 6,484 kept validation proteins. `val_vs_train.tsv` —
  `diamond blastp --more-sensitive --max-target-seqs 0 --evalue 1e-3` (2.0.9, 1 thread, ~90 min).
- `blast_ranks.py` → `blast_ranks_{ont}.csv`. Score = best bitscore to a term's training
  positives, self-hits removed. Same candidates as `*_rank_ont_filt_pf`. Ties are frequent,
  because one hit is a positive for many related terms, so both an optimistic and a
  tie-averaged ("mid") rank are kept. "No prediction" = no hit to the true term.

Result (Rank 1 / top 10 %, gated, `_pf` filter; BLAST mid-tie):

| | CAV | DeepGo-SE | BLAST |
|---|---|---|---|
| MF all (1,226) | 15.8 / 42.2 | 8.4 / 22.5 | 8.3 / 33.5 |
| MF, BLAST hit (553) | 32.0 / 76.9 | 17.9 / 43.4 | 18.4 / 74.3 |
| MF, no hit (673) | 2.5 / 13.7 | 0.6 / 5.3 | 0 / 0 |
| BP all (2,742) | 5.4 / 19.8 | 0.4 / 3.9 | 1.9 / 11.5 |
| BP, BLAST hit (1,100) | 11.7 / 38.9 | 0.6 / 5.9 | 4.6 / 28.7 |
| BP, no hit (1,642) | 1.2 / 6.9 | 0.2 / 2.5 | 0 / 0 |
| CC all (4,055) | 7.3 / 27.6 | 2.7 / 15.3 | 7.5 / 29.0 |
| CC, BLAST hit (1,628) | 13.9 / 52.1 | 5.5 / 25.2 | 18.8 / 72.3 |
| CC, no hit (2,427) | 2.8 / 11.2 | 0.9 / 8.6 | 0 / 0 |

Optimistic-tie BLAST Rank 1 is much higher (MF 22.0, BP 11.1, CC 20.8). Top 10 CAV ∪ BLAST
vs BLAST alone: MF 46.5 vs 33.5, BP 23.9 vs 11.5, CC 38.8 vs 29.0.

**2026-09-29, obsolete terms.** The table above includes the 484 obsolete GO terms as
candidates. With them excluded (now the default in `compute_ranks_all_terms.parse_obo`;
outputs in `noobsolete/`), see `METHODS_FIGURE_STEPS.md` §2 C "Obsolete terms" for the
before/after numbers. `preview_success_vs_identity.R` reads `noobsolete/`;
`preview_ranks_pf.R` still reads the older staged CSVs and `blast_ranks_*.csv` here.
