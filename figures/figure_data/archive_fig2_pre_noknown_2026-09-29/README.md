# Fig 2 data before the known-positive exclusion (archived 2026-09-29)

These are the Fig 2 CSVs as shipped (July eval, leaky validation pairs, obsolete GO terms
among the rank candidates). They were replaced by:
- temporal_* and go_specificity_ranks_{ont}.csv from figures/figure_data_staging/noknown/
  (go_eval/06_eval_exclude_known_positives.sh)
- go_specificity_ranks_allterms_{ont}.csv from go_eval/blast_baseline/noobsolete/
  (strict _pf filter, obsolete terms excluded)
Restore: copy these files back into figures/figure_data/ and delete the mAP=0.674/0.729/0.731
density and pr_curves files.
