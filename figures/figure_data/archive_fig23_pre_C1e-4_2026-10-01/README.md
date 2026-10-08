# Fig 2 / Fig 3 data before the C=1e-4 library swap (archived 2026-10-01)

These 29 files are the figure_data/ versions that were live until 2026-10-01. They were
built from the shipped (early-stopped lbfgs) GO/EC CAVs. Their replacements come from the
C=1e-4 converged shared-scaler library and were built by figures/rerun_fig23_C1e-4.sh
(staged in figure_data_staging/C1e-4/, see its README). The rendered figures from before
the swap are in figures/previews/archive_rendered_pre_C1e-4_2026-10-01/.
To restore: copy these files back into figure_data/ and rerun figures.R.

Not rebuilt (figures.R does not read them; they still hold old-library values):
ec_specificity_ranks.csv, ec_specificity_ranks_samelevel.csv.
