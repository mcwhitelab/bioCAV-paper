#!/usr/bin/env Rscript
# DRAFT / exploratory -- 3 candidate "second worked example" multidomain
# proteins for Figure 4 (EGFR, FAK1, PIK3CA), each scored against layer-25
# ("display layer 26") CAVs from the 20k-motif library. Shown side by side
# so a pick can be made before promoting one into the real Figure 4.
# Data from score_candidates_fig4b.py. Run standalone:
#   Rscript draft_fig4b_candidates.R

source("draft_fig_extra_proteins_common.R")

DATA <- "figure_data/vav_motif_repro/extra_proteins"
OUT  <- "figures"

base_theme <- function(...) theme_cowplot(font_size = 8, ...)

scores <- fromJSON(file.path(DATA, "candidates_fig4b_scores.json"), simplifyVector = FALSE)

ACCESSIONS <- c("P00533", "Q05397", "P42336")
TITLES <- c(P00533 = "P00533 - EGFR (5 domains)",
           Q05397 = "Q05397 - FAK1 (5 domains)",
           P42336 = "P42336 - PIK3CA (5 domains)")

panels <- map(ACCESSIONS, function(acc) {
  build_extra_protein_panel(acc, scores[[acc]], title = TITLES[acc], base_theme_fn = base_theme)
})

fig <- plot_grid(plotlist = panels, ncol = 1, labels = "AUTO", label_size = 10,
                 align = "hv", axis = "tblr")

ggsave(file.path(OUT, "draft_fig4b_candidates.pdf"), fig,
       width = 8, height = 10.5, bg = "white", device = cairo_pdf)
ggsave(file.path(OUT, "draft_fig4b_candidates.png"), fig,
       width = 8, height = 10.5, bg = "white", dpi = 200)
message("Saved draft_fig4b_candidates.pdf/.png")
