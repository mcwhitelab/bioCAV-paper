#!/usr/bin/env Rscript
# DRAFT / preview -- 3 "big receptor" candidates (MET, LDLR, Plexin A1) for
# a higher-domain-count second Figure 4 worked example, scored against
# layer-25 ("display layer 26") CAVs from the 20k-motif library.
# Data from score_bigreceptors.py. Run standalone:
#   Rscript draft_fig4b_bigreceptors.R

source("draft_fig_extra_proteins_common.R")

DATA <- "figure_data/vav_motif_repro/extra_proteins"
OUT  <- "figures"

base_theme <- function(...) theme_cowplot(font_size = 8, ...)

scores <- fromJSON(file.path(DATA, "candidates_bigreceptors_scores.json"), simplifyVector = FALSE)

ACCESSIONS <- c("P08581", "P01130", "Q9UIW2")
TITLES <- c(P08581 = "P08581 - MET (4 domain types, 6 instances)",
           P01130 = "P01130 - LDLR (4 domain types, 15 instances)",
           Q9UIW2 = "Q9UIW2 - Plexin A1 (8 domain types, 13 instances)")

panels <- map(ACCESSIONS, function(acc) {
  build_extra_protein_panel(acc, scores[[acc]], title = TITLES[acc], base_theme_fn = base_theme)
})

fig <- plot_grid(plotlist = panels, ncol = 1, labels = "AUTO", label_size = 10,
                 align = "hv", axis = "tblr")

ggsave(file.path(OUT, "draft_fig4b_bigreceptors.pdf"), fig,
       width = 9, height = 11, bg = "white", device = cairo_pdf)
ggsave(file.path(OUT, "draft_fig4b_bigreceptors.png"), fig,
       width = 9, height = 11, bg = "white", dpi = 200)
message("Saved draft_fig4b_bigreceptors.pdf/.png")
