#!/usr/bin/env Rscript
# DRAFT / exploratory -- reproduction of arxiv 2511.21614v1 Appendix A
# (7-protein grid of layer-26 CAV score profiles: P29317, P03001, Q4R4U2,
# P48960, Q9ZQ85, Q8CGX0, Q4JG17), plus FAK1 (moved here from Figure 4 --
# its "clean multidomain" point is already made by Plexin A1 there) and
# SET1_SCHPO/PLDZ_DICDI (also moved from Figure 4, to keep that figure to
# one row -- these two still make the "unannotated second occurrence"
# point on their own here). CAVs are layer-25 ("display layer 26") from
# the 20k-motif library at
# /xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs (per user
# instruction to use the new 20k CAVs). Data from score_extra_proteins.py
# (also covers SET1_SCHPO/PLDZ_DICDI) and score_candidates_fig4b.py (FAK1).
# Slated to become a supplemental figure once validated. Not yet part of
# the figures.R pipeline; run standalone with:
#   Rscript draft_fig_supp_appendixA.R

source("draft_fig_extra_proteins_common.R")

DATA <- "figure_data/vav_motif_repro/extra_proteins"
OUT  <- "figures"

base_theme <- function(...) theme_cowplot(font_size = 8, ...)

scores <- fromJSON(file.path(DATA, "extra_proteins_scores.json"), simplifyVector = FALSE)
cand_scores <- fromJSON(file.path(DATA, "candidates_fig4b_scores.json"), simplifyVector = FALSE)
scores[["Q05397"]] <- cand_scores[["Q05397"]]

ACCESSIONS <- c("P29317", "P03001", "Q4R4U2", "P48960", "Q9ZQ85", "Q8CGX0", "Q4JG17", "Q05397",
                "Q9Y7R4", "Q54SA1")
TITLES <- setNames(ACCESSIONS, ACCESSIONS)
TITLES["Q05397"] <- "Q05397 - FAK1"
TITLES["Q9Y7R4"] <- "Q9Y7R4 - SET1_SCHPO"
TITLES["Q54SA1"] <- "Q54SA1 - PLDZ_DICDI"

panels <- map(ACCESSIONS, function(acc) {
  build_extra_protein_panel(acc, scores[[acc]], title = TITLES[acc], base_theme_fn = base_theme)
})

fig <- plot_grid(plotlist = panels, ncol = 2, labels = "AUTO", label_size = 9,
                 align = "hv", axis = "tblr")

ggsave(file.path(OUT, "draft_fig_supp_appendixA.pdf"), fig,
       width = 9, height = 16, bg = "white", device = cairo_pdf)
ggsave(file.path(OUT, "draft_fig_supp_appendixA.png"), fig,
       width = 9, height = 16, bg = "white", dpi = 200)
message("Saved draft_fig_supp_appendixA.pdf/.png")
