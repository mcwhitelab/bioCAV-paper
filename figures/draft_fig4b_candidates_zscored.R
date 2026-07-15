#!/usr/bin/env Rscript
# DRAFT / preview only -- same 3 candidates as draft_fig4b_candidates.R, but
# with each domain's curve/peaks z-scored against ITS OWN negative-set
# distribution (neg_score_mean/neg_score_std from that domain's
# L25_report_v1.json concept_metrics), instead of raw CAV score. This is
# the same normalization convention already used by Figure 2's intro panel
# (z-score against each term's own negative/background distribution).
# Not saved anywhere permanent yet -- just a comparison preview.

source("draft_fig_extra_proteins_common.R")

DATA <- "figure_data/vav_motif_repro/extra_proteins"
CAV_DIR <- "/xdisk/clairemcwhite/shamail/tcav_outputs_esmplusplus_all/cavs"
OUT  <- "figures"

base_theme <- function(...) theme_cowplot(font_size = 8, ...)

scores <- fromJSON(file.path(DATA, "candidates_fig4b_scores.json"), simplifyVector = FALSE)

get_neg_stats <- function(motif_id) {
  report <- fromJSON(file.path(CAV_DIR, motif_id, "L25_report_v1.json"), simplifyVector = FALSE)
  cm <- report$concept_metrics
  list(mean = cm$neg_score_mean, sd = cm$neg_score_std)
}

zscore_protein_data <- function(protein_data) {
  for (motif_id in names(protein_data$domains)) {
    stats <- get_neg_stats(motif_id)
    dom <- protein_data$domains[[motif_id]]
    dom$curve <- lapply(dom$curve, function(v) if (is.null(v)) NULL else (v - stats$mean) / stats$sd)
    dom$peaks <- lapply(dom$peaks, function(p) { p$score <- (p$score - stats$mean) / stats$sd; p })
    protein_data$domains[[motif_id]] <- dom
  }
  protein_data
}

ACCESSIONS <- c("P00533", "Q05397", "P42336")
TITLES <- c(P00533 = "P00533 - EGFR (5 domains, z-scored)",
           Q05397 = "Q05397 - FAK1 (5 domains, z-scored)",
           P42336 = "P42336 - PIK3CA (5 domains, z-scored)")

scores_z <- map(scores[ACCESSIONS], zscore_protein_data)

panels <- map(ACCESSIONS, function(acc) {
  build_extra_protein_panel(acc, scores_z[[acc]], title = TITLES[acc], base_theme_fn = base_theme)
})

fig <- plot_grid(plotlist = panels, ncol = 1, labels = "AUTO", label_size = 10,
                 align = "hv", axis = "tblr")

ggsave(file.path(OUT, "draft_fig4b_candidates_zscored.pdf"), fig,
       width = 8, height = 10.5, bg = "white", device = cairo_pdf)
ggsave(file.path(OUT, "draft_fig4b_candidates_zscored.png"), fig,
       width = 8, height = 10.5, bg = "white", dpi = 200)
message("Saved draft_fig4b_candidates_zscored.pdf/.png")
