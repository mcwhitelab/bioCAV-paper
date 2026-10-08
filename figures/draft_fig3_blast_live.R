# Preview (2026-09-29), not part of figures.R. Fig 3 changes proposed for review:
#   previews/fig3A_v1_ecbench_blastp.png  panel A + EC-Bench's published blastp
#       (Correct/Incorrect, scored like CLEAN-SupCon)
#   previews/fig3A_v2_diamond_rank.png    panel A + DIAMOND nearest-neighbour rank
#       over the same candidates as CAV (rank tiers, like Fig 2C's BLAST)
#   previews/fig3E_umap_live.png          panel E without GO terms obsolete in go.obo
#       and EC numbers deleted/transferred in enzyme.dat (same UMAP coordinates)
# Both A versions use cav_rank_level_live (deleted/transferred ECs excluded; changes 1 of
# 209 ranks). Source: ec_eval/blast_baseline/ec_blast_ranks.py -> ec_blast_ranks.csv.
# Run from bioCAV-paper/figures/.
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggrastr)
})

FIG_FONT <- "Arial"
RASTER_DPI <- 150
OUT <- "previews"
oi <- c(orange = "#E69F00", sky_blue = "#56B4E9", green = "#009E73", yellow = "#F0E442",
        blue = "#0072B2", vermillion = "#D55E00", pink = "#CC79A7", black = "#000000")
base_theme <- function(...) {
  theme_cowplot(font_size = 8, font_family = FIG_FONT, ...) +
    theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8),
          strip.text = element_text(size = 8))
}
RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
RANK_COLORS <- c("Rank 1" = "#005A32", "Rank 2-3" = "#238B45", "Rank 4-10" = "#74C476",
                 "Rank >10" = "#C4B7DC", "No prediction" = "#8073AC")
tier_of <- function(rank, predicted) {
  case_when(!predicted ~ "No prediction", rank < 1.5 ~ "Rank 1", rank <= 3 ~ "Rank 2-3",
            rank <= 10 ~ "Rank 4-10", TRUE ~ "Rank >10")
}
EC_TIER_LEVELS <- c(RANK_LEVELS, "Correct", "Incorrect")
EC_TIER_COLORS <- c(RANK_COLORS, "Correct" = unname(oi["pink"]), "Incorrect" = "#4D4D4D")
EC_LEVEL_LABELS <- c("3" = "EC level 3", "4" = "EC level 4\n(fully specific)")

ec <- read_csv("../ec_eval/blast_baseline/ec_blast_ranks.csv", show_col_types = FALSE)

panel_a <- function(blast_rows, method_levels) {
  df <- bind_rows(
    ec |> transmute(level, method = "CAV", tier = tier_of(cav_rank_level_live, llr > 0)),
    blast_rows,
    ec |> transmute(level, method = "CLEAN-SupCon",
                    tier = if_else(clean_supcon_correct, "Correct", "Incorrect"))
  ) |>
    count(level, method, tier) |>
    group_by(level, method) |> mutate(prop = n / sum(n)) |> ungroup() |>
    mutate(tier = factor(tier, levels = EC_TIER_LEVELS),
           method = factor(method, levels = rev(method_levels)),
           level = factor(EC_LEVEL_LABELS[as.character(level)], levels = EC_LEVEL_LABELS[c("4", "3")]))
  ggplot(df, aes(x = prop, y = method, fill = tier)) +
    geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
    facet_grid(level ~ ., switch = "y") +
    scale_fill_manual(values = EC_TIER_COLORS, breaks = EC_TIER_LEVELS, name = NULL) +
    scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    base_theme() +
    labs(x = paste0("Share of protein-EC pairs (n=", nrow(ec), ")"), y = NULL) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "top", legend.location = "plot", legend.justification = "left",
          legend.key.size = unit(6, "pt"), legend.margin = margin(t = 0, b = 2, r = 0, l = 16),
          legend.spacing.x = unit(0, "pt"), strip.background = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 8),
          panel.spacing.y = unit(4, "pt"), plot.margin = margin(t = 2, r = 6, b = 2, l = 2))
}

# Width/height match panel A's slot in cav_fig3 (1.55/2.55 of 7.2in, 2.2in row).
save_a <- function(p, name) {
  ggsave(file.path(OUT, name), p, width = 4.4, height = 2.4, dpi = 300, bg = "white")
  ggsave(file.path(OUT, sub("png$", "pdf", name)), p, width = 4.4, height = 2.4, device = cairo_pdf)
}

v1 <- panel_a(ec |> transmute(level, method = "BLASTp",
                              tier = if_else(ecbench_blastp_correct, "Correct", "Incorrect")),
              c("CAV", "BLASTp", "CLEAN-SupCon"))
save_a(v1, "fig3A_v1_ecbench_blastp.png")

v2 <- panel_a(ec |> transmute(level, method = "BLAST", tier = tier_of(blast_rank_mid, blast_predicted)),
              c("CAV", "BLAST", "CLEAN-SupCon"))
save_a(v2, "fig3A_v2_diamond_rank.png")

# ---- Panel E without obsolete GO terms / deleted-or-transferred ECs -----------------
obo <- readLines("../go_eval/static/go_annotations/go.obo")
starts <- grep("^\\[", obo)
ids <- sub("^id: ", "", obo[grepl("^id: GO:", obo)])
id_line <- grep("^id: GO:", obo)
blk <- findInterval(seq_along(obo), starts)
id_of_blk <- setNames(ids, blk[id_line])
obsolete_go <- unname(id_of_blk[as.character(blk[grep("^is_obsolete: true", obo)])])

enz <- readLines("../ec_eval/static/ec_annotations/enzyme.dat")
enz_id <- cumsum(grepl("^ID   ", enz))
id_by_blk <- setNames(sub("^ID   ", "", enz[grepl("^ID   ", enz)]), unique(enz_id[grepl("^ID   ", enz)]))
dead_ec <- unname(id_by_blk[as.character(unique(enz_id[grepl("^DE   (Transferred|Deleted) entry", enz)]))])

umap <- read_csv("figure_data/go_ec_combined_umap.csv", show_col_types = FALSE)
live <- umap |> filter(!(category == "EC" & concept %in% dead_ec),
                       !(category != "EC" & concept %in% obsolete_go))
message(sprintf("UMAP: %d -> %d concepts (dropped %d obsolete GO, %d deleted/transferred EC)",
                nrow(umap), nrow(live),
                sum(umap$category != "EC" & umap$concept %in% obsolete_go),
                sum(umap$category == "EC" & umap$concept %in% dead_ec)))

COMBINED_CAT_COLORS <- c(EC = unname(oi["vermillion"]), MF = unname(oi["blue"]),
                         BP = unname(oi["green"]), CC = unname(oi["pink"]))
cat_counts <- live |> count(category)
cat_labels <- setNames(paste0(cat_counts$category, " (n=", cat_counts$n, ")"), cat_counts$category)
p_e <- live |>
  mutate(category = factor(category, levels = c("EC", "MF", "BP", "CC"))) |>
  ggplot(aes(D1, D2, color = category)) +
  geom_point_rast(size = 0.35, alpha = 0.55, raster.dpi = RASTER_DPI) +
  scale_color_manual(values = COMBINED_CAT_COLORS, labels = cat_labels, name = NULL) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.5, alpha = 1))) +
  base_theme() +
  labs(x = "UMAP 1", y = "UMAP 2", title = paste0("Combined CAV space (n=", nrow(live), ")")) +
  theme(plot.title = element_text(size = 8, face = "bold"), legend.position = "bottom",
        legend.key.size = unit(8, "pt"), axis.text = element_blank(), axis.ticks = element_blank())
ggsave(file.path(OUT, "fig3E_umap_live.png"), p_e, width = 2.4, height = 1.9 * 1.2, dpi = 300, bg = "white")
ggsave(file.path(OUT, "fig3E_umap_live.pdf"), p_e, width = 2.4, height = 1.9 * 1.2, device = cairo_pdf)
