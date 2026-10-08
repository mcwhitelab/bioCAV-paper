# Exploratory previews (2026-09-29): Fig 3A after known-positive exclusion, with the
# strict candidate filter and a DIAMOND nearest-neighbour baseline
# (ec_ranks_noknown.py). Not part of figures.R. Run from ec_eval/blast_baseline/.
#   preview_EC_pf.png          CAV / BLAST (ours, ranked) / CLEAN-SupCon / EC-Bench
#                              blastp (single prediction) x EC level
#   preview_EC_ownblast.png    CAV / BLAST (ours, ranked) / CLEAN-SupCon, by level
#   preview_EC_ecbench.png     CAV / CLEAN-SupCon / EC-Bench blastp only, by level
#   preview_EC_by_identity.png same, split by the protein's best identity to the
#                              true EC's CAV training positives
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

FIG_FONT <- "Arial"
RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
TIER_LEVELS <- c(RANK_LEVELS, "Correct", "Incorrect")
TIER_COLORS <- c("Rank 1" = "#005A32", "Rank 2-3" = "#238B45", "Rank 4-10" = "#74C476",
                 "Rank >10" = "#C4B7DC", "No prediction" = "#8073AC",
                 "Correct" = "#CC79A7", "Incorrect" = "#4D4D4D")
METHODS <- c("CAV", "BLAST (ranked)", "CLEAN-SupCon", "EC-Bench blastp")
LEVEL_LABELS <- c("3" = "EC level 3", "4" = "EC level 4")
ID_LEVELS <- c("No hit", "<30%", "30-50%", ">=50%")

tier_of <- function(rank, predicted) {
  case_when(!predicted ~ "No prediction", rank < 1.5 ~ "Rank 1", rank <= 3 ~ "Rank 2-3",
            rank <= 10 ~ "Rank 4-10", TRUE ~ "Rank >10")
}

ec <- read_csv("../../figures/figure_data_staging/noknown/ec_ranks_noknown.csv",
               show_col_types = FALSE) |>
  filter(!known_positive) |>
  mutate(id_bin = case_when(!blast_predicted ~ "No hit", blast_true_pident < 30 ~ "<30%",
                            blast_true_pident < 50 ~ "30-50%", TRUE ~ ">=50%") |>
           factor(levels = ID_LEVELS),
         all = factor("All leak-free pairs"),
         level = factor(LEVEL_LABELS[as.character(level)], levels = LEVEL_LABELS[c("4", "3")]))

long <- ec |>
  transmute(level, id_bin, all,
            `CAV` = tier_of(cav_rank_level_filt_pf, llr > 0),
            `BLAST (ranked)` = tier_of(blast_rank_mid, blast_predicted),
            `CLEAN-SupCon` = if_else(clean_supcon_correct, "Correct", "Incorrect"),
            `EC-Bench blastp` = if_else(ecbench_blastp_correct, "Correct", "Incorrect")) |>
  pivot_longer(all_of(METHODS), names_to = "method", values_to = "tier")

ec_plot <- function(split_col, fname, width, height, methods = METHODS) {
  df <- long |> filter(method %in% methods) |>
    count(split = .data[[split_col]], level, method, tier) |>
    group_by(split, level, method) |> mutate(prop = n / sum(n)) |> ungroup() |>
    mutate(tier = factor(tier, levels = TIER_LEVELS),
           method = factor(method, levels = rev(methods)))
  n_lab <- ec |> count(split = .data[[split_col]], level) |>
    mutate(method = factor("CAV", levels = rev(methods)))
  p <- ggplot(df, aes(x = prop, y = method, fill = tier)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    geom_text(data = n_lab, aes(x = 1.02, y = method, label = sprintf("n=%d", n)),
              inherit.aes = FALSE, hjust = 0, size = 2, family = FIG_FONT, color = "gray30") +
    facet_grid(level ~ split, switch = "y") +
    scale_fill_manual(values = TIER_COLORS, breaks = TIER_LEVELS, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1),
                       limits = c(0, 1.2), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_cowplot(font_size = 8, font_family = FIG_FONT) +
    labs(x = "Share of leak-free validation protein-EC pairs", y = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8),
          legend.position = "top", legend.justification = "left",
          legend.key.size = unit(6, "pt"), strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 8, margin = margin(l = 4, r = 6)),
          strip.text.x = element_text(size = 8), panel.spacing.x = unit(10, "pt"),
          panel.spacing.y = unit(8, "pt"))
  ggsave(fname, p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(sub("png$", "pdf", fname), p, width = width, height = height, device = cairo_pdf)
}

ec_plot("all", "preview_EC_pf.png", 6.6, 2.6)
ec_plot("id_bin", "preview_EC_by_identity.png", 7.2, 2.8)
ec_plot("all", "preview_EC_ecbench.png", 6.6, 2.2,
        methods = c("CAV", "CLEAN-SupCon", "EC-Bench blastp"))
ec_plot("all", "preview_EC_ownblast.png", 6.6, 2.2,
        methods = c("CAV", "BLAST (ranked)", "CLEAN-SupCon"))
