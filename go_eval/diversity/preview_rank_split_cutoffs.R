# Exploratory preview (2026-09-29): Fig 2C in the figures.R style, split by
# whether the GO term's CAV training positives span <K or >=K sequence
# families, for K in CUTOFFS (term_diversity.py). Known-positive-excluded data; gated ranks.
# Run from bioCAV-paper/go_eval/diversity/.
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

STAGE <- "../../figures/figure_data_staging/noknown"
FIG_FONT <- "Arial"
RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
RANK_COLORS <- c("Rank 1" = "#005A32", "Rank 2-3" = "#238B45", "Rank 4-10" = "#74C476",
                 "Rank >10" = "#C4B7DC", "No prediction" = "#8073AC")
ONTS <- c("mf", "bp", "cc")
ONT_DISPLAY <- c("MF", "BP", "CC")

tier_of <- function(rank, predicted) {
  case_when(!predicted ~ "No prediction", rank == 1 ~ "Rank 1", rank <= 3 ~ "Rank 2-3",
            rank <= 10 ~ "Rank 4-10", TRUE ~ "Rank >10")
}

CUTOFFS <- c(10, 15, 20, 25)

pairs <- map_dfr(ONTS, \(o) {
  read_csv(sprintf("%s/go_specificity_ranks_allterms_%s.csv", STAGE, o), show_col_types = FALSE) |>
    inner_join(read_csv(sprintf("term_diversity_%s.csv", o), show_col_types = FALSE),
               by = "go_term") |>
    mutate(ontology = toupper(o))
})

make_split <- function(k) {
  pairs <- pairs |>
    mutate(split = factor(if_else(n_families < k, sprintf("<%d families", k),
                                  sprintf(">=%d families", k)),
                          levels = c(sprintf("<%d families", k), sprintf(">=%d families", k))))
  counts <- pairs |> group_by(split, ontology) |>
    summarise(terms = n_distinct(go_term), pairs = n(), .groups = "drop")
  split_lab <- pairs |> group_by(split) |>
    summarise(lab = sprintf("%s\n(%d terms, %d pairs)", first(split), n_distinct(go_term), n()),
              .groups = "drop") |> deframe()

  rank_df <- pairs |>
    transmute(split, ontology,
              CAV      = tier_of(cav_rank_ont_filt, llr > 0),
              DeepGoSE = tier_of(tool_rank_ont_filt, tool_predicted == TRUE)) |>
    pivot_longer(c(CAV, DeepGoSE), names_to = "method", values_to = "tier") |>
    count(split, ontology, method, tier) |>
    group_by(split, ontology, method) |> mutate(prop = n / sum(n)) |> ungroup() |>
    left_join(counts, by = c("split", "ontology")) |>
    mutate(tier = factor(tier, levels = RANK_LEVELS),
           method = factor(method, levels = c("CAV", "DeepGoSE")),
           ont_lab = factor(sprintf("%s (n=%d)", ontology, pairs)),
           ontology = factor(ontology, levels = rev(ONT_DISPLAY)))

  fig <- ggplot(rank_df, aes(x = prop, y = ontology, fill = tier)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    facet_grid(method ~ split, switch = "y", labeller = labeller(split = split_lab)) +
    scale_fill_manual(values = RANK_COLORS, breaks = RANK_LEVELS, name = NULL) +
    scale_x_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = c(0, 0.02))) +
    theme_cowplot(font_size = 8, font_family = FIG_FONT) +
    labs(x = "Share of validation protein-GO pairs", y = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8),
          legend.position = "top", legend.justification = "left",
          legend.key.size = unit(6, "pt"), strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 8, margin = margin(l = 4, r = 6)),
          strip.text.x = element_text(size = 8), panel.spacing.x = unit(22, "pt"),
          panel.spacing.y = unit(10, "pt"))
  summ <- rank_df |> group_by(split, ontology, method) |>
    summarise(rank1 = round(100 * sum(prop[tier == "Rank 1"]), 1),
              top10 = round(100 * sum(prop[tier %in% RANK_LEVELS[1:3]]), 1), .groups = "drop") |>
    mutate(cutoff = k)
  list(fig = fig, summ = summ, counts = counts |> mutate(cutoff = k))
}

res <- map(CUTOFFS, make_split)
walk2(res, CUTOFFS, \(r, k) ggsave(sprintf("preview_rank_split%d.png", k), r$fig,
                                    width = 7.2, height = 2.9, dpi = 300, bg = "white"))
panels <- map2(res, seq_along(res), \(r, i) if (i == 1) r$fig else r$fig + theme(legend.position = "none"))
combined <- plot_grid(plotlist = panels, ncol = 1, labels = sprintf("cutoff %d", CUTOFFS),
                      label_size = 8, label_fontfamily = FIG_FONT, hjust = 0, label_x = 0.01,
                      rel_heights = c(1.1, 1, 1, 1))
ggsave("preview_rank_split_cutoffs.png", combined, width = 7.2, height = 10.5, dpi = 300, bg = "white")
ggsave("preview_rank_split_cutoffs.pdf", combined, width = 7.2, height = 10.5, device = cairo_pdf)

print(bind_rows(map(res, "counts")) |> arrange(cutoff, split, ontology), n = 30)
bind_rows(map(res, "summ")) |>
  pivot_wider(names_from = method, values_from = c(rank1, top10)) |>
  arrange(ontology, cutoff, split) |> print(n = 30)
