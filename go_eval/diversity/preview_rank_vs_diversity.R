# Exploratory preview (2026-09-29): Fig 2C-style ranks vs. sequence diversity of
# each GO term's CAV training positives. Known-positive-excluded data only.
#   A  per-term share of validation pairs ranked in the top 10 (gated, as in
#      panel C) vs. families per 100 training positives
#   B  same for Rank 1
#   C  rank composition of pairs, binned by their term's training diversity
#   D  per-pair CAV rank as a percentile of the ontology's candidate terms
#      (1 = top), ungated, vs. diversity
# Run from bioCAV-paper/go_eval/diversity/.
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

STAGE <- "../../figures/figure_data_staging/noknown"
FIG_FONT <- "Arial"
METHOD_COLORS <- c(CAV = "#0072B2", DeepGoSE = "#D55E00")
RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
RANK_COLORS <- c("Rank 1" = "#005A32", "Rank 2-3" = "#238B45", "Rank 4-10" = "#74C476",
                 "Rank >10" = "#C4B7DC", "No prediction" = "#8073AC")
DIV_LEVELS <- c("<15", "15-30", "30-45", ">=45")
base_theme <- function() {
  theme_cowplot(font_size = 8, font_family = FIG_FONT) +
    theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8),
          strip.text = element_text(size = 8), strip.background = element_blank(),
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3))
}
ONTS <- c("mf", "bp", "cc")
tier <- function(rank, nopred) {
  case_when(nopred ~ "No prediction", rank == 1 ~ "Rank 1", rank <= 3 ~ "Rank 2-3",
            rank <= 10 ~ "Rank 4-10", TRUE ~ "Rank >10")
}

pairs <- map_dfr(ONTS, \(o) {
  read_csv(sprintf("%s/go_specificity_ranks_allterms_%s.csv", STAGE, o), show_col_types = FALSE) |>
    inner_join(read_csv(sprintf("term_diversity_%s.csv", o), show_col_types = FALSE),
               by = "go_term") |>
    mutate(ontology = toupper(o))
}) |>
  mutate(ontology = factor(ontology, levels = toupper(ONTS)),
         div_bin = cut(families_per_100, c(-Inf, 15, 30, 45, Inf), right = FALSE,
                       labels = DIV_LEVELS),
         cav_tier = tier(cav_rank_ont_filt, llr <= 0),
         ds_tier = tier(tool_rank_ont_filt, !as.logical(tool_predicted)),
         cav_pct = 1 - (cav_rank_ont_filt - 1) / n_terms_ont)

long <- bind_rows(
  pairs |> transmute(ontology, go_term, families_per_100, div_bin, method = "CAV", tier = cav_tier),
  pairs |> transmute(ontology, go_term, families_per_100, div_bin, method = "DeepGoSE", tier = ds_tier)
) |> mutate(method = factor(method, levels = names(METHOD_COLORS)),
            tier = factor(tier, levels = rev(RANK_LEVELS)))

per_term <- long |>
  group_by(ontology, go_term, families_per_100, method) |>
  summarise(n = n(), top10 = mean(tier %in% c("Rank 1", "Rank 2-3", "Rank 4-10")),
            rank1 = mean(tier == "Rank 1"), .groups = "drop")

term_panel <- function(col, ylab) {
  ggplot(per_term, aes(families_per_100, .data[[col]], color = method)) +
    geom_point(aes(size = n), alpha = 0.3, stroke = 0) +
    geom_smooth(aes(weight = n), method = "loess", formula = y ~ x, se = FALSE,
                linewidth = 0.7, span = 0.9) +
    facet_wrap(~ontology, nrow = 1) +
    scale_color_manual(values = METHOD_COLORS, name = NULL) +
    scale_size_area(max_size = 2.5, guide = "none") +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    base_theme() +
    labs(x = "Sequence families per 100 CAV training positives", y = ylab)
}
p_a <- term_panel("top10", "Pairs in top 10\n(per term)") + theme(legend.position = "top")
p_b <- term_panel("rank1", "Pairs at Rank 1\n(per term)") + theme(legend.position = "none")

bin_n <- pairs |> count(ontology, div_bin)
p_c <- ggplot(long, aes(div_bin, fill = tier)) +
  geom_bar(position = "fill", width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(data = bin_n, aes(div_bin, 1.05, label = n), inherit.aes = FALSE,
            size = 2.2, family = FIG_FONT, color = "gray30") +
  facet_grid(method ~ ontology) +
  scale_fill_manual(values = RANK_COLORS, breaks = RANK_LEVELS, name = NULL) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1),
                     limits = c(0, 1.1), expand = c(0, 0)) +
  base_theme() +
  labs(x = "Term's sequence families per 100 CAV training positives", y = "Share of pairs") +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

p_d <- ggplot(pairs, aes(div_bin, cav_pct)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.3,
               color = METHOD_COLORS[["CAV"]], fill = METHOD_COLORS[["CAV"]], alpha = 0.3) +
  facet_wrap(~ontology, nrow = 1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  base_theme() +
  labs(x = "Term's sequence families per 100 CAV training positives",
       y = "CAV rank percentile\n(100% = top)")

fig <- plot_grid(p_a, p_b, p_c, p_d, ncol = 1, labels = LETTERS[1:4], label_size = 8,
                 label_fontfamily = FIG_FONT, rel_heights = c(1.05, 0.9, 1.5, 0.9))
ggsave("preview_rank_vs_diversity.png", fig, width = 7.2, height = 9, dpi = 300, bg = "white")
ggsave("preview_rank_vs_diversity.pdf", fig, width = 7.2, height = 9, device = cairo_pdf)

cat("\nRank composition by term-diversity bin (gated):\n")
long |> group_by(ontology, div_bin, method) |>
  summarise(n = n(), n_terms = n_distinct(go_term), rank1 = mean(tier == "Rank 1"),
            top10 = mean(tier %in% c("Rank 1", "Rank 2-3", "Rank 4-10")),
            nopred = mean(tier == "No prediction"), .groups = "drop") |>
  mutate(across(rank1:nopred, \(x) round(100 * x, 1))) |> print(n = 30)
cat("\nPer-term Spearman rho, families_per_100 vs top-10 share:\n")
per_term |> group_by(ontology, method) |>
  summarise(rho = cor(families_per_100, top10, method = "spearman"), .groups = "drop") |> print()
cat("\nCAV median rank percentile by bin:\n")
pairs |> group_by(ontology, div_bin) |>
  summarise(median_pct = round(100 * median(cav_pct), 1), .groups = "drop") |> print(n = 20)
