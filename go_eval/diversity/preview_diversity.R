# Exploratory preview (2026-09-29): Fig 2 performance vs. GO-term sequence
# diversity and per-protein sequence novelty. Known-positive-excluded data only
# (figures/figure_data_staging/noknown). Not part of figures.R.
#   A  per-term AUROC vs. families per 100 CAV training positives
#   B  per-term AUPR, same x
#   C  rank composition of (protein, term) pairs binned by the protein's best
#      identity to that term's CAV training positives
# Run from bioCAV-paper/go_eval/diversity/.
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

STAGE <- "../../figures/figure_data_staging/noknown"
FIG_FONT <- "Arial"
CAV_COLOR <- "#0072B2"; TOOL_COLOR <- "#D55E00"
METHOD_COLORS <- c(CAV = CAV_COLOR, DeepGoSE = TOOL_COLOR)
RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
RANK_COLORS <- c("Rank 1" = "#005A32", "Rank 2-3" = "#238B45", "Rank 4-10" = "#74C476",
                 "Rank >10" = "#C4B7DC", "No prediction" = "#8073AC")
base_theme <- function() {
  theme_cowplot(font_size = 8, font_family = FIG_FONT) +
    theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8),
          strip.text = element_text(size = 8), strip.background = element_blank())
}
ONTS <- c("mf", "bp", "cc")

terms <- map_dfr(ONTS, \(o) {
  read_csv(sprintf("term_diversity_%s.csv", o), show_col_types = FALSE) |>
    left_join(read_csv(sprintf("%s/temporal_tool_comparison_%s.csv", STAGE, o),
                       show_col_types = FALSE), by = "go_term") |>
    mutate(ontology = toupper(o))
}) |> mutate(ontology = factor(ontology, levels = toupper(ONTS)))

term_long <- terms |>
  transmute(ontology, go_term, families_per_100,
            CAV_AUROC = auc_val_vs_test_neg, DeepGoSE_AUROC = tool_auc,
            CAV_AUPR = aupr_val_vs_test_neg, DeepGoSE_AUPR = tool_aupr) |>
  pivot_longer(-c(ontology, go_term, families_per_100),
               names_to = c("method", "metric"), names_sep = "_") |>
  mutate(method = factor(method, levels = names(METHOD_COLORS)))

term_panel <- function(m) {
  term_long |> filter(metric == m) |>
    ggplot(aes(families_per_100, value, color = method)) +
    geom_point(size = 0.6, alpha = 0.35, stroke = 0) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.7, span = 0.9) +
    facet_wrap(~ontology, nrow = 1) +
    scale_color_manual(values = METHOD_COLORS, name = NULL) +
    scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    base_theme() +
    labs(x = "Sequence families per 100 CAV training positives", y = m) +
    theme(panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
          legend.position = "top")
}
p_a <- term_panel("AUROC")
p_b <- term_panel("AUPR") + theme(legend.position = "none")

ID_LEVELS <- c("no hit", "<30%", "30-50%", "50-90%", ">=90%")
pairs <- map_dfr(ONTS, \(o) {
  read_csv(sprintf("pair_novelty_%s.csv", o), show_col_types = FALSE) |>
    inner_join(read_csv(sprintf("%s/go_specificity_ranks_allterms_%s.csv", STAGE, o),
                        show_col_types = FALSE), by = c("protein_id", "go_term")) |>
    mutate(ontology = toupper(o))
}) |>
  mutate(ontology = factor(ontology, levels = toupper(ONTS)),
         id_bin = case_when(nn_identity == 0 ~ "no hit", nn_identity < 30 ~ "<30%",
                            nn_identity < 50 ~ "30-50%", nn_identity < 90 ~ "50-90%",
                            TRUE ~ ">=90%") |> factor(levels = ID_LEVELS))

tier <- function(rank, nopred) {
  case_when(nopred ~ "No prediction", rank == 1 ~ "Rank 1", rank <= 3 ~ "Rank 2-3",
            rank <= 10 ~ "Rank 4-10", TRUE ~ "Rank >10")
}
pair_long <- bind_rows(
  pairs |> transmute(ontology, id_bin, method = "CAV",
                     tier = tier(cav_rank_ont_filt, llr <= 0)),
  pairs |> transmute(ontology, id_bin, method = "DeepGoSE",
                     tier = tier(tool_rank_ont_filt, !as.logical(tool_predicted)))
) |> mutate(tier = factor(tier, levels = rev(RANK_LEVELS)),
            method = factor(method, levels = names(METHOD_COLORS)))
bin_n <- pairs |> count(ontology, id_bin)

p_c <- ggplot(pair_long, aes(id_bin, fill = tier)) +
  geom_bar(position = "fill", width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(data = bin_n, aes(id_bin, 1.04, label = n), inherit.aes = FALSE,
            size = 2.2, family = FIG_FONT, color = "gray30") +
  facet_grid(method ~ ontology) +
  scale_fill_manual(values = RANK_COLORS, breaks = RANK_LEVELS, name = NULL) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1),
                     limits = c(0, 1.08), expand = c(0, 0)) +
  base_theme() +
  labs(x = "Best identity of validation protein to the term's CAV training positives",
       y = "Share of pairs") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

fig <- plot_grid(p_a, p_b, p_c, ncol = 1, labels = c("A", "B", "C"), label_size = 8,
                 label_fontfamily = FIG_FONT, rel_heights = c(1.05, 0.9, 1.5))
ggsave("preview_fig2_diversity.png", fig, width = 7.2, height = 7.5, dpi = 300, bg = "white")
ggsave("preview_fig2_diversity.pdf", fig, width = 7.2, height = 7.5, device = cairo_pdf)

# console summary
cat("\nTerm-level Spearman rho (families_per_100 vs metric):\n")
term_long |> group_by(ontology, method, metric) |>
  summarise(rho = cor(families_per_100, value, method = "spearman", use = "complete.obs"),
            .groups = "drop") |> print(n = 20)
cat("\nPairs per identity bin and CAV / DeepGo-SE Rank-1 and top-10 (gated):\n")
pair_long |> group_by(ontology, id_bin, method) |>
  summarise(n = n(), rank1 = mean(tier == "Rank 1"),
            top10 = mean(tier %in% c("Rank 1", "Rank 2-3", "Rank 4-10")),
            nopred = mean(tier == "No prediction"), .groups = "drop") |>
  mutate(across(rank1:nopred, \(x) round(100 * x, 1))) |> print(n = 40)
