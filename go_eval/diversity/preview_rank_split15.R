# Exploratory preview (2026-09-29): Fig 2C in the figures.R style, split by
# whether the GO term's CAV training positives span <15 or >=15 sequence
# families (term_diversity.py). Known-positive-excluded data; gated ranks.
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

pairs <- map_dfr(ONTS, \(o) {
  read_csv(sprintf("%s/go_specificity_ranks_allterms_%s.csv", STAGE, o), show_col_types = FALSE) |>
    inner_join(read_csv(sprintf("term_diversity_%s.csv", o), show_col_types = FALSE),
               by = "go_term") |>
    mutate(ontology = toupper(o))
}) |>
  mutate(split = if_else(n_families < 15, "<15 families", ">=15 families"))

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

ggsave("preview_rank_split15.png", fig, width = 7.2, height = 2.9, dpi = 300, bg = "white")
ggsave("preview_rank_split15.pdf", fig, width = 7.2, height = 2.9, device = cairo_pdf)

print(counts)
rank_df |> group_by(split, ontology, method) |>
  summarise(rank1 = round(100 * sum(prop[tier == "Rank 1"]), 1),
            top10 = round(100 * sum(prop[tier %in% RANK_LEVELS[1:3]]), 1),
            nopred = round(100 * sum(prop[tier == "No prediction"]), 1), .groups = "drop") |>
  arrange(split, ontology, method) |> print(n = 20)
