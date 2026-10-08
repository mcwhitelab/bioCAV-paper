# Exploratory preview (2026-09-29): where does rank signal drop with the number of
# sequence families among a GO term's CAV training positives?
# Data: family_cutoff_screen.py (strict _pf filter, gated; BLAST tie-averaged).
#   A  top-10 rate per family-count bin (non-cumulative), 95% term-bootstrap CI
#   B  Rank-1 rate, same bins
#   C  top-10 below (solid) vs at-or-above (dashed) each cutoff K
# Run from bioCAV-paper/go_eval/blast_baseline/.
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

FIG_FONT <- "Arial"
METHOD_COLORS <- c(CAV = "#0072B2", DeepGoSE = "#D55E00", BLAST = "#6E6E6E")
BIN_LABELS <- c("1-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40-59", "60-99", ">=100")
ONT_LEVELS <- c("MF", "BP", "CC", "All")
base_theme <- function() {
  theme_cowplot(font_size = 8, font_family = FIG_FONT) +
    theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8),
          strip.text = element_text(size = 8), strip.background = element_blank(),
          panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3))
}

bins <- read_csv("family_screen_bins.csv", show_col_types = FALSE) |>
  mutate(bin = factor(bin, levels = BIN_LABELS), ontology = factor(ontology, levels = ONT_LEVELS),
         method = factor(method, levels = names(METHOD_COLORS)))
cut <- read_csv("family_screen_cutoff.csv", show_col_types = FALSE) |>
  mutate(ontology = factor(ontology, levels = ONT_LEVELS),
         method = factor(method, levels = names(METHOD_COLORS)),
         side = factor(side, levels = c("below", "at_or_above"),
                       labels = c("< K families", ">= K families")))

n_lab <- bins |> distinct(ontology, bin, n_terms)
dodge <- position_dodge(width = 0.5)

bin_panel <- function(m, ylab) {
  ggplot(filter(bins, metric == m), aes(bin, rate, color = method, group = method)) +
    geom_line(position = dodge, linewidth = 0.4) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = dodge, size = 0.15, linewidth = 0.3) +
    facet_wrap(~ontology, nrow = 1) +
    scale_color_manual(values = METHOD_COLORS, name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    base_theme() +
    labs(x = "Sequence families among the term's CAV training positives", y = ylab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
p_a <- bin_panel("top10", "Pairs in top 10") +
  geom_text(data = n_lab, aes(bin, 1, label = n_terms), inherit.aes = FALSE,
            size = 1.9, family = FIG_FONT, color = "gray40", vjust = 1) +
  theme(legend.position = "top")
p_b <- bin_panel("rank1", "Pairs at Rank 1") + theme(legend.position = "none")

p_c <- ggplot(filter(cut, metric == "top10"),
              aes(cutoff, rate, color = method, linetype = side, fill = method)) +
  geom_ribbon(data = filter(cut, metric == "top10", side == "< K families"),
              aes(ymin = lo, ymax = hi), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.6) +
  facet_wrap(~ontology, nrow = 1) +
  scale_color_manual(values = METHOD_COLORS, name = NULL) +
  scale_fill_manual(values = METHOD_COLORS, guide = "none") +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_x_continuous(breaks = c(10, 20, 40, 60, 100)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  base_theme() +
  labs(x = "Cutoff K (sequence families)", y = "Pairs in top 10") +
  theme(legend.position = "top") +
  guides(color = "none")

fig <- plot_grid(p_a, p_b, p_c, ncol = 1, labels = c("A", "B", "C"), label_size = 8,
                 label_fontfamily = FIG_FONT, rel_heights = c(1.15, 1, 1.05))
ggsave("preview_family_screen.png", fig, width = 7.2, height = 7.2, dpi = 300, bg = "white")
ggsave("preview_family_screen.pdf", fig, width = 7.2, height = 7.2, device = cairo_pdf)
