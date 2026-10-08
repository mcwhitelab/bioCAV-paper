#!/usr/bin/env Rscript
# draft_fig6_fib_modules.R -- PREVIEW ONLY (not wired into figures.R).
# Compact candidate replacement for the bottom half of Fig 6D (fibroblast continuum strips):
# module score vs L2 CAV score for the two fibroblast programs found in the top/bottom-20
# co-expression heatmap (previews/fib_coexpr_top20_stored.*): ECM / quiescent (9 annotated genes) and
# ADAMDEC1+ mucosal (6 genes), over 2000 paired-donor fibroblasts. Lines = per-condition mean
# module score in L2 deciles, ribbons = 95% CI, so the within-condition trend and the
# condition overlap are both visible. No statistics on the panel; they are in METHODS §6.4 / sc_fib_module_stats.tsv.
# Reads: figure_data_staging/singlecell_modules/sc_fib_module_scores.csv
#        (single_cell/convergence_check/10_export_module_scores.py)
# Writes: previews/fig6D_fib_modules.{png,pdf}  (2.5 x 2.2 in, the slot size)
suppressPackageStartupMessages({ library(tidyverse); library(cowplot); library(ggrastr) })
FIG_FONT <- "Arial"
SC_NORMAL_COLOR <- "#3a6fad"; SC_CANCER_COLOR <- "#c0392b"
base_theme <- function(...) theme_cowplot(font_size = 7, font_family = FIG_FONT, ...) +
  theme(axis.text = element_text(size = 6), axis.title = element_text(size = 6.5),
        plot.title = element_text(size = 6.5, face = "plain", hjust = 0, margin = margin(b = 0.5)),
        plot.subtitle = element_text(size = 5.5, face = "italic", color = "grey30", hjust = 0, margin = margin(b = 1.5)),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(2, 3, 1, 2))
d <- read_csv("figure_data_staging/singlecell_modules/sc_fib_module_scores.csv", show_col_types = FALSE) |>
  mutate(condition = factor(condition, levels = c("normal", "cancer")))
# gene lists for the subtitles, from the same stats file the scores came from
gene_list <- read_tsv("figure_data_staging/singlecell_modules/sc_fib_module_stats.tsv", show_col_types = FALSE) |>
  filter(cells == "all") |> mutate(col = paste0(module, "_score"), genes = str_replace_all(genes, ",", ", ")) |>
  select(col, genes) |> mutate(genes = str_wrap(genes, 34)) |> deframe()
MODS <- c(ecm_score = "ECM / quiescent program (9 genes)",
          adamdec1_score = "ADAMDEC1+ mucosal program (6 genes)")
panel <- function(col, show_x) {
  dd <- d |> transmute(condition, l2_score, score = .data[[col]])
  bins <- dd |> group_by(condition) |> mutate(bin = ntile(l2_score, 10)) |>
    group_by(condition, bin) |>
    summarise(x = mean(l2_score), y = mean(score), se = sd(score) / sqrt(n()), .groups = "drop")
  ggplot(bins, aes(x, y, color = condition, fill = condition)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.25) +
    geom_ribbon(aes(ymin = y - 1.96 * se, ymax = y + 1.96 * se), color = NA, alpha = 0.18, show.legend = FALSE) +
    geom_line(linewidth = 0.6, key_glyph = "path") +
    geom_point(size = 0.8, show.legend = FALSE) +
    scale_color_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR), name = NULL) +
    scale_fill_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR), name = NULL) +
    scale_x_continuous(limits = range(bins$x) + c(-0.15, 0.15)) +
    labs(title = MODS[[col]], subtitle = gene_list[[col]], y = "module score (z)",
         x = if (show_x) "L2 score (cancer-CAV direction)" else NULL) +
    base_theme() +
    theme(legend.position = if (show_x) "none" else c(0.78, 0.85),
          legend.text = element_text(size = 6), legend.key.size = unit(6, "pt"),
          legend.background = element_blank(),
          axis.text.x = if (show_x) element_text(size = 6) else element_blank(),
          axis.ticks.x = if (show_x) element_line(linewidth = 0.3) else element_blank())
}
p <- plot_grid(panel("ecm_score", FALSE), panel("adamdec1_score", TRUE),
               ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1.2))
ggsave("previews/fig6D_fib_modules.pdf", p, width = 2.5, height = 2.2, device = cairo_pdf)
ggsave("previews/fig6D_fib_modules.png", p, width = 2.5, height = 2.2, dpi = 400, bg = "white")
message("Saved previews/fig6D_fib_modules.{pdf,png}")
