#!/usr/bin/env Rscript
# draft_fig6_case_options.R -- PREVIEW / decision sheet only (not a paper figure).
# One card per usable single-cell pair, ordered by how well the L2 axis separates normal from
# cancer: top = L2 score distribution by condition (paired-donor cells; dashed = cells not used
# as CAV training positives), bottom = pseudobulk DESeq2 log2FC vs CAV gene r (DE genes in blue).
# Current Fig 6 cases are flagged. Numbers are shown because this is a decision aid.
# Reads: figure_data_staging/singlecell_cases/{case_metrics.tsv,case_cells.csv,case_de_vs_cav.csv}
#        (single_cell/convergence_check/14_allpair_separation.py; DE from 13a/13b)
# Writes: previews/fig6_case_options.{png,pdf}
suppressPackageStartupMessages({ library(tidyverse); library(cowplot); library(ggrastr) })
FIG_FONT <- "Arial"; NRM <- "#3a6fad"; CAN <- "#c0392b"
S <- "figure_data_staging/singlecell_cases"
met   <- read_tsv(file.path(S, "case_metrics.tsv"), show_col_types = FALSE)
cells <- read_csv(file.path(S, "case_cells.csv"), show_col_types = FALSE)
dev   <- read_csv(file.path(S, "case_de_vs_cav.csv"), show_col_types = FALSE)
CURRENT <- c("neutrophil__breast__normal_vs_breast_cancer", "epithelial_cell__lung10x__normal_vs_lung_cancer",
             "fibroblast__colorectum__normal_vs_colorectal_cancer")
met <- met |> mutate(sep = coalesce(auroc_heldout, auroc_all)) |> arrange(desc(sep))
th <- function() theme_cowplot(font_size = 6.5, font_family = FIG_FONT) +
  theme(axis.text = element_text(size = 5.5), axis.title = element_text(size = 6),
        axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
        plot.margin = margin(1, 3, 1, 2))
pretty <- function(p) p |> str_replace("__normal_vs_.*", "") |> str_replace("__", " / ") |>
  str_replace_all("_", " ") |> str_replace("lung10x", "lung (10x)")
card <- function(pid) {
  m <- filter(met, pair == pid); cc <- filter(cells, pair == pid)
  ho <- if (is.na(m$auroc_heldout)) sprintf("held-out n/a (%d cells, 1 class)", m$n_heldout) else sprintf("held-out %.2f", m$auroc_heldout)
  title <- paste0(pretty(pid), if (pid %in% CURRENT) "   [current Fig 6]" else "")
  sub1 <- sprintf("L2 AUROC %.2f, %s; random-dir. 95%% %.2f", m$auroc_all, ho, m$null95)
  sub2 <- sprintf("shift %.2f SD | donor AUROC %.2f-%.2f | %d donors, %d N / %d C",
                  m$shift_sd, m$donor_auroc_min, m$donor_auroc_max, m$n_donors, m$n_normal, m$n_cancer)
  top <- ggplot(cc, aes(l2_score, color = condition, fill = condition)) +
    geom_density(alpha = 0.15, linewidth = 0.4, adjust = 1.2) +
    geom_density(data = filter(cc, heldout), fill = NA, linewidth = 0.35, linetype = "22", adjust = 1.2) +
    scale_color_manual(values = c(normal = NRM, cancer = CAN)) + scale_fill_manual(values = c(normal = NRM, cancer = CAN)) +
    labs(x = "L2 score", y = "density", title = title, subtitle = paste(sub1, sub2, sep = "\n")) + th() +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(size = 6.8, face = if (pid %in% CURRENT) "bold" else "plain", margin = margin(b = 1)),
          plot.subtitle = element_text(size = 5.2, color = "grey25", margin = margin(b = 2), lineheight = 1.05))
  d <- filter(dev, pair == pid) |> mutate(sig = !is.na(padj) & padj < 0.05 & abs(log2fc_mle) >= 1)
  bot <- ggplot(d, aes(log2fc_mle, r)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.25) + geom_vline(xintercept = 0, color = "grey85", linewidth = 0.25) +
    geom_point_rast(data = filter(d, !sig), color = "grey65", size = 0.25, alpha = 0.3, raster.dpi = 200) +
    geom_point_rast(data = filter(d, sig), color = "#2277bb", size = 0.45, alpha = 0.7, raster.dpi = 200) +
    labs(x = "log2FC (pseudobulk)", y = "CAV r",
         subtitle = sprintf("Spearman %.2f over %s genes; %d DE genes", m$spearman_cav_de, format(m$n_de_tested, big.mark = ","), m$n_de_sig)) +
    th() + theme(plot.subtitle = element_text(size = 5.2, color = "grey25", margin = margin(b = 1)))
  plot_grid(top, bot, ncol = 1, rel_heights = c(1.15, 1), align = "v", axis = "lr")
}
cards <- map(met$pair, card)
key <- get_legend(ggplot(tibble(x = 1:2, c = factor(c("normal", "cancer"), levels = c("normal", "cancer"))), aes(x, x, color = c)) +
  geom_line(linewidth = 0.6) + scale_color_manual(values = c(normal = NRM, cancer = CAN), name = NULL) +
  theme_cowplot(font_size = 6.5, font_family = FIG_FONT) + theme(legend.position = "bottom"))
note <- ggdraw() + draw_label(paste(
  "Ordered by held-out L2 AUROC (cells not used as CAV training positives; in-sample AUROC where no held-out normal cells exist).",
  "Solid density = all paired-donor cells, dashed = held-out cells. 'random-dir. 95%' = 95th percentile AUROC of random directions",
  "orthogonalised against the same L0/L1 basis. Bottom: blue = DE genes (padj < 0.05, |log2FC| >= 1).", sep = "\n"),
  size = 6, fontfamily = FIG_FONT, x = 0.01, hjust = 0, lineheight = 1.1)
grid <- plot_grid(plotlist = cards, ncol = 5)
out <- plot_grid(grid, plot_grid(note, key, nrow = 1, rel_widths = c(4, 1)), ncol = 1, rel_heights = c(1, 0.06))
ggsave("previews/fig6_case_options.pdf", out, width = 13, height = 7.6, device = cairo_pdf)
ggsave("previews/fig6_case_options.png", out, width = 13, height = 7.6, dpi = 220, bg = "white")
message("Saved previews/fig6_case_options.{png,pdf}")
