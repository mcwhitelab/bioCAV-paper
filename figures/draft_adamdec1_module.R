#!/usr/bin/env Rscript
# DRAFT / exploratory — checking the claimed fibroblast subpopulation module
# "uniquely expressing ADAMDEC1, CXCL14, EDNRB, and PROCR" against the
# fibroblast/colorectum CAV data. Not part of the figures.R pipeline; run
# standalone with: Rscript draft_adamdec1_module.R
#
# Data from single_cell/scripts/draft_adamdec1_module_check.py.

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggrepel)
  library(ggrastr)
})

RASTER_DPI <- 150
DATA <- "figure_data"
OUT  <- "figures"

oi <- c(vermillion = "#D55E00", blue = "#0072B2")
SC_NORMAL_COLOR <- "#3a6fad"
SC_CANCER_COLOR <- "#c0392b"
base_theme <- function(...) theme_cowplot(font_size = 8, ...)

scatter <- read_csv(file.path(DATA, "draft_adamdec1_scatter.csv"), show_col_types = FALSE)
corr    <- read_csv(file.path(DATA, "draft_adamdec1_gene_corr.csv"), show_col_types = FALSE)
cells   <- read_csv(file.path(DATA, "draft_adamdec1_cell_scores.csv"), show_col_types = FALSE)

# --- Panel A: DE log2FC vs CAV r for the 4 candidate genes ---
p_a <- ggplot(scatter, aes(log2fc, cav_r, label = gene_name)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(color = unname(oi["blue"]), size = 2.5) +
  ggrepel::geom_text_repel(size = 3, fontface = "bold", seed = 42) +
  base_theme() +
  labs(x = "log2FC (DE, mixedlm)", y = "Pearson r (CAV, paired donors)",
       title = "Claimed module: all DE-null, mixed CAV direction") +
  theme(plot.title = element_text(size = 8, face = "bold"))

# --- Panel B: pairwise co-expression among the 4 genes ---
gene_order <- c("ADAMDEC1", "CXCL14", "EDNRB", "PROCR")
p_b <- corr |>
  mutate(gene1 = factor(gene1, levels = gene_order),
         gene2 = factor(gene2, levels = rev(gene_order))) |>
  ggplot(aes(gene1, gene2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.6) +
  scale_fill_gradient2(low = unname(oi["vermillion"]), mid = "white", high = unname(oi["blue"]),
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  base_theme() +
  labs(x = NULL, y = NULL, title = "Pairwise co-expression (single cell)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 8, face = "bold"))

# --- Panel C: per-cell module score (z-mean of all 4 genes) vs L2 CAV axis ---
p_c <- cells |>
  mutate(disease = if_else(str_detect(tolower(disease), "normal"), "normal", "cancer"),
         disease = factor(disease, levels = c("normal", "cancer"))) |>
  ggplot(aes(l2_score, module_score, color = disease)) +
  geom_point_rast(size = 0.6, alpha = 0.4, raster.dpi = RASTER_DPI) +
  geom_smooth(aes(group = 1), method = "loess", color = "black", linewidth = 0.6, se = TRUE) +
  scale_color_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR), name = NULL) +
  base_theme() +
  labs(x = "L2 score  (← normal    cancer →)", y = "Module score (ADAMDEC1/CXCL14/EDNRB/PROCR)",
       title = "Module score vs. continuous CAV axis") +
  theme(legend.position = "top", plot.title = element_text(size = 8, face = "bold"))

# --- Panel D: ADAMDEC1/CXCL14 pair (strongest co-expressors) colored by disease ---
p_d <- cells |>
  mutate(disease = if_else(str_detect(tolower(disease), "normal"), "normal", "cancer"),
         disease = factor(disease, levels = c("normal", "cancer"))) |>
  ggplot(aes(expr_ADAMDEC1, expr_CXCL14, color = disease)) +
  geom_point_rast(size = 0.6, alpha = 0.4, raster.dpi = RASTER_DPI) +
  scale_color_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR), name = NULL) +
  base_theme() +
  labs(x = "ADAMDEC1 (log1p CPM)", y = "CXCL14 (log1p CPM)",
       title = "Strongest co-expressing pair (r=0.58)") +
  theme(legend.position = "top", plot.title = element_text(size = 8, face = "bold"))

fig <- plot_grid(p_a, p_b, p_c, p_d, nrow = 2, labels = "AUTO", label_size = 9,
                 align = "hv", axis = "tblr")

ggsave(file.path(OUT, "draft_fig_adamdec1_module.pdf"), fig,
       width = 8.5, height = 8, bg = "white", device = cairo_pdf)
ggsave(file.path(OUT, "draft_fig_adamdec1_module.png"), fig,
       width = 8.5, height = 8, bg = "white", dpi = 200)
message("Saved draft_fig_adamdec1_module.pdf/.png")
