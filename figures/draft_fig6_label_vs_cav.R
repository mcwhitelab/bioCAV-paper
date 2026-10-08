#!/usr/bin/env Rscript
# draft_fig6_label_vs_cav.R -- PREVIEW ONLY (not wired into figures.R).
# Figure 6 case-study pairs: per-gene correlation with the condition label
# ("label r": Pearson r across cells between 0/1 normal/cancer and log2 CP10K)
# vs correlation with the L2 CAV score (CAV r). Same cells (paired donors),
# same normalisation, same genes. Genes on the diagonal are ones the CAV axis
# recovers only as "cancer vs normal"; genes far off it carry per-cell
# structure the label does not.
#
# Reads:  single_cell/.../results/label_correlation_paired_donors/<pair>.tsv
#         (scripts/paired_donor_label_corr.py), figure_data/sc_continuum_genes.csv
# Writes: previews/fig6_label_vs_cav_fibroblast.png (single panel)
#         previews/fig6_label_vs_cav_3pairs.png     (all three side by side)
#
# Usage (from figures/): Rscript draft_fig6_label_vs_cav.R

suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(ggrepel); library(ggrastr)
})

FIG_FONT <- "Arial"
update_geom_defaults("text_repel", list(family = FIG_FONT))
RES  <- "../single_cell/cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d/results/label_correlation_paired_donors"
SC_BACKGROUND_COLOR <- "#2277bb"; SC_CONTINUUM_FILL <- "#f0a500"; SC_CONTINUUM_EDGE <- "#7a5300"
PAIRS <- c("neutrophil__breast__normal_vs_breast_cancer"         = "Neutrophil — breast",
           "epithelial_cell__lung10x__normal_vs_lung_cancer"     = "Epithelial cell — lung",
           "fibroblast__colorectum__normal_vs_colorectal_cancer" = "Fibroblast — colorectum")

read_pair <- function(p) read_tsv(file.path(RES, paste0(p, ".tsv")), show_col_types = FALSE)
rho_of    <- function(d) cor(d$label_r, d$cav_r, method = "spearman")
cont      <- read_csv("figure_data/sc_continuum_genes.csv", show_col_types = FALSE)

build_panel <- function(pair, key = FALSE, subtitle_extra = NULL) {
  df   <- read_pair(pair)
  cids <- cont %>% filter(pair == !!pair) %>% pull(gene) %>% unique()
  fit  <- lm(cav_r ~ label_r, data = df)
  df   <- df %>% mutate(resid = cav_r - predict(fit, df), is_continuum = gene %in% cids)

  # Label: the continuum genes, plus the largest residuals from the fitted
  # CAV r ~ label r trend on each side (CAV-correlated beyond what the label
  # explains).
  off <- bind_rows(df %>% filter(!is_continuum) %>% slice_min(resid, n = 5),
                   df %>% filter(!is_continuum) %>% slice_max(resid, n = 5))
  df  <- df %>% mutate(category = case_when(is_continuum ~ "continuum",
                                            gene %in% off$gene ~ "off", TRUE ~ "other"))
  lab <- df %>% filter(category != "other")

  lim <- range(c(df$label_r, df$cav_r)) * 1.05
  span <- diff(lim)
  p <- ggplot(df, aes(label_r, cav_r)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "grey55", linewidth = 0.3, linetype = "22") +
    geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1], color = "grey40", linewidth = 0.3) +
    geom_point_rast(data = filter(df, category == "other"), color = "grey70",
                    size = 0.4, alpha = 0.25, raster.dpi = 300) +
    geom_point(data = filter(df, category == "off"), color = SC_BACKGROUND_COLOR,
               size = 1.1, alpha = 0.85) +
    geom_point(data = filter(df, category == "continuum"), shape = 21, fill = SC_CONTINUUM_FILL,
               color = SC_CONTINUUM_EDGE, size = 1.8, stroke = 0.3) +
    geom_text_repel(data = lab, aes(label = gene_name), size = 2, family = FIG_FONT,
                    color = "grey20", max.overlaps = 40, segment.size = 0.25,
                    segment.color = "grey60", box.padding = 0.3, min.segment.length = 0.05,
                    force = 3, seed = 42) +
    annotate("text", x = lim[2], y = lim[2], label = "y = x", hjust = 1.1, vjust = -0.4,
             size = 2, color = "grey45", family = FIG_FONT) +
    coord_equal(xlim = lim, ylim = lim) +
    theme_cowplot(font_size = 8, font_family = FIG_FONT) +
    theme(plot.title = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 7)) +
    labs(x = "Pearson r with condition label (normal = 0, cancer = 1)",
         y = "Pearson r with L2 CAV score",
         title = PAIRS[pair],
         subtitle = paste0(sprintf("Spearman r = %.2f (n = %s genes)", rho_of(df),
                                   format(nrow(df), big.mark = ",")), subtitle_extra))
  if (key) {
    # Key in the (empty) upper-left corner, positioned relative to the limits
    kx <- lim[1] + 0.04 * span; ky <- lim[2] - c(0.05, 0.11, 0.17) * span
    p <- p +
      annotate("point", x = kx, y = ky[1:2], shape = 21, size = c(1.8, 1.1), stroke = 0.3,
               fill = c(SC_CONTINUUM_FILL, SC_BACKGROUND_COLOR),
               color = c(SC_CONTINUUM_EDGE, SC_BACKGROUND_COLOR)) +
      annotate("segment", x = kx - 0.02 * span, xend = kx + 0.02 * span, y = ky[3], yend = ky[3],
               color = "grey40", linewidth = 0.3) +
      annotate("text", x = kx + 0.04 * span, y = ky, hjust = 0, size = 2, color = "grey20",
               family = FIG_FONT,
               label = c("continuum-panel gene", "largest residual from trend", "fitted trend"))
  }
  p
}

dir.create("previews", showWarnings = FALSE)

# Single fibroblast panel (as previewed first)
rhos <- sapply(names(PAIRS), function(p) rho_of(read_pair(p)))
p_fib <- build_panel(names(PAIRS)[3], key = TRUE,
                     subtitle_extra = sprintf("\nfor comparison: neutrophil %.2f, lung epithelial %.2f",
                                              rhos[1], rhos[2]))
ggsave("previews/fig6_label_vs_cav_fibroblast.png", p_fib, width = 3.6, height = 3.9, dpi = 300, bg = "white")
cat("wrote previews/fig6_label_vs_cav_fibroblast.png\n")

# All three side by side
panels <- lapply(seq_along(PAIRS), function(i) build_panel(names(PAIRS)[i], key = (i == 1)))
fig <- plot_grid(plotlist = panels, nrow = 1, align = "h", labels = c("A", "B", "C"), label_size = 8)
ggsave("previews/fig6_label_vs_cav_3pairs.png", fig, width = 10, height = 3.8, dpi = 300, bg = "white")
cat("wrote previews/fig6_label_vs_cav_3pairs.png\n")
