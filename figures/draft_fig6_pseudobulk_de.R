#!/usr/bin/env Rscript
# draft_fig6_pseudobulk_de.R -- PREVIEW ONLY (not wired into figures.R).
# Figure 6 DE-vs-CAV scatter, second version: conventional pseudobulk DESeq2
# log2FC (paired design ~ donor + condition, apeglm-shrunk) alongside the
# current mixed-model version, for side-by-side review.
#
# Reads (read-only): figure_data/sc_de_vs_cav.csv, figure_data/sc_continuum_genes.csv,
#   single_cell/.../results/de_pseudobulk_paired_donors/<MIN>/<pair>.tsv,
#   single_cell/.../results/gene_correlation_paired_donors/<pair>.tsv
# Writes: previews/fig6_pseudobulk_de_preview_<MIN>.png
#
# Usage (from figures/): Rscript draft_fig6_pseudobulk_de.R [min10|min1] [mle|apeglm]

suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(ggrepel); library(ggrastr)
})

a    <- commandArgs(trailingOnly = TRUE)
MIN  <- if (length(a) >= 1) a[1] else "min10"
LFC  <- if (length(a) >= 2) a[2] else "mle"   # "mle" (DESeq2 log2FoldChange) or "apeglm" (shrunk)
FIG_FONT <- "Arial"
update_geom_defaults("text_repel", list(family = FIG_FONT))
RES  <- "../single_cell/cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d/results"
base_theme <- function(...) theme_cowplot(font_size = 8, font_family = FIG_FONT, ...) +
  theme(axis.text = element_text(size = 8), plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 7))

PAIRS <- c("neutrophil__breast__normal_vs_breast_cancer"         = "Neutrophil — breast",
           "epithelial_cell__lung10x__normal_vs_lung_cancer"     = "Epithelial cell — lung",
           "fibroblast__colorectum__normal_vs_colorectal_cancer" = "Fibroblast — colorectum")
SC_BACKGROUND_COLOR <- "#2277bb"; SC_CONTINUUM_FILL <- "#f0a500"; SC_CONTINUUM_EDGE <- "#7a5300"

mixed <- read_csv("figure_data/sc_de_vs_cav.csv", show_col_types = FALSE)
cont  <- read_csv("figure_data/sc_continuum_genes.csv", show_col_types = FALSE) %>% distinct(pair, gene)

scatter <- function(df, xlab, title, subtitle, sig_col, xthr) {
  df <- df %>% mutate(category = case_when(is_continuum ~ "continuum",
                                           .data[[sig_col]] ~ "sig", TRUE ~ "other"))
  s  <- df %>% filter(category == "sig") %>%
    mutate(d = sqrt((log2fc / max(abs(log2fc)))^2 + (cav_r / max(abs(cav_r)))^2)) %>%
    slice_max(d, n = 8)
  lab <- df %>% filter(category == "continuum" | gene %in% s$gene)
  ggplot(df, aes(log2fc, cav_r)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = c(-xthr, xthr), color = "grey75", linewidth = 0.25, linetype = "22") +
    geom_point_rast(data = filter(df, category == "other"), color = "grey70",
                    size = 0.4, alpha = 0.25, raster.dpi = 150) +
    geom_point_rast(data = filter(df, category == "sig"), color = SC_BACKGROUND_COLOR,
                    size = 0.8, alpha = 0.6, raster.dpi = 150) +
    geom_point(data = filter(df, category == "continuum"), shape = 21, fill = SC_CONTINUUM_FILL,
               color = SC_CONTINUUM_EDGE, size = 1.6, stroke = 0.3) +
    geom_text_repel(data = lab, aes(label = gene_name), size = 2, max.overlaps = 30,
                    segment.size = 0.25, segment.color = "grey60", color = "grey20",
                    box.padding = 0.35, min.segment.length = 0.1, seed = 42) +
    base_theme() + labs(x = xlab, y = "Pearson r (CAV)", title = title, subtitle = subtitle)
}

top <- list(); bot <- list()
for (p in names(PAIRS)) {
  cav <- read_tsv(file.path(RES, "gene_correlation_paired_donors", paste0(p, ".tsv")),
                  show_col_types = FALSE) %>% select(gene, cav_r = r)
  cids <- cont$gene[cont$pair == p]

  # Row 1: current mixed-model version (as in figures.R; |log2FC| >= 1.5 styling)
  m <- mixed %>% filter(pair == p) %>% mutate(is_continuum = gene %in% cids,
                                              above = abs(log2fc) >= 1.5)
  rm <- cor(m$log2fc, m$cav_r, method = "spearman")
  top[[p]] <- scatter(m, "Mixed-model effect (Δ mean log2 CP10K)", PAIRS[p],
                      sprintf("Spearman r = %.2f (n = %s)", rm, format(nrow(m), big.mark = ",")),
                      "above", 1.5)

  # Row 2: pseudobulk DESeq2, apeglm-shrunk log2FC; sig = padj < 0.05 & |LFC| >= 1
  pb <- read_tsv(file.path(RES, "de_pseudobulk_paired_donors", MIN, paste0(p, ".tsv")),
                 show_col_types = FALSE) %>%
    inner_join(cav, by = "gene") %>%
    mutate(log2fc_apeglm = log2fc, log2fc = if (LFC == "mle") log2fc_mle else log2fc_apeglm) %>%
    filter(is.finite(log2fc), is.finite(cav_r)) %>%
    mutate(is_continuum = gene %in% cids,
           sig = !is.na(padj) & padj < 0.05 & abs(log2fc) >= 1)
  rp  <- cor(pb$log2fc, pb$cav_r, method = "spearman")
  ralt <- cor(if (LFC == "mle") pb$log2fc_apeglm else pb$log2fc_mle, pb$cav_r, method = "spearman")
  same <- m %>% filter(gene %in% pb$gene)
  rms <- cor(same$log2fc, same$cav_r, method = "spearman")
  bot[[p]] <- scatter(pb, sprintf("log2FC (pseudobulk DESeq2%s)", if (LFC == "mle") "" else ", apeglm"), PAIRS[p],
                      sprintf("Spearman r = %.2f (n = %s genes)\n%s LFC: r = %.2f; %d donors; %d DE genes\nmixed model, same genes: r = %.2f",
                              rp, format(nrow(pb), big.mark = ","), if (LFC == "mle") "shrunk" else "MLE", ralt,
                              pb$n_donors[1], sum(pb$sig), rms),
                      "sig", 1)
}

row_lab <- function(txt) ggdraw() + draw_label(txt, fontface = "bold", size = 8, x = 0.01, hjust = 0,
                                              fontfamily = FIG_FONT)
fig <- plot_grid(
  row_lab("Current: per-cell mixed model (donor random intercept); highlighted |effect| ≥ 1.5"),
  plot_grid(plotlist = top, nrow = 1, labels = c("B", "C", "D"), label_size = 8),
  row_lab(sprintf("Proposed: pseudobulk DESeq2 ~ donor + condition (paired), %s LFC; highlighted padj < 0.05 & |log2FC| ≥ 1  [%s cells/sample]",
                  if (LFC == "mle") "MLE" else "apeglm", sub("min", "≥", MIN))),
  plot_grid(plotlist = bot, nrow = 1, labels = c("B'", "C'", "D'"), label_size = 8),
  ncol = 1, rel_heights = c(0.07, 1, 0.07, 1))

dir.create("previews", showWarnings = FALSE)
out <- sprintf("previews/fig6_pseudobulk_de_preview_%s_%s.png", MIN, LFC)
ggsave(out, fig, width = 7.2, height = 5.6, dpi = 300, bg = "white")
cat("wrote", out, "\n")
