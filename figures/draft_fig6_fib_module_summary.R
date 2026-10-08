#!/usr/bin/env Rscript
# draft_fig6_fib_module_summary.R -- PREVIEW ONLY (not wired into figures.R).
# Two compact summaries of the fibroblast top/bottom-20 co-expression heatmap (current CAVs;
# full matrix = previews/fib_coexpr_top20_stored.*), sized for the lower half of Fig 6D:
#   A  module x module mean co-expression (diagonal = within-module, excluding self)
#   B  'little clusters': each module's own gene x gene block, 2 x 2
# ECM module = the 9 annotated ECM genes. Labels coloured by CAV direction (red = up in cancer,
# blue = down); colour scale = Pearson co-expression (purple - white - green), as in the heatmap.
# Reads: figure_data_staging/singlecell_modules/sc_fib_module_{corr,genes}.csv
#        (single_cell/convergence_check/11_export_module_corr.py)
# Writes: previews/fig6D_module_summary_{A_matrix,B_clusters}.{png,pdf}  (2.5 x 2.2 in)
suppressPackageStartupMessages({ library(tidyverse); library(cowplot) })
FIG_FONT <- "Arial"
UP <- "#c0392b"; DOWN <- "#3a6fad"
fill_scale <- scale_fill_gradient2(low = "#762a83", mid = "white", high = "#1b7837", midpoint = 0,
                                   limits = c(-1, 1), breaks = c(-1, 0, 1), name = "co-expression r")
base_theme <- function(...) theme_cowplot(font_size = 7, font_family = FIG_FONT, ...) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), plot.margin = margin(2, 2, 2, 2))
S <- "figure_data_staging/singlecell_modules"
genes <- read_csv(file.path(S, "sc_fib_module_genes.csv"), show_col_types = FALSE)
corr  <- read_csv(file.path(S, "sc_fib_module_corr.csv"), show_col_types = FALSE)
SHORT <- c("ECM / quiescent fibroblast" = "ECM / quiescent", "ADAMDEC1+ mucosal fibroblast" = "ADAMDEC1+",
           "T cell" = "T cell", "B / plasma" = "B / plasma")
MOD_ORDER <- c("ECM / quiescent fibroblast", "ADAMDEC1+ mucosal fibroblast", "T cell", "B / plasma")
mods <- genes |> group_by(module) |> summarise(n = n(), cav = mean(cav_r), .groups = "drop") |>
  mutate(module = factor(module, levels = MOD_ORDER), lab = sprintf("%s (%d)", SHORT[as.character(module)], n),
         col = if_else(cav > 0, UP, DOWN)) |> arrange(module)

# ---- A: module x module mean co-expression ----
cm <- corr |> left_join(select(genes, gene1 = gene, m1 = module), by = "gene1") |>
  left_join(select(genes, gene2 = gene, m2 = module), by = "gene2") |>
  filter(gene1 != gene2) |> group_by(m1, m2) |> summarise(r = mean(r), .groups = "drop") |>
  mutate(m1 = factor(m1, levels = MOD_ORDER), m2 = factor(m2, levels = rev(MOD_ORDER)))
labA <- setNames(mods$lab, as.character(mods$module)); colA <- setNames(mods$col, as.character(mods$module))
pA <- ggplot(cm, aes(m1, m2, fill = r)) + geom_tile(color = "white", linewidth = 0.8) + fill_scale +
  scale_x_discrete(labels = labA, position = "top") + scale_y_discrete(labels = labA) + coord_equal() +
  base_theme() + labs(x = NULL, y = NULL) +
  theme(plot.margin = margin(2, 14, 2, 2),
        axis.text.x = element_text(size = 6, angle = 30, hjust = 0, color = colA[levels(cm$m1)]),
        axis.text.y = element_text(size = 6, color = colA[rev(levels(cm$m1))]),
        legend.position = "bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
        legend.key.height = unit(4, "pt"), legend.key.width = unit(12, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.title.position = "top", legend.title.align = 0.5)
ggsave("previews/fig6D_module_summary_A_matrix.pdf", pA, width = 2.5, height = 2.2, device = cairo_pdf)
ggsave("previews/fig6D_module_summary_A_matrix.png", pA, width = 2.5, height = 2.2, dpi = 400, bg = "white")

# ---- B: little clusters (each module's own block) ----
block <- function(mod) {
  g <- genes |> filter(module == mod) |> arrange(order)
  d <- corr |> filter(gene1 %in% g$gene, gene2 %in% g$gene) |>
    mutate(gene1 = factor(gene1, levels = g$gene), gene2 = factor(gene2, levels = rev(g$gene)))
  gc <- setNames(if_else(g$cav_r > 0, UP, DOWN), g$gene)
  ggplot(d, aes(gene1, gene2, fill = r)) + geom_tile(color = "white", linewidth = 0.3) + fill_scale +
    coord_equal() + base_theme() + labs(x = NULL, y = NULL, title = SHORT[[mod]]) +
    theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 5, color = gc[rev(levels(d$gene2))]),
          plot.title = element_text(size = 6.5, face = "plain", hjust = 0.5, margin = margin(b = 1)),
          legend.position = "none")
}
legB <- get_legend(pA + theme(legend.position = "bottom"))
pB <- plot_grid(plot_grid(block(MOD_ORDER[1]), block(MOD_ORDER[2]), block(MOD_ORDER[3]), block(MOD_ORDER[4]),
                          ncol = 2, rel_heights = c(9, 7)),
                legB, ncol = 1, rel_heights = c(1, 0.13))
ggsave("previews/fig6D_module_summary_B_clusters.pdf", pB, width = 2.5, height = 2.2, device = cairo_pdf)
ggsave("previews/fig6D_module_summary_B_clusters.png", pB, width = 2.5, height = 2.2, dpi = 400, bg = "white")

# ---- C: shrunken heatmap of the module genes only (block structure + between-module contrast) ----
gC <- genes |> mutate(module = factor(module, levels = MOD_ORDER)) |> arrange(module, order) |> mutate(pos = row_number())
dC <- corr |> inner_join(select(gC, gene1 = gene, x = pos), by = "gene1") |>
  inner_join(select(gC, gene2 = gene, y = pos), by = "gene2")
bC <- gC |> group_by(module) |> summarise(a = min(pos), b = max(pos), cav = mean(cav_r), .groups = "drop") |>
  mutate(lab = SHORT[as.character(module)], col = if_else(cav > 0, UP, DOWN), mid = (a + b) / 2)
nC <- nrow(gC)
pC <- ggplot(dC, aes(x, y, fill = r)) + geom_tile() + fill_scale +
  geom_vline(xintercept = bC$a[-1] - 0.5, color = "white", linewidth = 0.6) +
  geom_hline(yintercept = bC$a[-1] - 0.5, color = "white", linewidth = 0.6) +
  # module names to the right of the blocks, coloured by CAV direction
  geom_text(data = bC, aes(x = nC + 1.2, y = mid, label = lab, color = I(col)), inherit.aes = FALSE,
            hjust = 0, size = 6 / .pt, family = FIG_FONT) +
  scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  coord_equal(clip = "off") + base_theme() + labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(), plot.margin = margin(2, 52, 2, 6),
        legend.position = "bottom", legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
        legend.key.height = unit(4, "pt"), legend.key.width = unit(12, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.title.position = "top")
ggsave("previews/fig6D_module_summary_C_blocks.pdf", pC, width = 2.5, height = 2.2, device = cairo_pdf)
ggsave("previews/fig6D_module_summary_C_blocks.png", pC, width = 2.5, height = 2.2, dpi = 400, bg = "white")
message("Saved previews/fig6D_module_summary_{A_matrix,B_clusters,C_blocks}.{png,pdf}")
