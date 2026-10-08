#!/usr/bin/env Rscript
# draft_supp_fib_heatmaps.R -- PREVIEW ONLY (not wired into figures.R).
# Supplementary figure: fibroblast/colorectum, top/bottom 20 genes by paired-donor CAV r,
# clustered by Pearson co-expression across the 2,000 paired-donor fibroblasts.
#   A  C = 1, PCA-128, converged (convergence_check arm conv1; the Fig 6 axis)
#   B  same genes ranked by pseudobulk DESeq2 Wald statistic (min10; de_stat), CAV r from conv1
#   C  same frame, C = 1e-4 (shared scaler, PCA-128, converged; convergence_check arm conv4)
# Gene order, module calls and the weak-member rule come from
# single_cell/convergence_check/08_fib_module_heatmap.py (average linkage on 1 - r, cut at
# 0.85; a cluster is named when >= 2 of its genes hit a marker set; a gene whose mean
# co-expression with its module is < 0.5 x the module median is moved out unlabelled).
# Gene labels: gene (CAV r) or, in B, gene (Wald; CAV r); red = up in cancer, blue = down (sign of the
# ranking statistic, so in B it is the DE direction). Fill: co-expression r.
# Reads: figure_data_staging/singlecell_heatmaps/fib_coexpr_{conv1,de_stat,conv4}_{corr,genes}.csv
# Writes: previews/supp_fib_heatmaps.{png,pdf}
suppressPackageStartupMessages({ library(tidyverse); library(cowplot) })
FIG_FONT <- "Arial"; UP <- "#c0392b"; DOWN <- "#3a6fad"
S <- "figure_data_staging/singlecell_heatmaps"
fill_scale <- scale_fill_gradient2(low = "#762a83", mid = "white", high = "#1b7837", midpoint = 0,
                                   limits = c(-1, 1), breaks = c(-1, 0, 1), name = "co-expression r")
heat <- function(arm, title) {
  g <- read_csv(file.path(S, sprintf("fib_coexpr_%s_genes.csv", arm)), show_col_types = FALSE) |>
    mutate(module = replace_na(module, ""),
           # short module names, matching Fig 6D
           module = recode(module, "ECM / quiescent fibroblast" = "ECM / quiescent",
                           "ADAMDEC1+ mucosal fibroblast" = "ADAMDEC1+ mucosal"))
  d <- read_csv(file.path(S, sprintf("fib_coexpr_%s_corr.csv", arm)), show_col_types = FALSE) |>
    left_join(select(g, gene1 = gene, x = pos), by = "gene1") |>
    left_join(select(g, gene2 = gene, y = pos), by = "gene2")
  n <- nrow(g)
  blocks <- g |> mutate(run = cumsum(module != lag(module, default = "<start>"))) |>
    group_by(run, module) |> summarise(a = min(pos), b = max(pos), .groups = "drop") |>
    mutate(w = b - a + 1, mid = (a + b) / 2,
           lab = map2_chr(module, w, \(m, w) if (nchar(m) <= w * 2.4) m else str_wrap(m, max(8, round(w * 2.4)))))
  named <- filter(blocks, module != "")
  g <- g |> mutate(lab = if_else(key == "CAV r", sprintf("%s (%+.2f)", gene, cav_r),
                                 sprintf("%s (%+.1f; %+.2f)", gene, value, cav_r)),
                  col = if_else(value > 0, UP, DOWN))
  ggplot(d, aes(x, y, fill = r)) +
    geom_tile() + fill_scale +
    geom_vline(xintercept = blocks$a[-1] - 0.5, color = "white", linewidth = 0.7) +
    geom_hline(yintercept = blocks$a[-1] - 0.5, color = "white", linewidth = 0.7) +
    # module brackets + names above the blocks
    geom_segment(data = named, aes(x = a - 0.35, xend = b + 0.35, y = -0.4, yend = -0.4),
                 inherit.aes = FALSE, color = "#666666", linewidth = 0.35) +
    geom_text(data = named, aes(x = mid, y = -0.9, label = lab), inherit.aes = FALSE,
              vjust = 0, size = 6 / .pt, lineheight = 0.9, family = FIG_FONT, color = "#333333") +
    # gene labels on the right, coloured by CAV direction
    geom_text(data = g, aes(x = n + 0.9, y = pos, label = lab, color = I(col)), inherit.aes = FALSE,
              hjust = 0, size = 5.2 / .pt, family = FIG_FONT) +
    scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
    coord_equal(clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_family = FIG_FONT) +
    theme(plot.title = element_text(size = 7.5, face = "bold", hjust = 0, margin = margin(b = 24)),
          legend.position = "none", plot.margin = margin(4, 70, 4, 4))
}
pA <- heat("conv1", "C = 1, PCA-128, converged (the Fig 6 axis); label = gene (CAV r)")
pB <- heat("de_stat", "Pseudobulk DE Wald statistic; label = gene (Wald; CAV r, C = 1)")
pC <- heat("conv4", "C = 10⁻⁴, PCA-128, converged; label = gene (CAV r)")
leg <- get_legend(pA + theme(legend.position = "bottom", legend.title = element_text(size = 6.5),
                             legend.text = element_text(size = 6), legend.key.height = unit(5, "pt"),
                             legend.key.width = unit(22, "pt"), legend.title.position = "top"))
fig <- plot_grid(pA, pB, pC, leg, ncol = 1, rel_heights = c(1, 1, 1, 0.1),
                 labels = c("A", "B", "C", ""), label_size = 9, label_fontfamily = FIG_FONT)
ggsave("previews/supp_fib_heatmaps.pdf", fig, width = 6.5, height = 14, device = cairo_pdf)
ggsave("previews/supp_fib_heatmaps.png", fig, width = 6.5, height = 14, dpi = 300, bg = "white")
message("Saved previews/supp_fib_heatmaps.{pdf,png}")
