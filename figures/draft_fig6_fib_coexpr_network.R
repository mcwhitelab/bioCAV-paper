#!/usr/bin/env Rscript
# draft_fig6_fib_coexpr_network.R -- PREVIEW ONLY (not wired into figures.R).
# Network view of the top/bottom-20 fibroblast/colorectum co-expression heatmap
# (previews/fib_coexpr_top20_stored.*): nodes = the two fibroblast programs of the Fig 6D module
# panel (ECM9 known ECM/quiescent markers + ADAMDEC1+ mucosal, read from sc_fib_module_stats.tsv)
# plus the heatmap's T cell and B / plasma clusters; the 11 other ECM-block genes, RPL39 and
# NSG1 are dropped. Edges = Pearson co-expression
# r >= EDGE_R across the 2000 paired-donor fibroblasts, igraph/graphlayouts stress layout
# on edge length 1 - r. Node fill = module (same colours as the heatmap module bar), hull = module,
# label = gene name, label colour = CAV direction (red up in cancer, blue down).
# Reads: figure_data_staging/singlecell_heatmaps/fib_coexpr_stored_{corr,genes}.csv
#        figure_data_staging/singlecell_modules/sc_fib_module_stats.tsv
#        (single_cell/convergence_check/08_fib_module_heatmap.py)
# Writes: previews/fib_coexpr_network_stored.{png,pdf}
# Run: R_LIBS_USER=~/R/4.5 Rscript draft_fig6_fib_coexpr_network.R
suppressPackageStartupMessages({ library(tidyverse); library(igraph); library(ggraph); library(ggforce) })
FIG_FONT <- "Arial"
EDGE_R <- 0.25
SC_NORMAL_COLOR <- "#3a6fad"; SC_CANCER_COLOR <- "#c0392b"
MOD_COLOR <- c("T cell" = "#D62728", "B / plasma" = "#E377C2",
               "ECM / quiescent fibroblast" = "#8C564B", "ADAMDEC1+ mucosal fibroblast" = "#E6C229")
D <- "figure_data_staging/singlecell_heatmaps"
fig6_genes <- read_tsv("figure_data_staging/singlecell_modules/sc_fib_module_stats.tsv", show_col_types = FALSE) |>
  filter(cells == "all") |> pull(genes) |> str_split(",") |> unlist()
genes <- read_csv(file.path(D, "fib_coexpr_stored_genes.csv"), show_col_types = FALSE) |>
  filter(gene %in% fig6_genes | module %in% c("T cell", "B / plasma")) |>
  transmute(name = gene, module = factor(module, levels = names(MOD_COLOR)), cav_r)
edges <- read_csv(file.path(D, "fib_coexpr_stored_corr.csv"), show_col_types = FALSE) |>
  filter(gene1 < gene2, gene1 %in% genes$name, gene2 %in% genes$name, r >= EDGE_R) |>
  rename(from = gene1, to = gene2)
g <- graph_from_data_frame(edges, vertices = genes, directed = FALSE)
cat(sprintf("%d nodes, %d edges (r >= %.2f); isolated: %s\n", vcount(g), ecount(g), EDGE_R,
            paste(V(g)$name[degree(g) == 0], collapse = ", ")))
set.seed(7)
# stress majorisation on distance 1 - r: keeps the dense ECM block spread enough to label, and
# graphlayouts packs disconnected components next to each other instead of flinging them apart
lay <- create_layout(g, layout = "stress", weights = 1 - E(g)$r)
# module name placed above each hull
mod_lab <- lay |> group_by(module) |> summarise(x = mean(x), y = max(y) + 0.17 * diff(range(lay$y)))
p <- ggraph(lay) +
  geom_mark_hull(aes(x, y, group = module, fill = module), colour = NA, alpha = 0.12,
                 expand = unit(4, "mm"), radius = unit(3, "mm"), concavity = 4) +
  geom_edge_link(aes(width = r, alpha = r), colour = "grey45") +
  scale_edge_width(range = c(0.15, 1.1), breaks = c(0.3, 0.5, 0.7), name = "co-expression r") +
  scale_edge_alpha(range = c(0.25, 0.8), guide = "none") +
  geom_node_point(aes(fill = module), shape = 21, size = 3.2, colour = "white", stroke = 0.4) +
  geom_node_text(aes(label = name, colour = cav_r > 0), size = 2.1, family = FIG_FONT, fontface = "bold",
                 vjust = -0.9, show.legend = FALSE) +
  geom_text(data = mod_lab, aes(x, y, label = module), inherit.aes = FALSE, size = 2.6,
            family = FIG_FONT, colour = "grey20", fontface = "italic") +
  scale_fill_manual(values = MOD_COLOR, name = "module") +
  scale_colour_manual(values = c(`TRUE` = SC_CANCER_COLOR, `FALSE` = SC_NORMAL_COLOR),
                      guide = "none") +
  coord_equal(clip = "off") +
  labs(title = "Fibroblast — colorectum: co-expression of Fig 6 fibroblast programs + T / B clusters (top/bottom 20 by CAV r)",
       subtitle = sprintf("Edges: Pearson r ≥ %.2f across 2000 paired-donor fibroblasts; stress layout; gene label red = up in cancer (CAV r > 0), blue = down", EDGE_R)) +
  theme_void(base_family = FIG_FONT, base_size = 7) +
  theme(plot.title = element_text(size = 7.5), plot.subtitle = element_text(size = 6, colour = "grey30"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6.5),
        plot.margin = margin(6, 6, 6, 6))
ggsave("previews/fib_coexpr_network_stored.png", p, width = 8.5, height = 4, dpi = 250, bg = "white")
ggsave("previews/fib_coexpr_network_stored.pdf", p, width = 8.5, height = 4, device = cairo_pdf)
