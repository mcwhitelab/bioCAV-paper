#!/usr/bin/env Rscript
# DRAFT / exploratory -- (A) 3-D CAV direction-space UMAP colored by Pfam
# clan (7688 concepts from the 7692-motif library), single azimuth-330deg
# view with a clan legend and visible x/y/z axis lines, no title. Coordinates
# decoded directly from the embedded Plotly data in the precomputed
# motif_clustering/results/figures/pfam_umap_3d.html (not recomputed).
# The 2D projection panel and the other rotation angles were dropped after
# review -- 330 deg alone showed the branching structure most clearly.
# (B) blank panel reserved for a manually-assembled illustration of new
# zinc fingers found by guilt-by-association in the same embedding.
# Not part of the figures.R pipeline; run standalone with:
#   Rscript draft_fig_pfam_clan_umap.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(png)
  library(grid)
})

OUT <- "figures"

# draw_image() needs the `magick` package, which isn't installed here --
# use a plain grid rasterGrob instead (png + grid only).
umap3d_png <- readPNG(file.path(OUT, "pfam_clan_umap_3d_mosaic.png"))
p_umap3d <- ggdraw() + draw_grob(rasterGrob(umap3d_png, interpolate = TRUE))

p_blank <- ggdraw() +
  draw_label("(reserved for manual guilt-by-association\nzinc finger illustration)",
             size = 8, color = "grey70", fontface = "italic")

fig <- plot_grid(p_umap3d, p_blank, nrow = 1, labels = "AUTO", label_size = 10,
                 rel_widths = c(1.3, 0.7))

ggsave(file.path(OUT, "draft_fig_pfam_clan_umap.pdf"), fig,
       width = 15, height = 7, bg = "white", device = cairo_pdf)
ggsave(file.path(OUT, "draft_fig_pfam_clan_umap.png"), fig,
       width = 15, height = 7, bg = "white", dpi = 200)
message("Saved draft_fig_pfam_clan_umap.pdf/.png")
