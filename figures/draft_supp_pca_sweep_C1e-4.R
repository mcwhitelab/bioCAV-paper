#!/usr/bin/env Rscript
# draft_supp_pca_sweep_C1e-4.R -- PREVIEW ONLY (not wired into figures.R; figure_data/ untouched).
# C = 1e-4 converged version of fig_supp_pca_sweep, judged by CACGTG (canonical E-box) binding
# instead of bHLH family counts, plus a preview of the Fig 5B MYC top-20 (PCA-128).
# Reads:  figure_data_staging/pca_sweep_C1e-4/{pca_sweep_cacgtg.csv, fig5B_myc_top20.csv}
#         (interface_search/fastas/pairs/myc_mad/pca_sweep_supp/export_pca_sweep_cacgtg_staging.py)
# Writes: previews/supp_pca_sweep_C1e-4.{png,pdf}, previews/fig5B_myc_top20_C1e-4.{png,pdf}
# Uses the live figures.R header (packages, palette, base_theme), extracted at run time.
src <- readLines("figures.R")
eval(parse(text = src[1:(grep("^load_ont <- function", src)[1] - 1)]))
STG <- "figure_data_staging/pca_sweep_C1e-4"
LABEL_SET <- "db+lit+inf"                       # database motifs + literature fact-check
d <- read_csv(file.path(STG, "pca_sweep_cacgtg.csv"), show_col_types = FALSE)
CON_COL <- c(MYC = unname(oi["blue"]), MAD1 = unname(oi["orange"]))
KS <- sort(unique(d$k[d$arm == "C1e-4_converged"]))
KB <- c(8, 16, 32, 64, 128, 256, 512, 1152)
k_lab <- function(k) ifelse(k == 1152, "1152\n(no PCA)", as.character(k))
dd <- d |> filter(arm == "C1e-4_converged") |> mutate(concept = factor(concept, c("MYC", "MAD1")))
main <- filter(dd, label_set == LABEL_SET)
dbonly <- filter(dd, label_set == "db")
band <- annotate("rect", xmin = 96, xmax = 160, ymin = -Inf, ymax = Inf, fill = "grey92")
vline <- geom_vline(xintercept = 128, linetype = "dotted", linewidth = 0.3, color = "grey45")
xscale <- scale_x_continuous(trans = "log2", breaks = KB, labels = k_lab(KB))
panel <- function(y, ylab, ylim, extra = NULL) {
  ggplot(main, aes(k, .data[[y]], color = concept)) + band + vline + extra +
    geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
    scale_color_manual(values = CON_COL, labels = c(MYC = "MYC 366–420", MAD1 = "MAD1 61–94"), name = NULL) +
    xscale + scale_y_continuous(limits = ylim) +
    labs(x = "PCA components", y = ylab) + base_theme() +
    theme(legend.position = "none", axis.text.x = element_text(size = 7), plot.margin = margin(4, 12, 2, 4))
}
pA <- panel("auroc_within_bhlh", "CACGTG binders vs other bHLH\n(AUROC)", c(0.4, 0.92),
            list(geom_hline(yintercept = 0.5, linewidth = 0.3, color = "grey60"),
                 geom_line(data = dbonly, linewidth = 0.35, linetype = "22")))
pB <- panel("binders_top100", "CACGTG binders\nin top 100", c(0, max(main$binders_top100) + 2))
pC <- panel("bhlh_top100", "bHLH proteins\nin top 100", c(0, 102))
leg <- get_legend(pA + theme(legend.position = "top", legend.text = element_text(size = 8, margin = margin(l = 3, r = 10)), legend.key.width = unit(14, "pt")))
note <- ggdraw() + draw_label("shaded: k = 96–160;  dashed (A): database motifs only", size = 7,
                              color = "grey40", fontfamily = FIG_FONT)
fig <- plot_grid(plot_grid(leg, note, nrow = 1, rel_widths = c(1, 1.2)),
                 plot_grid(pA, pB, pC, nrow = 1, align = "h", labels = LETTERS[1:3],
                           label_size = 9, label_fontfamily = FIG_FONT),
                 ncol = 1, rel_heights = c(0.1, 1))
dir.create("previews", showWarnings = FALSE)
ggsave("previews/supp_pca_sweep_C1e-4.pdf", fig, width = 7, height = 2.6, bg = "white", device = cairo_pdf)
ggsave("previews/supp_pca_sweep_C1e-4.png", fig, width = 7, height = 2.6, bg = "white", dpi = 300)

# --- Fig 5B preview: MYC top 20, C = 1e-4, PCA-128 ---
t20 <- read_csv(file.path(STG, "fig5B_myc_top20.csv"), show_col_types = FALSE) |>
  mutate(label = paste0(gene, c(literature = " \u2020", inferred = " \u2021", database = "", none = "")[call_source]),
         label = factor(label, rev(label)),
         cacgtg_call = factor(cacgtg_call, c("CACGTG", "other motif", "no motif")))
CALL_COL <- c("CACGTG" = unname(oi["blue"]), "other motif" = unname(oi["vermillion"]), "no motif" = "grey60")
p5b <- ggplot(t20, aes(cav_score, label, color = cacgtg_call)) +
  geom_point(size = 2) +
  scale_color_manual(values = CALL_COL, name = NULL, drop = FALSE) +
  labs(x = "CAV score", y = NULL,
       caption = "\u2020 literature, \u2021 inferred; others HOCOMOCO/JASPAR") +
  base_theme() +
  theme(legend.position = "top", legend.text = element_text(size = 7, margin = margin(l = 1, r = 6)),
        axis.text.y = element_text(size = 7), plot.caption = element_text(size = 6, color = "grey40", hjust = 0),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))
ggsave("previews/fig5B_myc_top20_C1e-4.pdf", p5b, width = 2.6, height = 3.6, bg = "white", device = cairo_pdf)
ggsave("previews/fig5B_myc_top20_C1e-4.png", p5b, width = 2.6, height = 3.6, bg = "white", dpi = 300)
message("Saved previews/supp_pca_sweep_C1e-4.* and previews/fig5B_myc_top20_C1e-4.*")
