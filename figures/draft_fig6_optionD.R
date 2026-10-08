#!/usr/bin/env Rscript
# draft_fig6_optionD.R -- PREVIEW ONLY (not wired into figures.R).
# Full Fig 6 mock for "option 2" (2026-10-05): keep fibroblast/colorectum as the deliberate
# contrast case (no disease separation, no DE; METHODS §6.4 items 22b/22c) and replace the lower
# half of D (continuum strips) with a fibroblast-program panel whose x axis is labelled
# "L2 score (cancer-CAV direction)" rather than "<- normal  cancer ->".
# Live figure_data/ and the live figures.R Fig 6 code (extracted verbatim), plus the preview
# edits already requested: no Spearman subtitle (stats in METHODS §6.4 item 21), scatter labels
# only the continuum (top/bottom 5) genes.
# Usage (from figures/): Rscript draft_fig6_optionD.R [lines|blocks]
#   lines  = module score vs L2 (draft_fig6_fib_modules.R)
#   blocks = shrunken co-expression heatmap of the module genes (draft_fig6_fib_module_summary.R, C)
# Writes: previews/fig6_optionD_<variant>.{png,pdf}
a <- commandArgs(trailingOnly = TRUE); variant <- if (length(a)) a[1] else "lines"
FIB_PAIR <- "fibroblast__colorectum__normal_vs_colorectal_cancer"
penv <- new.env()   # build the D panel in its own environment so its theme helpers don't leak
if (variant == "lines") { sys.source("draft_fig6_fib_modules.R", envir = penv); D_BOTTOM <- penv$p
} else { sys.source("draft_fig6_fib_module_summary.R", envir = penv); D_BOTTOM <- penv$pC }
sub_once <- function(pat, rep, x) { hit <- grepl(pat, x, fixed = TRUE); stopifnot(sum(hit) == 1)
  x[hit] <- sub(pat, rep, x[hit], fixed = TRUE); x }
src <- readLines("figures.R")
hdr_end   <- grep("^load_ont <- function", src)[1] - 1
blk_start <- grep("^sc_meta_path <- ", src)[1]
blk_end   <- grep('message\\("Saved cav_fig6.pdf', src)[1]
code <- c(src[1:hdr_end], src[blk_start:blk_end], "}")
code <- sub('^OUT  <- "figures"', 'OUT  <- "previews"', code)
code <- gsub('"cav_fig6\\.(pdf|png)"', sprintf('"fig6_optionD_%s.\\1"', variant), code)
code <- sub_once('label_ids  <- union(auto_label, df_stat %>% filter(category != "background") %>% pull(gene))',
                 'label_ids  <- df_stat %>% filter(category == "continuum") %>% pull(gene)', code)
code <- sub_once('subtitle = sprintf("Spearman r = %.2f (n = %s; DE %d donors)", rho,',
                 'subtitle = NULL, dummy = sprintf("%.2f %s %d", rho,', code)
# D: lines are aligned with the scatter above; the square block heatmap is not (alignment would
# squeeze the scatter to the heatmap's width)
D_ALIGN <- if (variant == "lines") "v" else "none"
code <- sub_once("plot_grid(sc_build_scatter(pid), sc_build_continuum(pid, brackets = (i == 1)),",
                 "if (pid == FIB_PAIR) plot_grid(sc_build_scatter(pid), D_BOTTOM, ncol = 1, align = D_ALIGN, axis = \"lr\", rel_heights = c(0.62, 1)) else plot_grid(sc_build_scatter(pid), sc_build_continuum(pid, brackets = (i == 1)),", code)
eval(parse(text = code))
