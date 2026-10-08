#!/usr/bin/env Rscript
# draft_fig6_C1e-4.R -- PREVIEW ONLY (not wired into figures.R).
# Renders Figure 6 with the live figures.R code (header + Fig 6 block, extracted
# verbatim at run time) but reading sc_*.csv from a staging dir, so the current
# and the library-consistent (shared scaler, no PCA, C=1e-4, converged) CAVs can be
# compared side by side.
#
# Reads: figure_data_staging/singlecell_C1e-4/<arm>/sc_*.csv
#        (single_cell/convergence_check/07_export_staging.sh; <arm>=stored reproduces
#        figure_data/ to rounding)
# Writes: previews/fig6_C1e-4_<arm>.{png,pdf}
# Preview edits applied to the extracted Fig 6 code (figures.R itself untouched):
#   - scatter labels only the continuum (top/bottom 5) genes; DE genes stay blue, unlabelled
#   - no Spearman subtitle; the stats go to previews/fig6_C1e-4_<arm>_stats.tsv (and METHODS §6.4)
#   - "ecm" as 2nd arg: TEMPORARY highlight of the quiescent-fibroblast ECM program in the scatters
# Usage (from figures/): Rscript draft_fig6_C1e-4.R [stored|np_conv4] [ecm]

a   <- commandArgs(trailingOnly = TRUE)
arm <- if (length(a) >= 1) a[1] else "np_conv4"
ECM <- length(a) >= 2 && a[2] == "ecm"
SC_ECM_GENES <- c("MGP", "DCN", "OGN", "C3", "CCDC80", "MFAP5", "SFRP2", "DPT", "LUM")
SC_ECM_COLOR <- "#009E73"
sc_stats <- list()
sub_once <- function(pat, rep, x) {        # fail loudly if figures.R drifts
  hit <- grepl(pat, x, fixed = TRUE); stopifnot(sum(hit) == 1)
  x[hit] <- sub(pat, rep, x[hit], fixed = TRUE); x
}
src <- readLines("figures.R")
hdr_end   <- grep("^load_ont <- function", src)[1] - 1
blk_start <- grep("^sc_meta_path <- ", src)[1]
blk_end   <- grep('message\\("Saved cav_fig6.pdf', src)[1]
code <- c(src[1:hdr_end], src[blk_start:blk_end], "}")
code <- sub('^DATA <- "figure_data"', sprintf('DATA <- "figure_data_staging/singlecell_C1e-4/%s"', arm), code)
code <- sub('^OUT  <- "figures"', 'OUT  <- "previews"', code)
code <- gsub('"cav_fig6\\.(pdf|png)"', sprintf('"fig6_C1e-4_%s.\\1"', arm), code)
tag  <- paste0(arm, if (ECM) "_ecm" else "")
code <- gsub(sprintf("fig6_C1e-4_%s.", arm), sprintf("fig6_C1e-4_%s.", tag), code, fixed = TRUE)
# 1. labels: continuum genes only
code <- sub_once("label_ids  <- union(auto_label, df_stat %>% filter(category != \"background\") %>% pull(gene))",
                 "label_ids  <- c(df_stat %>% filter(category == \"continuum\") %>% pull(gene), if (ECM) df$gene[df$category == \"ecm\"])", code)
# 2. subtitle -> stats table
code <- sub_once('subtitle = sprintf("Spearman r = %.2f (n = %s; DE %d donors)", rho,',
                 'subtitle = NULL, caption = { sc_stats[[pair_id]] <<- tibble(pair = pair_id, spearman_rho = rho, n_genes = n_rho, de_donors = df$n_donors[1], p = suppressWarnings(cor.test(df_all$log2fc, df_all$cav_r, method = "spearman", exact = FALSE)$p.value)); NULL }, dummy = sprintf("%.2f %s %d", rho,', code)
# 3. ECM highlight category, layer and label colour
code <- sub_once('TRUE                      ~ "below_threshold"',
                 'ECM & gene_name %in% SC_ECM_GENES ~ "ecm", TRUE ~ "below_threshold"', code)
code <- sub_once("      geom_text_repel(",
                 "      geom_point(data = filter(df, category == \"ecm\"), shape = 21, fill = SC_ECM_COLOR, color = \"#00513c\", size = 1.6, stroke = 0.3) +\n      geom_text_repel(", code)
code <- sub_once("continuum  = SC_CONTINUUM_EDGE)) +", "continuum  = SC_CONTINUUM_EDGE, ecm = \"#00513c\")) +", code)
stopifnot(any(grepl(arm, code, fixed = TRUE)))
eval(parse(text = code))
write_tsv(bind_rows(sc_stats), file.path("previews", sprintf("fig6_C1e-4_%s_stats.tsv", arm)))
