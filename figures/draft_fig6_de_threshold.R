#!/usr/bin/env Rscript
# draft_fig6_de_threshold.R -- PREVIEW ONLY. Live Fig 6 code with a stricter DE call for the
# blue scatter points (default padj < 0.01 and |log2FC| >= 2; live is padj < 0.05, |log2FC| >= 1).
# Usage (from figures/): Rscript draft_fig6_de_threshold.R [PADJ=0.01] [LFC=2]
# Writes: previews/fig6_de_padj<PADJ>_lfc<LFC>.{png,pdf}
a <- commandArgs(trailingOnly = TRUE)
PADJ <- if (length(a) >= 1) a[1] else "0.01"; LFC <- if (length(a) >= 2) a[2] else "2"
sub_once <- function(pat, rep, x) { hit <- grepl(pat, x, fixed = TRUE); stopifnot(sum(hit) == 1)
  x[hit] <- sub(pat, rep, x[hit], fixed = TRUE); x }
src <- readLines("figures.R")
hdr_end   <- grep("^load_ont <- function", src)[1] - 1
blk_start <- grep("^sc_meta_path <- ", src)[1]
blk_end   <- grep('message\\("Saved cav_fig6.pdf', src)[1]
code <- c(src[1:hdr_end], src[blk_start:blk_end], "}")
code <- sub('^OUT  <- "figures"', 'OUT  <- "previews"', code)
code <- gsub('"cav_fig6\\.(pdf|png)"', sprintf('"fig6_de_padj%s_lfc%s.\\1"', PADJ, LFC), code)
code <- sub_once("padj < 0.05 & abs(log2fc) >= 1", sprintf("padj < %s & abs(log2fc) >= %s", PADJ, LFC), code)
code <- gsub("padj < 0.05, |log2FC| \\u2265 1", sprintf("padj < %s, |log2FC| \\u2265 %s", PADJ, LFC), code, fixed = TRUE)
eval(parse(text = code))
