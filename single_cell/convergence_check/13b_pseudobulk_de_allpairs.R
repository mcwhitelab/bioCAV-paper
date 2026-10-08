# COPY of scripts/paired_donor_pseudobulk_de.R with all 10 pairs + staging output (2026-10-05)
#!/usr/bin/env Rscript
# paired_donor_pseudobulk_de.R -- conventional pseudobulk DE for Figure 6's
# three case-study pairs (second, true-log2FC version of the DE axis).
#
# Protocol (pseudobulk, paired design; Squair et al. 2021 Nat Commun,
# Crowell et al. 2020 Nat Commun / muscat):
#   * input: raw counts summed per donor x condition
#     (scripts/paired_donor_pseudobulk.py), paired donors only
#   * drop pseudobulk samples with < MIN_CELLS cells; a donor left with only
#     one condition is then dropped too (paired design needs both)
#   * gene filter: edgeR::filterByExpr(group = condition), same universe for
#     both engines
#   * primary:   DESeq2, design ~ donor + condition (normal = reference),
#                Wald test, log2FC shrunk with lfcShrink(type = "apeglm")
#   * robustness: edgeR quasi-likelihood, TMM, same design,
#                glmQLFit(robust = TRUE) + glmQLFTest
#
# Output: results/de_pseudobulk_paired_donors/min<MIN_CELLS>/<pair>.tsv
#
# Usage (from single_cell/):
#   Rscript scripts/paired_donor_pseudobulk_de.R [MIN_CELLS=10]

suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(edgeR); library(data.table)
})

args      <- commandArgs(trailingOnly = TRUE)
MIN_CELLS <- if (length(args)) as.integer(args[1]) else 10L
LIB       <- "cav_library/b617ee1b-f8c8-4de9-b82b-e803ab93550d"
IN_DIR    <- "convergence_check/de_pseudobulk_allpairs/pseudobulk_counts"
OUT_DIR   <- file.path("convergence_check/de_pseudobulk_allpairs", sprintf("min%d", MIN_CELLS))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PAIRS <- read.csv("static/paired_donor_pairs.csv")$pair_label

genes   <- fread(file.path(IN_DIR, "genes.tsv"))
summary <- list()

for (pair in PAIRS) {
  cts  <- fread(cmd = paste("zcat", shQuote(file.path(IN_DIR, paste0(pair, ".counts.tsv.gz")))))
  gid  <- cts[[1]]
  cts  <- as.matrix(cts[, -1]); rownames(cts) <- gid
  smp  <- as.data.frame(fread(file.path(IN_DIR, paste0(pair, ".samples.tsv"))))
  rownames(smp) <- smp$sample
  stopifnot(identical(colnames(cts), smp$sample))

  # sample QC: min cells, then keep only donors still paired
  smp  <- smp[smp$n_cells >= MIN_CELLS, ]
  ok   <- names(which(table(smp$donor) == 2))
  dropped <- setdiff(unique(fread(file.path(IN_DIR, paste0(pair, ".samples.tsv")))$donor), ok)
  smp  <- smp[smp$donor %in% ok, ]
  cts  <- cts[, smp$sample]
  smp$donor     <- factor(smp$donor)
  smp$condition <- factor(smp$condition, levels = c("normal", "cancer"))
  n_don <- nlevels(smp$donor)
  cat(sprintf("%s: %d donors (dropped: %s), %d samples\n", pair, n_don,
              if (length(dropped)) paste(dropped, collapse = ",") else "none", ncol(cts)))
  if (n_don < 2) { cat("  <2 paired donors -- skipped\n"); next }

  keep <- filterByExpr(cts, group = smp$condition)
  cts  <- cts[keep, ]
  cat(sprintf("  %d genes pass filterByExpr\n", nrow(cts)))

  # --- DESeq2 (primary) ---
  dds <- DESeqDataSetFromMatrix(cts, smp, design = ~ donor + condition)
  dds <- DESeq(dds, quiet = TRUE)
  coef_name <- "condition_cancer_vs_normal"
  res  <- results(dds, name = coef_name, alpha = 0.05)
  shr  <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  out  <- data.table(gene          = rownames(res),
                     baseMean      = res$baseMean,
                     log2fc_mle    = res$log2FoldChange,
                     lfcSE         = res$lfcSE,
                     log2fc        = shr$log2FoldChange,   # apeglm-shrunk; plotted
                     pval          = res$pvalue,
                     padj          = res$padj)

  # --- edgeR QL (robustness) ---
  y  <- DGEList(cts, group = smp$condition)
  y  <- calcNormFactors(y)
  mm <- model.matrix(~ donor + condition, data = smp)
  y  <- estimateDisp(y, mm)
  fit <- glmQLFit(y, mm, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "conditioncancer")
  tt  <- topTags(qlf, n = Inf, sort.by = "none")$table
  out[, `:=`(edgeR_log2fc = tt[gene, "logFC"],
             edgeR_pval   = tt[gene, "PValue"],
             edgeR_padj   = tt[gene, "FDR"])]

  out <- merge(genes, out, by = "gene", sort = FALSE)
  setorder(out, -baseMean)
  out[, `:=`(n_donors = n_don, min_cells = MIN_CELLS)]
  fwrite(out, file.path(OUT_DIR, paste0(pair, ".tsv")), sep = "\t")

  summary[[pair]] <- data.table(
    pair = pair, n_donors = n_don, dropped_donors = paste(dropped, collapse = ","),
    n_samples = ncol(cts), resid_df = ncol(cts) - n_don - 1L, n_genes = nrow(out),
    deseq2_sig    = sum(out$padj < 0.05 & abs(out$log2fc) >= 1, na.rm = TRUE),
    edgeR_sig     = sum(out$edgeR_padj < 0.05 & abs(out$edgeR_log2fc) >= 1, na.rm = TRUE),
    spearman_deseq2_edgeR = cor(out$log2fc_mle, out$edgeR_log2fc, method = "spearman"))
}

summary <- rbindlist(summary)
print(summary)
fwrite(summary, file.path(OUT_DIR, "summary.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
