library(tidyverse)
library(jsonlite)

# Builds a self-contained interactive HTML from a cached circlepack layout:
# free-text search with highlighting, zoom-to-node, and hover tooltips instead
# of printed labels. Run pfam_circlepack.R first to produce the layout cache.

infile <- "/groups/clairemcwhite/cav_workspace/bioCAV-paper/motif_clustering/results/figures/circlepack_input_diffusion.csv"
outdir <- "/groups/clairemcwhite/cav_workspace/bioCAV-paper/motif_clustering/results/figures"

# One cached layout over the full data, built into several HTMLs that differ only
# in how many of the largest top-cut clusters are drawn without outlines. 0 draws
# everything. Hiding keeps the families and their sub-clusters, and suppresses the
# outline of those clusters plus their children at the next cut -- same rule as
# pfam_circlepack.R.
layout_tag <- "_k00314"
hide_counts <- c(0, 3)

cut_cols <- c("cut_k00314", "cut_k02769", "cut_k04629", "cut_k09172", "cut_k18984")
cut_names <- sprintf("k = %s", formatC(as.integer(sub("cut_k0*", "", cut_cols)), big.mark = ",", format = "d"))

pfam_all <- read_csv(infile, show_col_types = FALSE)

# Same colour recipe as the static figures, so the two agree
set.seed(42)
clan_ids <- sort(unique(na.omit(pfam_all$clan_acc)))
clan_colours <- setNames(
    hcl(h = sample(seq(0, 360, length.out = length(clan_ids))),
        c = runif(length(clan_ids), 70, 120),
        l = runif(length(clan_ids), 35, 65)),
    clan_ids)

# The JS splits label fields on "|", so make sure the data has none
stopifnot(!any(grepl("|", with(pfam_all,
    paste(accession, pfam_short_name, pfam_long_name, clan_acc, clan_name)), fixed = TRUE)))

template <- readLines("/groups/clairemcwhite/cav_workspace/bioCAV-paper/motif_clustering/circlepack_scripts/interactive_template.html", warn = FALSE) %>%
    paste(collapse = "\n")

lay <- readRDS(file.path(outdir, paste0("circlepack_layout", layout_tag, ".rds")))

for (hide_top_n in hide_counts) {
    variant <- paste0(layout_tag, if (hide_top_n > 0) paste0("_hidetop", hide_top_n) else "")

    # Cluster nodes left undrawn: the n largest at the top cut and their children
    # at the next cut. Names are clusternum pasted to the column name, as built
    # by pfam_circlepack.R.
    hidden_nodes <- character(0)
    if (hide_top_n > 0) {
        biggest <- pfam_all %>% count(.data[[cut_cols[1]]], sort = TRUE) %>%
            slice(seq_len(hide_top_n))
        inside <- pfam_all %>% filter(.data[[cut_cols[1]]] %in% biggest[[1]])
        hidden_nodes <- c(paste0(biggest[[1]], cut_cols[1]),
                          paste0(sort(unique(inside[[cut_cols[2]]])), cut_cols[2]))
    }

    clusters <- lay %>% filter(depth >= 1, depth <= length(cut_cols),
                               !name %in% hidden_nodes)
    tips <- lay %>% filter(depth == length(cut_cols) + 1) %>%
        select(name, x, y, r) %>%
        left_join(pfam_all, by = c("name" = "accession"))

    tips <- tips %>%
        mutate(clan_acc = replace_na(clan_acc, ""),
               clan_name = replace_na(clan_name, ""),
               colour_idx = match(clan_acc, clan_ids) - 1L,          # 0-based, NA when no clan
               colour_idx = replace_na(colour_idx, -1L),
               label = paste(name, pfam_short_name, pfam_long_name,
                             clan_acc, clan_name, sep = "|"))

    # cluster ids at each cut, flattened row-major for the tooltip
    cids <- tips %>% select(all_of(cut_cols)) %>% as.matrix() %>% t() %>% as.vector()

    payload <- list(
        meta = list(cutNames = cut_names, nTips = nrow(tips), nClusters = nrow(clusters)),
        palette = unname(clan_colours),
        cl = list(x = round(clusters$x, 2), y = round(clusters$y, 2),
                  r = round(clusters$r, 3), d = as.integer(clusters$depth)),
        tp = list(x = round(tips$x, 2), y = round(tips$y, 2),
                  r = round(tips$r, 3), c = as.integer(tips$colour_idx)),
        lab = tips$label,
        cid = as.integer(cids)
    )

    subtitle <- paste0(format(nrow(tips), big.mark = ","), " Pfam families &middot; ",
                       length(unique(tips$clan_acc[tips$clan_acc != ""])), " clans",
                       if (hide_top_n > 0)
                           paste0(" &middot; ", hide_top_n, " largest ", cut_names[1],
                                  " clusters drawn without outlines") else "")

    # split rather than sub(): the payload is megabytes and must not be read as a
    # regex replacement, where backslashes would be reinterpreted
    pieces <- strsplit(template, "/*__DATA__*/null", fixed = TRUE)[[1]]
    stopifnot(length(pieces) == 2)
    html <- paste0(pieces[1],
                   toJSON(payload, auto_unbox = TRUE, digits = NA),
                   pieces[2]) %>%
        gsub("__TITLE__", "Pfam circle pack", ., fixed = TRUE) %>%
        gsub("__SUBTITLE__", subtitle, ., fixed = TRUE)

    out <- file.path(outdir, paste0("pfam_circlepack_interactive", variant, ".html"))
    writeLines(html, out, useBytes = TRUE)
    cat(sprintf("%s: %d tips, %d clusters, %.1f MB\n",
                basename(out), nrow(tips), nrow(clusters), file.size(out) / 1024^2))
}
