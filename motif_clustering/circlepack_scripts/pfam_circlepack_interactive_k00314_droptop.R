library(tidyverse)
library(jsonlite)

# Builds a self-contained interactive HTML from a cached circlepack layout:
# free-text search with highlighting, zoom-to-node, and hover tooltips instead
# of printed labels. Run pfam_circlepack.R first to produce the layout cache.

infile <- "/groups/clairemcwhite/cav_workspace/bioCAV-paper/motif_clustering/results/figures/circlepack_input_diffusion.csv"
outdir <- "/groups/clairemcwhite/cav_workspace/bioCAV-paper/motif_clustering/results/figures"

# Reads the droptop layout, in which the largest top-cut clusters were removed
# from the hierarchy and their next-cut children promoted to the top level --
# see pfam_circlepack_k00314_droptop.R, which must be run first. Nothing is
# hidden at draw time here: every node in the layout is drawn.
layout_tag <- "_k00314_droptop3"

cut_cols <- c("cut_k00314", "cut_k02769", "cut_k04629", "cut_k09172", "cut_k18984")
cut_names <- sprintf("k = %s", format(as.integer(sub("cut_k0*", "", cut_cols)), big.mark = ",", trim = TRUE))

pfam_all <- read_csv(infile, show_col_types = FALSE)

# Same colour recipe as the static figures, so the two agree
# Clan colours live in a table so a clan keeps its colour across figures and
# across Pfam releases: existing entries are reused and only clans new to the
# data get a fresh hue appended.  812 clans is far past any palette that stays
# distinguishable, so hues are sampled rather than assigned: the point is seeing
# which branches share a clan, not looking any one clan up.
clan_colour_table <- file.path(outdir, "clan_colours_R.csv")
clan_ids <- sort(unique(na.omit(pfam_all$clan_acc)))
colour_tab <- if (file.exists(clan_colour_table)) {
    read_csv(clan_colour_table, show_col_types = FALSE)
} else {
    tibble(clan_acc = character(), colour = character())
}
new_ids <- setdiff(clan_ids, colour_tab$clan_acc)
if (length(new_ids) > 0) {
    set.seed(42 + nrow(colour_tab))
    colour_tab <- bind_rows(colour_tab, tibble(
        clan_acc = new_ids,
        colour   = hcl(h = sample(seq(0, 360, length.out = length(new_ids))),
                       c = runif(length(new_ids), 70, 120),
                       l = runif(length(new_ids), 35, 65))))
    write_csv(colour_tab, clan_colour_table)
    cat(sprintf("added %d clan colour(s) new to this release -> %s\n",
                length(new_ids), basename(clan_colour_table)))
}
clan_colours <- setNames(colour_tab$colour[match(clan_ids, colour_tab$clan_acc)], clan_ids)
stopifnot(!any(is.na(clan_colours)))

# The JS splits label fields on "|", so make sure the data has none
stopifnot(!any(grepl("|", with(pfam_all,
    paste(accession, pfam_short_name, pfam_long_name, clan_acc, clan_name)), fixed = TRUE)))

template <- readLines("/groups/clairemcwhite/cav_workspace/bioCAV-paper/motif_clustering/circlepack_scripts/interactive_template.html", warn = FALSE) %>%
    paste(collapse = "\n")

lay <- readRDS(file.path(outdir, paste0("circlepack_layout", layout_tag, ".rds")))

{
    variant <- layout_tag

    # Promotion breaks the depth == cut level correspondence, so split tips from
    # clusters by name, not by depth, and read each cluster's level off its id.
    stopifnot(!any(grepl("^[0-9]", pfam_all$accession)))
    lay <- lay %>% mutate(is_tip  = name %in% pfam_all$accession,
                          cut_idx = match(sub("^[0-9]+", "", name), cut_cols))

    clusters <- lay %>% filter(depth >= 1, !is_tip)
    tips <- lay %>% filter(is_tip) %>%
        select(name, x, y, r) %>%
        left_join(pfam_all, by = c("name" = "accession"))
    stopifnot(nrow(tips) == nrow(pfam_all), !any(is.na(clusters$cut_idx)))

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
                  r = round(clusters$r, 3), d = as.integer(clusters$cut_idx)),
        tp = list(x = round(tips$x, 2), y = round(tips$y, 2),
                  r = round(tips$r, 3), c = as.integer(tips$colour_idx)),
        lab = tips$label,
        cid = as.integer(cids)
    )

    subtitle <- paste0(format(nrow(tips), big.mark = ","), " Pfam families &middot; ",
                       length(unique(tips$clan_acc[tips$clan_acc != ""])), " clans",
                       " &middot; 3 largest ", cut_names[1], " clusters removed from ",
                       "the hierarchy, their ", cut_names[2],
                       " sub-clusters promoted to the top level")

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
