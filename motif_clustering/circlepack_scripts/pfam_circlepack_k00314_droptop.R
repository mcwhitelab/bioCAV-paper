library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)

# Circle-packed view of the Pfam dendrogram, cut at several heights.
# Adapted from https://gist.github.com/clairemcwhite/12965e7d221281b32ac7832dbbc8ff80
# The cut_* columns are already present in the input, so the cutree step is skipped.

infile <- "/groups/clairemcwhite/cav_workspace/bioCAV-paper/motif_clustering/results/figures/circlepack_input_diffusion.csv"
outdir <- "/groups/clairemcwhite/cav_workspace/bioCAV-paper/motif_clustering/results/figures"

# ---- the two things to change ----
cuts <- c("cut_k00314", "cut_k02769", "cut_k04629", "cut_k09172", "cut_k18984")
# The top cut has a few clusters far larger than the rest, and their outlines
# swallow the figure. How many of the largest to DROP FROM THE HIERARCHY (0
# keeps all). Unlike the hidetop variant, which only suppressed outlines and
# left the packing untouched, this removes the top-cut node itself and promotes
# its next-cut children to the top level, where they are drawn like any other
# cluster. Every family is kept; the sub-clusters stay. Because the hierarchy
# changes, the whole figure repacks -- circles move everywhere, not only inside
# the dropped clusters.
drop_top_n <- 3
# ----------------------------------

drop_largest_top <- drop_top_n > 0
# Dropping a level changes the hierarchy, so this variant needs its own pack --
# it must NOT share the hidetop layout cache.
layout_tag <- paste0("_", sub("cut_", "", cuts[1]),
                     if (drop_largest_top) paste0("_droptop", drop_top_n) else "")
tag <- layout_tag

pfam_all <- read_csv(infile, show_col_types = FALSE) %>%
    rename(ID = accession) %>%
    select(ID, clan_acc, clan_name, pfam_short_name, all_of(cuts))

pfam <- pfam_all   # no families are dropped any more; the two stay in step

# Top-cut clusters to drop from the hierarchy: the n largest. Their next-cut
# children are promoted to the top level and drawn there.
dropped <- character(0)
promoted <- character(0)
if (drop_largest_top) {
    biggest <- pfam %>% count(.data[[cuts[1]]], sort = TRUE) %>% slice(seq_len(drop_top_n))
    inside <- pfam %>% filter(.data[[cuts[1]]] %in% biggest[[1]])
    promoted_kids <- sort(unique(inside[[cuts[2]]]))
    dropped <- paste0(biggest[[1]], cuts[1])
    promoted <- paste0(promoted_kids, cuts[2])
    cat(sprintf("dropping %d largest %s cluster(s) [%s] from the hierarchy; their %d %s children are promoted to the top level; %d of %d families (%.1f%%) sit inside and are still drawn\n",
                drop_top_n, cuts[1], paste(biggest[[1]], collapse = ", "),
                length(promoted_kids), cuts[2],
                sum(biggest$n), nrow(pfam_all), 100 * sum(biggest$n) / nrow(pfam_all)))
}

# Each cluster needs an ID unique across cut levels, not a plain number.
# Borrow the uniquifier from the column header, as in the tutorial.
clusters_uniqued <- pfam %>%
    select(ID, all_of(cuts)) %>%
    pivot_longer(all_of(cuts), names_to = "clusterset", values_to = "clusternum") %>%
    mutate(clusterid = paste0(clusternum, clusterset)) %>%
    select(-clusternum) %>%
    pivot_wider(names_from = clusterset, values_from = clusterid)

# Reparent the dropped clusters' rows straight onto the origin, so the top-cut
# node vanishes from the hierarchy and its children become top level. Done after
# the pivot, on the uniqued ids.
if (drop_largest_top) {
    clusters_uniqued[[cuts[1]]] <- ifelse(clusters_uniqued[[cuts[1]]] %in% dropped,
                                          "origin", clusters_uniqued[[cuts[1]]])
}

# All nodes hang off one origin so the layout has a single root
clusters_uniqued$link <- "origin"

# Edges in from/to order, one level at a time, tips last
levels_ordered <- c("link", cuts, "ID")
edge_list <- map_dfr(seq_len(length(levels_ordered) - 1), function(i) {
    clusters_uniqued %>%
        select(from = all_of(levels_ordered[i]), to = all_of(levels_ordered[i + 1])) %>%
        distinct()
}) %>%
    filter(from != to)          # the origin -> origin self-loop left by reparenting

# The promotion is the whole point of this variant, so check it landed: the top
# level should hold the surviving top-cut clusters plus the promoted children.
n_top <- n_distinct(pfam[[cuts[1]]]) - length(dropped) + length(promoted)
stopifnot(!any(edge_list$from == edge_list$to),
          sum(edge_list$from == "origin") == n_top)
cat("top-level nodes:", n_top, "(", n_distinct(pfam[[cuts[1]]]) - length(dropped),
    "kept top-cut clusters +", length(promoted), "promoted )\n")

mygraph <- as_tbl_graph(edge_list)

# Attach family metadata to the tip nodes
mygraph <- mygraph %>%
    activate("nodes") %>%
    left_join(pfam %>% mutate(name = ID), by = "name") %>%
    mutate(in_clan = case_when(
        is.na(ID)       ~ NA_character_,          # internal cluster node
        is.na(clan_acc) ~ "No clan assigned",
        TRUE            ~ "In a Pfam clan"
    ))

# The pack itself takes a few minutes on 60k+ nodes, so cache it per cut set
layout_cache <- file.path(outdir, paste0("circlepack_layout", layout_tag, ".rds"))
if (file.exists(layout_cache)) {
    lay <- readRDS(layout_cache)
} else {
    lay <- create_layout(mygraph, layout = "circlepack")
    saveRDS(lay, layout_cache)
}
cat("nodes:", nrow(lay), " max depth:", max(lay$depth), "\n")

# depth 0 = origin (not drawn), 1..n = the cuts, n+1 = individual families
n_cut <- length(cuts)
cut_labels <- sprintf("k = %s", format(as.integer(sub("cut_k0*", "", cuts)), big.mark = ",", trim = TRUE))
depth_labels <- c(cut_labels, "Pfam family")

# Sequential blue for the nested cuts, orange tips so single families stay visible
blue_steps <- if (n_cut == 5) {
    c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")
} else {
    colorRampPalette(c("#86b6ef", "#104281"))(n_cut)
}
depth_fills <- setNames(c(blue_steps, "#eb6834"), depth_labels)

# Promotion breaks the depth == cut level correspondence (a promoted cut_k02769
# cluster now sits at depth 1), so read the level off the node name instead --
# cluster ids are the cluster number pasted to the cut column name.
stopifnot(!any(grepl("^[0-9]", pfam$ID)))   # or the regex below eats an accession
w_ramp <- c(0.35, 0.25, 0.18, 0.12, 0.08, 0.04)
lay <- lay %>%
    mutate(is_tip  = name %in% pfam$ID,
           cut_idx = match(sub("^[0-9]+", "", name), cuts),
           lvl_idx = ifelse(is_tip, n_cut + 1L, cut_idx),
           level   = factor(ifelse(depth >= 1, depth_labels[lvl_idx], NA),
                            levels = depth_labels),
           edge_w  = w_ramp[pmin(pmax(lvl_idx, 1), length(w_ramp))])

drawn    <- lay %>% filter(depth >= 1)      # hides the single origin circle only
clusters <- drawn %>% filter(!is_tip)
tips     <- drawn %>% filter(is_tip)
stopifnot(nrow(tips) == nrow(pfam), !any(is.na(clusters$cut_idx)))

# Size the canvas from the packing's own bounding box, so coord_fixed doesn't
# leave gutters down the sides. Different cut sets pack to different aspects.
bb     <- with(drawn, c(min(x - r), max(x + r), min(y - r), max(y + r)))
aspect <- (bb[2] - bb[1]) / (bb[4] - bb[3])
panel_h  <- 12
panel_w  <- panel_h * aspect
fig_h    <- panel_h + 1.5                    # room for title and the 3-line subtitle
fig_w    <- panel_w + 0.5                    # margins only
fig_w_leg <- panel_w + 2.2                   # margins plus the legend column
cat(sprintf("content aspect %.3f -> figure %.1f x %.1f in\n", aspect, fig_w, fig_h))

# Each figure goes out as a 300 dpi PNG for quick viewing and a vector PDF for
# the paper. cairo_pdf embeds the fonts, so the text stays selectable and the
# circles stay vector at any zoom.
devices <- c("png", "pdf")
# Set only_figs (e.g. "clanid") before sourcing to render just those figures --
# this node has one core, so the three-figure run is strictly serial and the
# expensive vector PDFs are worth being able to queue by hand.
fig_filter <- if (exists("only_figs")) only_figs else NULL
save_fig <- function(plot, stem, width, height) {
    # match the figure key exactly -- a substring test would let "clan" select
    # "clanid" too
    if (!is.null(fig_filter) &&
        !any(stem == paste0("pfam_circlepack_", fig_filter, tag))) {
        return(invisible(NULL))
    }
    for (dev in devices) {
        f <- file.path(outdir, paste0(stem, ".", dev))
        if (dev == "pdf") {
            ggsave(f, plot, width = width, height = height,
                   device = cairo_pdf, limitsize = FALSE)
        } else {
            ggsave(f, plot, width = width, height = height,
                   dpi = 300, limitsize = FALSE)
        }
        cat(sprintf("%s: %.1f MB\n", basename(f), file.size(f) / 1024^2))
    }
}

base_theme <- theme_void(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 20, hjust = 0,
                                    margin = margin(b = 4)),
          plot.subtitle = element_text(colour = "#52514e", size = 13, hjust = 0,
                                       margin = margin(b = 12)),
          plot.margin = margin(18, 18, 18, 18),
          plot.background = element_rect(fill = "#fcfcfb", colour = NA))

cut_list  <- paste(cut_labels, collapse = ", ")
drop_note <- if (drop_largest_top) {
    paste0("\n", if (drop_top_n == 1) "Largest " else paste(drop_top_n, "largest "),
           cut_labels[1], " cluster", if (drop_top_n > 1) "s" else "",
           " (", format(sum(biggest$n), big.mark = ",", trim = TRUE),
           " families) removed from the hierarchy;\ntheir ",
           format(length(promoted), big.mark = ",", trim = TRUE), " ", cut_labels[2],
           " sub-clusters are promoted to the top level")
} else ""

# ---- 1. filled by cut level ----
p_depth <- ggraph(lay) +
    geom_node_circle(data = drawn, aes(fill = level), colour = "#ffffff", linewidth = 0.06) +
    scale_fill_manual(values = depth_fills, name = "Dendrogram cut") +
    coord_fixed() + base_theme + theme(legend.position = "right") +
    labs(title = paste0(format(nrow(pfam), big.mark = ","),
                        " Pfam families, nested by dendrogram cut height"),
         subtitle = paste0("Each orange dot is one family; enclosing circles are clusters at cuts ",
                           cut_list, drop_note))

save_fig(p_depth, paste0("pfam_circlepack_depth", tag), fig_w_leg, fig_h)

# ---- 2. tips by clan membership ----
p_clan <- ggraph(lay) +
    geom_node_circle(data = clusters, aes(linewidth = I(edge_w)),
                     fill = "#ffffff", colour = "#a9a8a1") +
    geom_node_circle(data = tips, aes(fill = in_clan), colour = NA) +
    scale_fill_manual(values = c("In a Pfam clan" = "#2a78d6",
                                 "No clan assigned" = "#b8b7b0"),
                      name = "Pfam family") +
    coord_fixed() + base_theme + theme(legend.position = "right") +
    labs(title = "Clan membership across the same nested clusters",
         subtitle = paste0("Cluster circles are outlines only; each dot is one Pfam family",
                           drop_note))

save_fig(p_clan, paste0("pfam_circlepack_clan", tag), fig_w_leg, fig_h)

# ---- 3. tips by clan identity, one random colour each ----
# 812 clans is far past any palette that stays distinguishable, so hues are
# sampled rather than assigned: the point is seeing which branches share a clan,
# not looking any one clan up. Tips get no stroke -- at this size a white
# outline swallows the fill.
# Colours come from the full clan list, not just the clans in this variant, so a
# given clan keeps its colour whichever cut set or subset is plotted.
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
clans_shown <- length(unique(na.omit(pfam$clan_acc)))

p_clanid <- ggraph(lay) +
    geom_node_circle(data = clusters, aes(linewidth = I(edge_w)),
                     fill = "#ffffff", colour = "#cbcac3") +
    geom_node_circle(data = tips, aes(fill = clan_acc), colour = NA) +
    scale_fill_manual(values = clan_colours, na.value = "#c9c8c0", guide = "none") +
    coord_fixed() + base_theme +
    labs(title = "Pfam clans across the nested clusters",
         subtitle = paste0("One random colour per clan (", clans_shown,
                           " clans shown); families with no clan assignment are grey",
                           drop_note))

save_fig(p_clanid, paste0("pfam_circlepack_clanid", tag), fig_w, fig_h)

# Cluster counts per level, for the write-up
pfam %>% summarise(across(all_of(cuts), n_distinct)) %>% print()
