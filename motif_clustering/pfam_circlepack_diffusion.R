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
# swallow the figure. How many of the largest to leave undrawn (0 draws all).
# The families inside them are kept and packed exactly as before -- only the
# outline of each such cluster, and of its children at the next cut, is
# suppressed, so the sub-structure reads without the giant rings around it.
hide_top_n <- 3
# ----------------------------------

hide_largest_top <- hide_top_n > 0
# The packing covers the full data whatever is hidden, so every variant of a cut
# set shares one layout cache; only the output figures get the hide tag.
layout_tag <- paste0("_", sub("cut_", "", cuts[1]))
tag <- paste0(layout_tag, if (hide_largest_top) paste0("_hidetop", hide_top_n) else "")

pfam_all <- read_csv(infile, show_col_types = FALSE) %>%
    rename(ID = accession) %>%
    select(ID, clan_acc, clan_name, pfam_short_name, all_of(cuts))

pfam <- pfam_all   # no families are dropped any more; the two stay in step

# Cluster nodes to leave undrawn: the n largest at the top cut, plus their
# children at the next cut down. Names match the clusterid built below.
hidden_nodes <- character(0)
if (hide_largest_top) {
    biggest <- pfam %>% count(.data[[cuts[1]]], sort = TRUE) %>% slice(seq_len(hide_top_n))
    inside <- pfam %>% filter(.data[[cuts[1]]] %in% biggest[[1]])
    hidden_kids <- sort(unique(inside[[cuts[2]]]))
    hidden_nodes <- c(paste0(biggest[[1]], cuts[1]), paste0(hidden_kids, cuts[2]))
    cat(sprintf("hiding %d largest %s cluster(s) [%s] plus their %d %s children; %d of %d families (%.1f%%) sit inside and are still drawn\n",
                hide_top_n, cuts[1], paste(biggest[[1]], collapse = ", "),
                length(hidden_kids), cuts[2],
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

# All nodes hang off one origin so the layout has a single root
clusters_uniqued$link <- "origin"

# Edges in from/to order, one level at a time, tips last
levels_ordered <- c("link", cuts, "ID")
edge_list <- map_dfr(seq_len(length(levels_ordered) - 1), function(i) {
    clusters_uniqued %>%
        select(from = all_of(levels_ordered[i]), to = all_of(levels_ordered[i + 1])) %>%
        distinct()
})

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
cut_labels <- sprintf("k = %s", formatC(as.integer(sub("cut_k0*", "", cuts)), big.mark = ",", format = "d"))
depth_labels <- c(cut_labels, "Pfam family")

# Sequential blue for the nested cuts, orange tips so single families stay visible
blue_steps <- if (n_cut == 5) {
    c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281")
} else {
    colorRampPalette(c("#86b6ef", "#104281"))(n_cut)
}
depth_fills <- setNames(c(blue_steps, "#eb6834"), depth_labels)

# Outlines thin out as circles get smaller
w_ramp <- c(0.35, 0.25, 0.18, 0.12, 0.08, 0.04)
lay <- lay %>%
    mutate(level  = factor(ifelse(depth >= 1, depth_labels[depth], NA), levels = depth_labels),
           edge_w = w_ramp[pmin(pmax(depth, 1), length(w_ramp))])

drawn    <- lay %>% filter(depth >= 1,      # hides the single origin circle
                           !name %in% hidden_nodes)
clusters <- drawn %>% filter(depth <= n_cut)
tips     <- drawn %>% filter(depth == n_cut + 1)

# Size the canvas from the packing's own bounding box, so coord_fixed doesn't
# leave gutters down the sides. Different cut sets pack to different aspects.
bb     <- with(drawn, c(min(x - r), max(x + r), min(y - r), max(y + r)))
aspect <- (bb[2] - bb[1]) / (bb[4] - bb[3])
panel_h  <- 12
panel_w  <- panel_h * aspect
fig_h    <- panel_h + 1.2                    # room for title and subtitle
fig_w    <- panel_w + 0.5                    # margins only
fig_w_leg <- panel_w + 2.2                   # margins plus the legend column
cat(sprintf("content aspect %.3f -> figure %.1f x %.1f in\n", aspect, fig_w, fig_h))

base_theme <- theme_void(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 20, hjust = 0,
                                    margin = margin(b = 4)),
          plot.subtitle = element_text(colour = "#52514e", size = 13, hjust = 0,
                                       margin = margin(b = 12)),
          plot.margin = margin(18, 18, 18, 18),
          plot.background = element_rect(fill = "#fcfcfb", colour = NA))

cut_list  <- paste(sub("cut ", "", cut_labels), collapse = ", ")
drop_note <- if (hide_largest_top) {
    paste0("\n", if (hide_top_n == 1) "Largest " else paste(hide_top_n, "largest "),
           cut_labels[1], " cluster", if (hide_top_n > 1) "s" else "",
           " (", sum(biggest$n), " families) drawn without their outlines, ",
           "and without those of their ", cut_labels[2], " sub-clusters")
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

ggsave(file.path(outdir, paste0("pfam_circlepack_depth", tag, ".png")), p_depth,
       width = fig_w_leg, height = fig_h, dpi = 300, limitsize = FALSE)

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

ggsave(file.path(outdir, paste0("pfam_circlepack_clan", tag, ".png")), p_clan,
       width = fig_w_leg, height = fig_h, dpi = 300, limitsize = FALSE)

# ---- 3. tips by clan identity, one random colour each ----
# 812 clans is far past any palette that stays distinguishable, so hues are
# sampled rather than assigned: the point is seeing which branches share a clan,
# not looking any one clan up. Tips get no stroke -- at this size a white
# outline swallows the fill.
# Colours come from the full clan list, not just the clans in this variant, so a
# given clan keeps its colour whichever cut set or subset is plotted.
set.seed(42)
clan_ids <- sort(unique(na.omit(pfam_all$clan_acc)))
clans_shown <- length(unique(na.omit(pfam$clan_acc)))
clan_colours <- setNames(
    hcl(h = sample(seq(0, 360, length.out = length(clan_ids))),
        c = runif(length(clan_ids), 70, 120),
        l = runif(length(clan_ids), 35, 65)),
    clan_ids)

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

ggsave(file.path(outdir, paste0("pfam_circlepack_clanid", tag, ".png")), p_clanid,
       width = fig_w, height = fig_h, dpi = 300, limitsize = FALSE)

# Cluster counts per level, for the write-up
pfam %>% summarise(across(all_of(cuts), n_distinct)) %>% print()
