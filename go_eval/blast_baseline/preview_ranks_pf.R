# Exploratory previews (2026-09-29): Fig 2C rank composition with the strict
# candidate filter (*_rank_ont_filt_pf: ancestors + GPAD-known + proFAB-known
# terms removed) and a DIAMOND nearest-neighbour baseline (blast_ranks.py,
# tie-averaged rank). Known-positive-excluded pairs. Not part of figures.R.
#   preview_C_pf.png          overall, CAV / DeepGo-SE / BLAST x ontology
#   preview_C_by_identity.png split by best identity of the protein to the true
#                             term's CAV training positives
#   preview_C_by_families.png split at 20 sequence families (term_diversity.py)
#   preview_C_by_families3.png split <20 / 20-59 / >=60 families (from the cutoff screen)
#   preview_C_by_families3_grouped.png  same, methods grouped within each ontology
# Run from bioCAV-paper/go_eval/blast_baseline/.
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

STAGE <- "../../figures/figure_data_staging/noknown"
FIG_FONT <- "Arial"
RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
RANK_COLORS <- c("Rank 1" = "#005A32", "Rank 2-3" = "#238B45", "Rank 4-10" = "#74C476",
                 "Rank >10" = "#C4B7DC", "No prediction" = "#8073AC")
METHODS <- c("CAV", "DeepGoSE", "BLAST")
ONTS <- c("mf", "bp", "cc")
ONT_DISPLAY <- c("MF", "BP", "CC")
ID_LEVELS <- c("No hit", "<30%", "30-50%", ">=50%")

tier_of <- function(rank, predicted) {
  case_when(!predicted ~ "No prediction", rank < 1.5 ~ "Rank 1", rank <= 3 ~ "Rank 2-3",
            rank <= 10 ~ "Rank 4-10", TRUE ~ "Rank >10")
}

pairs <- map_dfr(ONTS, \(o) {
  read_csv(sprintf("%s/go_specificity_ranks_allterms_%s.csv", STAGE, o), show_col_types = FALSE) |>
    inner_join(read_csv(sprintf("blast_ranks_%s.csv", o), show_col_types = FALSE),
               by = c("protein_id", "go_term")) |>
    inner_join(read_csv(sprintf("../diversity/term_diversity_%s.csv", o), show_col_types = FALSE),
               by = "go_term") |>
    mutate(ontology = toupper(o))
}) |>
  mutate(id_bin = case_when(!blast_predicted ~ "No hit", blast_true_pident < 30 ~ "<30%",
                            blast_true_pident < 50 ~ "30-50%", TRUE ~ ">=50%") |>
           factor(levels = ID_LEVELS),
         fam3 = cut(n_families, c(-Inf, 20, 60, Inf), right = FALSE,
                    labels = c("<20 families", "20-59 families", ">=60 families")),
         fam_bin = factor(if_else(n_families < 20, "<20 families", ">=20 families"),
                          levels = c("<20 families", ">=20 families")),
         all = factor("All pairs"))

long <- pairs |>
  transmute(ontology, id_bin, fam_bin, fam3, all, go_term,
            CAV      = tier_of(cav_rank_ont_filt_pf, llr > 0),
            DeepGoSE = tier_of(tool_rank_ont_filt_pf, tool_predicted == TRUE),
            BLAST    = tier_of(blast_rank_mid, blast_predicted)) |>
  pivot_longer(all_of(METHODS), names_to = "method", values_to = "tier")

rank_plot <- function(split_col, height, fname, x_lab = "Share of validation protein-GO pairs") {
  df <- long |>
    count(split = .data[[split_col]], ontology, method, tier) |>
    group_by(split, ontology, method) |> mutate(prop = n / sum(n)) |> ungroup() |>
    mutate(tier = factor(tier, levels = RANK_LEVELS),
           method = factor(method, levels = METHODS),
           ontology = factor(ontology, levels = rev(ONT_DISPLAY)))
  n_lab <- pairs |> count(split = .data[[split_col]], ontology) |>
    mutate(ontology = factor(ontology, levels = rev(ONT_DISPLAY)),
           method = factor("CAV", levels = METHODS))
  strip <- pairs |> group_by(split = .data[[split_col]]) |>
    summarise(lab = sprintf("%s\n(%d pairs)", first(as.character(split)), n()), .groups = "drop") |>
    deframe()
  p <- ggplot(df, aes(x = prop, y = ontology, fill = tier)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    geom_text(data = n_lab, aes(x = 1.01, y = ontology, label = n), inherit.aes = FALSE,
              hjust = 0, size = 2, family = FIG_FONT, color = "gray30") +
    facet_grid(method ~ split, switch = "y", labeller = labeller(split = strip)) +
    scale_fill_manual(values = RANK_COLORS, breaks = RANK_LEVELS, name = NULL) +
    scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1),
                       limits = c(0, 1.14), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_cowplot(font_size = 8, font_family = FIG_FONT) +
    labs(x = x_lab, y = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8),
          legend.position = "top", legend.justification = "left",
          legend.key.size = unit(6, "pt"), strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 8, margin = margin(l = 4, r = 6)),
          strip.text.x = element_text(size = 8), panel.spacing.x = unit(10, "pt"),
          panel.spacing.y = unit(8, "pt"))
  width <- max(4.6, 1.6 + 1.6 * n_distinct(df$split))
  ggsave(fname, p, width = min(width, 7.2), height = height, dpi = 300, bg = "white")
  ggsave(sub("png$", "pdf", fname), p, width = min(width, 7.2), height = height, device = cairo_pdf)
  df
}

d_all <- rank_plot("all", 3.4, "preview_C_pf.png")
d_id  <- rank_plot("id_bin", 3.4, "preview_C_by_identity.png")
d_fam <- rank_plot("fam_bin", 3.4, "preview_C_by_families.png")
d_fam3 <- rank_plot("fam3", 3.4, "preview_C_by_families3.png")
# Grouped by ontology: each MF/BP/CC block holds its three methods, so methods are
# compared side by side within an ontology.
rank_plot_by_ont <- function(split_col, height, fname) {
  df <- long |>
    count(split = .data[[split_col]], ontology, method, tier) |>
    group_by(split, ontology, method) |> mutate(prop = n / sum(n)) |> ungroup() |>
    mutate(tier = factor(tier, levels = RANK_LEVELS),
           method = factor(method, levels = rev(METHODS)),
           ontology = factor(ontology, levels = ONT_DISPLAY))
  n_lab <- pairs |> count(split = .data[[split_col]], ontology) |>
    mutate(ontology = factor(ontology, levels = ONT_DISPLAY),
           method = factor("CAV", levels = rev(METHODS)))
  strip <- pairs |> group_by(split = .data[[split_col]]) |>
    summarise(lab = sprintf("%s\n(%d pairs)", first(as.character(split)), n()), .groups = "drop") |>
    deframe()
  p <- ggplot(df, aes(x = prop, y = method, fill = tier)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    geom_text(data = n_lab, aes(x = 1.01, y = method, label = sprintf("n=%d", n)),
              inherit.aes = FALSE, hjust = 0, size = 2, family = FIG_FONT, color = "gray30") +
    facet_grid(ontology ~ split, switch = "y", labeller = labeller(split = strip)) +
    scale_fill_manual(values = RANK_COLORS, breaks = RANK_LEVELS, name = NULL) +
    scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1),
                       limits = c(0, 1.17), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_cowplot(font_size = 8, font_family = FIG_FONT) +
    labs(x = "Share of validation protein-GO pairs", y = NULL) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8),
          legend.position = "top", legend.justification = "left",
          legend.key.size = unit(6, "pt"), strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 8, face = "bold",
                                            margin = margin(l = 4, r = 6)),
          strip.text.x = element_text(size = 8), panel.spacing.x = unit(10, "pt"),
          panel.spacing.y = unit(8, "pt"))
  ggsave(fname, p, width = 7.2, height = height, dpi = 300, bg = "white")
  ggsave(sub("png$", "pdf", fname), p, width = 7.2, height = height, device = cairo_pdf)
}
rank_plot_by_ont("fam3", 3.4, "preview_C_by_families3_grouped.png")


summ <- function(d) d |> group_by(split, ontology, method) |>
  summarise(rank1 = round(100 * sum(prop[tier == "Rank 1"]), 1),
            top10 = round(100 * sum(prop[tier %in% RANK_LEVELS[1:3]]), 1), .groups = "drop") |>
  pivot_wider(names_from = method, values_from = c(rank1, top10)) |> arrange(split, desc(ontology))
print(summ(d_id), n = 20, width = 200)
print(summ(d_fam), n = 20, width = 200)
print(summ(d_fam3), n = 30, width = 200)
print(pairs |> count(fam3, ontology, name = "pairs") |>
  left_join(pairs |> distinct(fam3, ontology, go_term) |> count(fam3, ontology, name = "terms"),
            by = c("fam3", "ontology")))
