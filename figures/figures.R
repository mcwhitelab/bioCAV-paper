#!/usr/bin/env Rscript
# figures.R — Paper figures from pipeline-exported CSVs
#
# Prerequisites:
#   install.packages(c("tidyverse", "cowplot", "ggridges", "ggrepel"))
#
# Usage:
#   Rscript figures.R
#   # or interactively: source("figures.R")
#
# Data source:
#   Run pipeline scripts with --figure-data-dir figure_data --label <mf|bp|cc>
#   to populate figure_data/ before running this script.

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggridges)
  library(ggrepel)
  library(RColorBrewer)
  library(ggrastr)
  library(png)
  library(grid)
})

RASTER_DPI <- 150

DATA <- "figure_data"
OUT  <- "figures"
dir.create(OUT, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# House font: Arial isn't installed on this system (no Arial in `fc-list`),
# so we use Liberation Sans -- a metrically-compatible Arial substitute
# (same glyph widths) -- for every text element. Setting the geom defaults
# here covers every geom_text()/geom_text_repel()/annotate("text", ...) call
# in this file (and in draft_fig_extra_proteins_common.R, sourced later into
# the same session) without having to touch each call site individually.
# ---------------------------------------------------------------------------
FIG_FONT <- "Arial"
update_geom_defaults("text",       list(family = FIG_FONT))
update_geom_defaults("text_repel", list(family = FIG_FONT))

# Okabe-Ito palette (colorblind-friendly)
oi <- c(
  orange     = "#E69F00",
  sky_blue   = "#56B4E9",
  green      = "#009E73",
  yellow     = "#F0E442",
  blue       = "#0072B2",
  vermillion = "#D55E00",
  pink       = "#CC79A7",
  black      = "#000000"
)

CAV_COLOR  <- unname(oi["blue"])
TOOL_COLOR <- unname(oi["vermillion"])

# ---------------------------------------------------------------------------
# Helper: base theme
# ---------------------------------------------------------------------------
# theme_cowplot scales axis/legend/strip text down from `font_size` (axis.text
# lands at 6.86pt, legend/strip at rel(0.857)); override them back to a flat 8pt
# so tick labels match the axis titles. Panels that set their own sizes after
# base_theme() still win.
base_theme <- function(...) {
  theme_cowplot(font_size = 8, font_family = FIG_FONT, ...) +
    theme(
      axis.text   = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text  = element_text(size = 8)
    )
}

# ---------------------------------------------------------------------------
# Helper: load a per-ontology file, return NULL if missing
# ---------------------------------------------------------------------------
load_ont <- function(stem, ont) {
  f <- file.path(DATA, paste0(stem, "_", ont, ".csv"))
  if (!file.exists(f)) {
    message("Missing: ", f)
    return(NULL)
  }
  read_csv(f, show_col_types = FALSE) |> mutate(ontology = toupper(ont))
}

# ===========================================================================
# Figure 2 — GO evaluation
# ===========================================================================

# ---------------------------------------------------------------------------
# Panel A: intro panel -- CAV score (z-scored against the negative/
# background distribution) for held-out positive vs. negative examples,
# pooled across ontologies. This is the "does the method work at all"
# panel: negatives sit centered near z=0 (they define the reference
# frame), positives shift well above it. Note the z uses ONE pooled
# negative distribution per source file, not each GO term's own: upstream
# (compare_tool_temporal.py) pools raw per-term cav_score values across
# terms and the CSVs carry no term column, so a per-term z is not
# recoverable here. Source: figure_data/
# temporal_pos_neg_density_*.csv (one file per ontology, written by
# summarize_temporal_eval.py; filenames carry an external-tool mAP label,
# not the ontology name, so we just glob and pool rather than matching by
# ontology).
# ---------------------------------------------------------------------------
pos_neg_files <- list.files(DATA, pattern = "^temporal_pos_neg_density_.*\\.csv$",
                            full.names = TRUE)

p_2a_zscore <- NULL
if (length(pos_neg_files) > 0) {
  zscore_df <- map_dfr(pos_neg_files, function(f) {
    d <- read_csv(f, show_col_types = FALSE)
    neg <- d |> filter(label == "negative") |> pull(cav_score)
    d |> mutate(z = (cav_score - mean(neg)) / sd(neg))
  }) |>
    mutate(label = factor(label, levels = c("negative", "positive"),
                          labels = c("Background (negative)", "Positive")))

  ZSCORE_COLORS <- c("Background (negative)" = "gray50", "Positive" = CAV_COLOR)

  p_2a_zscore <- zscore_df |>
    ggplot(aes(x = z, fill = label, color = label)) +
    geom_histogram(aes(y = after_stat(density)), position = "identity",
                   bins = 60, alpha = 0.55, linewidth = 0) +
    geom_vline(xintercept = 0, color = "gray40", linetype = "dashed", linewidth = 0.4) +
    scale_fill_manual(values = ZSCORE_COLORS, name = NULL) +
    scale_color_manual(values = ZSCORE_COLORS, name = NULL) +
    scale_x_continuous(limits = c(-5, 15)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    base_theme() +
    labs(x = "CAV score (z, vs. pooled\nnegative background)", y = "Density") +
    guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1)) +
    theme(
      legend.position = "top",
      legend.key.size = unit(7, "pt"),
      legend.text     = element_text(size = 8),
      legend.margin   = margin(t = 0, b = 0),
      plot.margin     = margin(t = 2, r = 4, b = 2, l = 2)
    )
} else {
  message("Skipping Figure 2 intro panel: no figure_data/temporal_pos_neg_density_*.csv found")
}

# Load and combine all three ontologies (used by panels C, D, E)
ont_comp <- bind_rows(
  load_ont("temporal_tool_comparison", "mf"),
  load_ont("temporal_tool_comparison", "bp"),
  load_ont("temporal_tool_comparison", "cc")
)

# Ordered factor for ontology rows: MF on top, CC middle, BP bottom in ggridges
# (ggridges maps factor levels bottom→top on the y-axis)
ONT_LEVELS <- c("BP", "CC", "MF")

ridges_opts <- list(
  scale          = 0.85,
  rel_min_height = 0.01,
  linewidth      = 0.4,
  alpha          = 0.55
)

if (nrow(ont_comp) > 0) {

  # Method identity (violins) and rank-tier ramp (Panel A)
  METHOD_LEVELS <- c("CAV", "DeepGoSE")
  METHOD_COLORS <- c(CAV = CAV_COLOR, DeepGoSE = TOOL_COLOR)
  ONT_DISPLAY   <- c("MF", "BP", "CC")

  # Rank tiers are an ordered composition -> single-hue sequential green ramp
  # (dark green = best rank), except the two "barely found" tiers (Rank >10,
  # No prediction), which shade into purple instead of pale green so they
  # read as a distinct "not really found" family from the ranked greens.
  RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
  RANK_COLORS <- c(
    "Rank 1"        = "#005A32",
    "Rank 2-3"      = "#238B45",
    "Rank 4-10"     = "#74C476",
    "Rank >10"      = "#C4B7DC",
    "No prediction" = "#8073AC"
  )

  # ---------------------------------------------------------------------------
  # Panel A: rank-tier composition (100% stacked bars, one per method x ontology)
  # For each (protein, true GO term) pair, classify where the true term ranks
  # among all trained CAVs. CAV "not found" = LLR <= 0 (low confidence);
  # DeepGoSE "not found" = tool did not predict the term. The rank data is
  # bimodal (mostly rank 1 or not found), so a composition bar reads far better
  # than a histogram (which was mostly empty in the middle).
  # ---------------------------------------------------------------------------
  ranks_all <- bind_rows(
    load_ont("go_specificity_ranks", "mf"),
    load_ont("go_specificity_ranks", "bp"),
    load_ont("go_specificity_ranks", "cc")
  )

  if (nrow(ranks_all) > 0) {
    tier_of <- function(rank, predicted) {
      case_when(
        !predicted ~ "No prediction",
        rank == 1  ~ "Rank 1",
        rank <= 3  ~ "Rank 2-3",
        rank <= 10 ~ "Rank 4-10",
        TRUE       ~ "Rank >10"
      )
    }

    rankA <- ranks_all |>
      transmute(
        ontology,
        CAV      = tier_of(cav_rank,  llr > 0),
        DeepGoSE = tier_of(tool_rank, tool_predicted == TRUE)
      ) |>
      pivot_longer(c(CAV, DeepGoSE), names_to = "method", values_to = "tier") |>
      count(ontology, method, tier) |>
      group_by(ontology, method) |>
      mutate(prop = n / sum(n)) |>
      ungroup() |>
      mutate(
        tier     = factor(tier,     levels = RANK_LEVELS),
        method   = factor(method,   levels = METHOD_LEVELS),        # CAV facet on top
        ontology = factor(ontology, levels = rev(ONT_DISPLAY))      # MF on top within facet
      )

    p_2a <- rankA |>
      ggplot(aes(x = prop, y = ontology, fill = tier)) +
      geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
      facet_grid(method ~ ., switch = "y") +
      scale_fill_manual(values = RANK_COLORS, breaks = RANK_LEVELS, name = NULL) +
      scale_x_continuous(labels = scales::percent,
                         expand = expansion(mult = c(0, 0.02))) +
      base_theme() +
      labs(x = "Share of validation protein-GO pairs", y = NULL) +
      guides(fill = guide_legend(nrow = 1)) +
      theme(
        legend.position   = "top",
        legend.location   = "plot",
        legend.justification = "left",
        legend.key.size   = unit(6, "pt"),
        legend.text       = element_text(size = 8),
        legend.margin     = margin(t = 0, b = 2, r = 0, l = 16),
        legend.spacing.x  = unit(0, "pt"),
        legend.box.spacing = unit(2, "pt"),
        strip.background  = element_blank(),
        strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0, size = 8, margin = margin(l = 4, r = 6)),
        axis.text.y       = element_text(margin = margin(r = 4)),
        panel.spacing.y   = unit(10, "pt"),
        plot.margin       = margin(t = 2, r = 2, b = 2, l = 2)
      )
  } else {
    message("Skipping rank composition: no go_specificity_ranks_*.csv files found")
    p_2a <- NULL
  }

  # ---------------------------------------------------------------------------
  # Panels B & C: dodged violins (CAV vs DeepGoSE) per ontology, y fixed [0,1]
  # ---------------------------------------------------------------------------
  make_violin <- function(cav_col, tool_col, y_lab) {
    ont_comp |>
      transmute(ontology,
                CAV      = .data[[cav_col]],
                DeepGoSE = .data[[tool_col]]) |>
      pivot_longer(c(CAV, DeepGoSE), names_to = "method", values_to = "val") |>
      drop_na(val) |>
      mutate(ontology = factor(ontology, levels = ONT_DISPLAY),
             method   = factor(method,   levels = METHOD_LEVELS)) |>
      ggplot(aes(x = ontology, y = val, fill = method, color = method)) +
      geom_violin(position = position_dodge(width = 0.8), width = 0.75,
                  alpha = 0.5, linewidth = 0.3, scale = "width",
                  draw_quantiles = 0.5) +
      scale_fill_manual(values  = METHOD_COLORS, name = NULL) +
      scale_color_manual(values = METHOD_COLORS, name = NULL) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                         breaks = seq(0, 1, 0.25)) +
      base_theme() +
      labs(x = NULL, y = y_lab) +
      theme(
        legend.position    = "none",
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
      )
  }

  p_2b <- make_violin("auc_val_vs_test_neg",  "tool_auc",  "AUC")
  p_2c <- make_violin("aupr_val_vs_test_neg", "tool_aupr", "AUPR")

  # ---------------------------------------------------------------------------
  # Panel D: CAV vs CLEAN-SupCon on EC specificity (protein, EC term) pairs.
  # CLEAN-SupCon emits a single predicted EC per protein rather than a scored
  # candidate list, so it has no "rank" of the true term -- only Correct /
  # Incorrect. CAV keeps its rank tiers. These are different concepts from
  # "No prediction" (CAV made no confident call at all, LLR <= 0), so
  # Correct/Incorrect get their own colors rather than reusing the rank ramp
  # or the "No prediction" purple.
  # ---------------------------------------------------------------------------
  ec_ranks_path <- file.path(DATA, "ec_specificity_ranks.csv")

  p_2d <- NULL
  if (file.exists(ec_ranks_path)) {
    ec_ranks <- read_csv(ec_ranks_path, show_col_types = FALSE)

    EC_TIER_LEVELS <- c(RANK_LEVELS, "Correct", "Incorrect")
    EC_TIER_COLORS <- c(
      RANK_COLORS,
      # Not blue/orange -- those already mean CAV/DeepGoSE in panels B & C,
      # so reusing them here for CLEAN-SupCon would misleadingly imply the
      # same encoding.
      "Correct"   = unname(oi["pink"]),
      "Incorrect" = "#4D4D4D"
    )

    EC_LEVEL_LABELS <- c("3" = "EC level 3", "4" = "EC level 4\n(fully specific)")

    rankD <- bind_rows(
      ec_ranks |>
        transmute(level, method = "CAV",
                  tier   = tier_of(cav_rank, llr > 0)),
      ec_ranks |>
        transmute(level, method = "CLEAN-SupCon",
                  tier   = if_else(clean_supcon_correct, "Correct", "Incorrect"))
    ) |>
      count(level, method, tier) |>
      group_by(level, method) |>
      mutate(prop = n / sum(n), n_pairs = sum(n)) |>
      ungroup() |>
      mutate(
        tier   = factor(tier,   levels = EC_TIER_LEVELS),
        method = factor(method, levels = c("CLEAN-SupCon", "CAV")),  # CAV on top
        level  = factor(EC_LEVEL_LABELS[as.character(level)],
                        levels = EC_LEVEL_LABELS[c("4", "3")])        # level 4 on top
      )

    n_ec_pairs <- nrow(ec_ranks)

    p_2d <- rankD |>
      ggplot(aes(x = prop, y = method, fill = tier)) +
      geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
      facet_grid(level ~ ., switch = "y") +
      scale_fill_manual(values = EC_TIER_COLORS, breaks = EC_TIER_LEVELS, name = NULL) +
      scale_x_continuous(labels = scales::percent,
                         expand = expansion(mult = c(0, 0.02))) +
      base_theme() +
      labs(x = paste0("Share of protein-EC pairs (n=", n_ec_pairs, ")"),
           y = NULL) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      theme(
        legend.position   = "top",
        legend.location   = "plot",
        legend.justification = "left",
        legend.key.size   = unit(6, "pt"),
        legend.text       = element_text(size = 8),
        legend.margin     = margin(t = 0, b = 2, r = 0, l = 16),
        legend.spacing.x  = unit(0, "pt"),
        strip.background  = element_blank(),
        strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0, size = 8),
        panel.spacing.y   = unit(4, "pt"),
        plot.margin       = margin(t = 2, r = 6, b = 2, l = 2)
      )
  } else {
    message("Skipping EC rank composition: figure_data/ec_specificity_ranks.csv not found")
  }

  # ---------------------------------------------------------------------------
  # Panel E: cross-ontology validation — cosine similarity between a GO
  # term's CAV and its EC-annotated partner's CAV (curated GO<->EC mapping,
  # ec2go), vs. the background similarity to all other trained EC CAVs.
  # GO and EC CAVs are trained completely independently (different label
  # sources, no shared supervision), so this tests whether they converge on
  # the same biological structure.
  # ---------------------------------------------------------------------------
  ge_pairs_path <- file.path(DATA, "go_ec_cosine_pairs.csv")
  ge_bg_path    <- file.path(DATA, "go_ec_cosine_background_sample.csv")

  p_2e <- NULL
  if (file.exists(ge_pairs_path) && file.exists(ge_bg_path)) {
    ge_pairs <- read_csv(ge_pairs_path, show_col_types = FALSE)
    ge_bg    <- read_csv(ge_bg_path,    show_col_types = FALSE)

    n_pairs   <- nrow(ge_pairs)
    n_ec_cavs <- ge_pairs$n_ec_cavs[1]
    pct_rank1 <- mean(ge_pairs$rank == 1) * 100

    GE_MATCH_COLOR <- CAV_COLOR
    GE_BG_COLOR    <- "gray60"

    hist_df <- bind_rows(
      ge_bg    |> transmute(category = "Background",     sim = sim),
      ge_pairs |> transmute(category = "Annotated match", sim = matched_sim)
    ) |>
      mutate(category = factor(category, levels = c("Background", "Annotated match")))

    med_bg    <- median(ge_bg$sim)
    med_match <- median(ge_pairs$matched_sim)

    GE_LEGEND_THEME <- theme(
      legend.position = "top",
      legend.key.size = unit(7, "pt"),
      legend.text     = element_text(size = 8),
      legend.margin   = margin(t = 0, b = 0),
      plot.margin     = margin(t = 2, r = 4, b = 2, l = 2)
    )

    # Panel E: cosine-similarity histogram.
    p_2e_solo <- hist_df |>
      ggplot(aes(x = sim, fill = category, color = category)) +
      geom_histogram(aes(y = after_stat(density)), position = "identity",
                     bins = 40, alpha = 0.6, linewidth = 0) +
      geom_vline(xintercept = med_bg,    color = GE_BG_COLOR,    linetype = "dashed", linewidth = 0.4) +
      geom_vline(xintercept = med_match, color = GE_MATCH_COLOR, linetype = "dashed", linewidth = 0.4) +
      scale_fill_manual(values  = c("Background" = GE_BG_COLOR, "Annotated match" = GE_MATCH_COLOR), name = NULL) +
      scale_color_manual(values = c("Background" = GE_BG_COLOR, "Annotated match" = GE_MATCH_COLOR), name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      base_theme() +
      labs(x = "Cosine similarity (GO-CAV vs. EC-CAV)", y = "Density") +
      GE_LEGEND_THEME

    max_rank_show <- 20
    rank_df <- ge_pairs |>
      mutate(rank_bin = if_else(rank > max_rank_show, max_rank_show + 1L, as.integer(rank))) |>
      count(rank_bin) |>
      complete(rank_bin = 1:(max_rank_show + 1), fill = list(n = 0)) |>
      mutate(rank_label = if_else(rank_bin > max_rank_show, paste0(">", max_rank_show), as.character(rank_bin)),
             rank_label = factor(rank_label, levels = c(as.character(1:max_rank_show), paste0(">", max_rank_show))),
             category   = "Annotated match")

    # Panel F: rank histogram. Same fill/alpha/legend treatment as E (single
    # "Annotated match" category, same legend row position/size) so the two
    # read as a matched pair; the rank-1 stat moves to an in-plot annotation
    # instead of a title, so it doesn't compete with cowplot's "F" label.
    p_2f_solo <- rank_df |>
      ggplot(aes(x = rank_label, y = n, fill = category)) +
      geom_col(alpha = 0.6, width = 0.75) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("Rank 1: %.0f%% of\n%d GO terms", pct_rank1, n_pairs),
               size = 8 / .pt, hjust = 1.05, vjust = 1.3) +
      scale_fill_manual(values = c("Annotated match" = GE_MATCH_COLOR), name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      base_theme() +
      labs(x = paste0("Rank of annotated EC partner\n(out of ", n_ec_cavs, " trained EC CAVs)"),
           y = "Number of GO terms") +
      scale_x_discrete(breaks = c("1", "5", "10", "15", paste0(">", max_rank_show))) +
      GE_LEGEND_THEME +
      theme(axis.text.x = element_text(size = 8))

    # align = "h", axis = "tb": match the top/bottom plot-panel edges of E
    # and F so their x-axes sit at the same height despite F's rotated tick
    # labels taking more vertical space than E's horizontal ones.
    p_2e <- plot_grid(p_2e_solo, p_2f_solo, nrow = 1, align = "h", axis = "tb",
                      labels = c("E", "F"), label_size = 8, label_fontfamily = FIG_FONT,
                      rel_widths = c(1, 1))
  } else {
    message("Skipping GO-EC cosine validation: figure_data/go_ec_cosine_*.csv not found")
  }

  # ---------------------------------------------------------------------------
  # Figure 3, panel B: CAV projection score along the EC hierarchy (level
  # 1->4) for a handful of validation proteins with a full lineage. Not a
  # cosine-similarity metric -- shows the score sharpening as the CAV gets
  # more specific down the hierarchy. Source: hierarchy_decay/results/
  # ec_hierarchy_decay.csv (copied to figure_data/). Protein IDs are labeled
  # via the "invisible 1D scatter" trick: a second, axis-free ggplot with
  # one tiny point per protein at its level-4 (rightmost) score, y-aligned
  # to the main panel via matching coord_cartesian(ylim=), with
  # geom_text_repel pushing the label off to the right of that point. This
  # keeps the labels out of the line plot itself (no crowding among the
  # lines) while still reading as "labels for the right end of each line."
  # ---------------------------------------------------------------------------
  ec_decay_path <- file.path(DATA, "ec_hierarchy_decay.csv")

  p_3b <- NULL
  if (file.exists(ec_decay_path)) {
    ec_decay <- read_csv(ec_decay_path, show_col_types = FALSE)

    decay_level_labels <- c("1", "2", "3", "4")
    decay_ends <- ec_decay |> group_by(protein_id) |> slice_max(level, n = 1, with_ties = FALSE) |> ungroup()
    decay_y_rng  <- range(ec_decay$score, na.rm = TRUE)
    decay_y_pad  <- diff(decay_y_rng) * 0.06
    decay_ylim   <- c(decay_y_rng[1] - decay_y_pad, decay_y_rng[2] + decay_y_pad)
    # Okabe-Ito rather than ggplot's default hue_pal, and a solid/dashed
    # alternation so the lines stay separable where they cross (and in
    # grayscale). Yellow and black are skipped: yellow is too low-contrast
    # on white for a 0.5pt line.
    decay_ids       <- sort(unique(ec_decay$protein_id))
    decay_colors    <- setNames(
      # vermillion before orange, so the alternation below puts the two most
      # confusable hues on different line styles rather than both on dashed.
      rep_len(unname(oi[c("blue", "vermillion", "orange", "green", "pink", "sky_blue")]),
              length(decay_ids)),
      decay_ids
    )
    decay_linetypes <- setNames(rep_len(c("solid", "dashed"), length(decay_ids)), decay_ids)

    p_3b_main <- ec_decay |>
      ggplot(aes(x = level, y = score, group = protein_id, color = protein_id,
                 linetype = protein_id)) +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      geom_point(size = 1.2) +
      scale_color_manual(values = decay_colors, guide = "none") +
      scale_linetype_manual(values = decay_linetypes, guide = "none") +
      # A small expansion, not expand = c(0, 0): with no padding the level-1
      # and level-4 markers sit exactly on the panel edge and get clipped in half.
      scale_x_continuous(breaks = 1:4, labels = decay_level_labels,
                         expand = expansion(mult = 0.045)) +
      coord_cartesian(ylim = decay_ylim) +
      base_theme() +
      labs(x = "EC hierarchy depth\n(4 = fully specific)", y = "CAV projection score") +
      theme(legend.position = "none", plot.margin = margin(t = 4, r = 2, b = 4, l = 4))

    p_3b_labels <- decay_ends |>
      ggplot(aes(x = 0, y = score, color = protein_id)) +
      geom_point(alpha = 0) +
      geom_text_repel(aes(label = protein_id),
                      direction = "y", hjust = 0, nudge_x = 0.06, xlim = c(0.03, NA),
                      segment.color = NA,
                      size = 8 / .pt, fontface = "bold", show.legend = FALSE,
                      force = 8, box.padding = 0.3, seed = 42) +
      scale_color_manual(values = decay_colors, guide = "none") +
      scale_x_continuous(limits = c(-0.05, 1), expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      coord_cartesian(ylim = decay_ylim, clip = "off") +
      theme_void() +
      theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 0))

    p_3b <- plot_grid(p_3b_main, p_3b_labels, nrow = 1, align = "h", axis = "tb",
                       rel_widths = c(1, 0.58))
    # Spacer below so B's plot is shorter than A's without changing the row height.
    p_3b <- plot_grid(p_3b, NULL, ncol = 1, rel_heights = c(1, 0.16))
  } else {
    message("Skipping Figure 3 EC-hierarchy panel: figure_data/ec_hierarchy_decay.csv not found")
  }

  # ---------------------------------------------------------------------------
  # Figure 3, panel E: combined 2D embedding of all trained CAV concept
  # vectors -- EC + GO's three namespaces (MF/BP/CC) -- in one shared space,
  # colored by category. Source: figure_data/go_ec_combined_umap.csv, written
  # by motif_clustering/combined_cav_umap_3d.py's --coords-csv-out.
  # ---------------------------------------------------------------------------
  combined_umap_path <- file.path(DATA, "go_ec_combined_umap.csv")

  p_3e <- NULL
  if (file.exists(combined_umap_path)) {
    combined_umap <- read_csv(combined_umap_path, show_col_types = FALSE)

    COMBINED_CAT_COLORS <- c(EC = unname(oi["vermillion"]), MF = unname(oi["blue"]),
                             BP = unname(oi["green"]), CC = unname(oi["pink"]))
    COMBINED_CAT_LEVELS <- c("EC", "MF", "BP", "CC")

    cat_counts <- combined_umap |> count(category)
    cat_labels <- setNames(
      paste0(cat_counts$category, " (n=", cat_counts$n, ")"),
      cat_counts$category
    )

    p_3e <- combined_umap |>
      mutate(category = factor(category, levels = COMBINED_CAT_LEVELS)) |>
      ggplot(aes(D1, D2, color = category)) +
      geom_point_rast(size = 0.35, alpha = 0.55, raster.dpi = RASTER_DPI) +
      scale_color_manual(values = COMBINED_CAT_COLORS, labels = cat_labels, name = NULL) +
      guides(color = guide_legend(nrow = 2, override.aes = list(size = 2.5, alpha = 1))) +
      base_theme() +
      labs(x = "UMAP 1", y = "UMAP 2",
           title = paste0("Combined CAV space (n=", nrow(combined_umap), ")")) +
      theme(plot.title = element_text(size = 8, face = "bold"),
            legend.position = "bottom", legend.key.size = unit(8, "pt"),
            legend.text = element_text(size = 8), axis.text = element_blank(),
            axis.ticks = element_blank())
  } else {
    message("Skipping Figure 3 combined-embedding panel: figure_data/go_ec_combined_umap.csv not found ",
            "(run motif_clustering/combined_cav_umap_3d.py with --coords-csv-out first)")
  }

  # ---------------------------------------------------------------------------
  # Assemble Figure 2 — GO evaluation. Stacked rows: A blank (reserved), B
  # (intro z-score panel); C (rank-tier composition); D & E (AUC/AUPR violins).
  # ---------------------------------------------------------------------------
  intro_row <- plot_grid(
    NULL,
    if (!is.null(p_2a_zscore)) p_2a_zscore else ggplot() + theme_void(),
    nrow       = 1,
    labels     = c("A", "B"),
    label_size = 8,
    label_fontfamily = FIG_FONT
  )

  legend_method <- get_legend(
    p_2b + theme(legend.position = "bottom", legend.key.size = unit(9, "pt"))
  )

  de_row <- plot_grid(
    p_2b, p_2c,
    nrow       = 1,
    labels     = c("D", "E"),
    label_size = 8,
    label_fontfamily = FIG_FONT
  )

  # Legend sits under panel D specifically (not centered under D+E), shifted
  # right within D's column via a leading spacer.
  legend_under_d <- plot_grid(NULL, legend_method, nrow = 1, rel_widths = c(0.3, 0.7))
  legend_row     <- plot_grid(legend_under_d, NULL, nrow = 1, rel_widths = c(1, 1))

  violin_block <- plot_grid(
    de_row, legend_row,
    ncol        = 1,
    rel_heights = c(1, 0.12)
  )

  rows        <- list(intro_row)
  row_labels  <- c("")   # intro_row already carries its own A/B labels
  rel_heights <- c(0.8)

  if (!is.null(p_2a)) {
    rows        <- c(rows, list(p_2a))
    row_labels  <- c(row_labels, "C")
    rel_heights <- c(rel_heights, 0.85)
  }
  rows        <- c(rows, list(violin_block))
  row_labels  <- c(row_labels, "")
  rel_heights <- c(rel_heights, 0.72)

  fig2 <- plot_grid(
    plotlist    = rows,
    ncol        = 1,
    labels      = row_labels,
    label_size  = 8,
    label_fontfamily = FIG_FONT,
    rel_heights = rel_heights
  )

  # 3.6in x 5.56in -- the figure at 4/5 of its previous 4.5 x 7.71in
  # footprint, with the D/E violin row given less vertical space than the
  # rows above it (rel_height 1 -> 0.72; the height above is the 0.8-scaled
  # canvas reduced to match the smaller rel_heights sum, so only the violin
  # row shrinks). Font sizes stay in points, so text reads relatively larger
  # on the smaller canvas.
  ggsave(file.path(OUT, "cav_fig2.pdf"), fig2, width = 3.6, height = 5.13, device = cairo_pdf)
  ggsave(file.path(OUT, "cav_fig2.png"), fig2, width = 3.6, height = 5.13, dpi = 300, bg = "white")
  message("Saved cav_fig2.pdf / cav_fig2.png")

  # ---------------------------------------------------------------------------
  # Assemble Figure 3 — EC evaluation & cross-ontology validation. Two rows:
  #   row 1: A (EC rank composition, 2/3 width) | B (EC hierarchy depth, 1/3)
  #   row 2: C (GO-EC cosine) | D (rank of EC partner) | E (combined CAV space)
  # Missing panels become void placeholders rather than being dropped, so the
  # rel_widths keep meaning if a source CSV is absent.
  # ---------------------------------------------------------------------------
  void_panel <- ggplot() + theme_void()

  fig3_row1 <- plot_grid(
    if (!is.null(p_2d)) p_2d else void_panel,
    if (!is.null(p_3b)) p_3b else void_panel,
    nrow       = 1,
    rel_widths = c(1.55, 1),
    labels     = c("A", "B"),
    label_size = 8,
    label_fontfamily = FIG_FONT
  )

  fig3_row2 <- plot_grid(
    if (!is.null(p_2e_solo)) p_2e_solo else void_panel,
    if (!is.null(p_2f_solo)) p_2f_solo else void_panel,
    if (!is.null(p_3e))      p_3e      else void_panel,
    nrow       = 1,
    rel_widths = c(1, 1, 1),
    labels     = c("C", "D", "E"),
    label_size = 8,
    label_fontfamily = FIG_FONT
  )

  if (!is.null(p_2d) || !is.null(p_3b) || !is.null(p_2e) || !is.null(p_3e)) {
    fig3 <- plot_grid(
      fig3_row1, fig3_row2,
      ncol        = 1,
      labels      = c("", ""),   # rows carry their own panel labels
      rel_heights = c(2.2, 1.9)
    )

    # 7.2in ~ 183mm -- double-column journal width. A two-row layout with
    # three panels across row 2 needs it: at the old 4.0in single-column
    # width, D's rank axis and E's UMAP legend have nowhere to go. Height is
    # set explicitly rather than via the old sum(rel_heights) * 1.826
    # multiplier, which was derived for the previous 4-row single-column
    # stack and does not carry over.
    ggsave(file.path(OUT, "cav_fig3.pdf"), fig3, width = 7.2, height = 4.1,
           bg = "white", device = cairo_pdf)
    ggsave(file.path(OUT, "cav_fig3.png"), fig3, width = 7.2, height = 4.1,
           bg = "white", dpi = 300)
    message("Saved cav_fig3.pdf / cav_fig3.png")
  } else {
    message("Skipping Figure 3: none of its source CSVs were found")
  }

} else {
  message("Skipping Figure 2/3: no temporal_tool_comparison_*.csv files found")
}

# ===========================================================================
# EC figures
# ===========================================================================

# ---------------------------------------------------------------------------
# EC histogram: per-EC-term recall distribution, one curve per tool
# (restricted to tools with full overlapping coverage with CAV ECs)
# ---------------------------------------------------------------------------
ec_per_term_path <- file.path(DATA, "ec_per_term_recall.csv")

p_ec_hist <- NULL
if (file.exists(ec_per_term_path)) {
  ec_per_term <- read_csv(ec_per_term_path, show_col_types = FALSE)

  # Assign Okabe-Ito colors: CAV gets the blue slot, others get remaining colors
  tools_ordered <- c("CAV", sort(setdiff(unique(ec_per_term$tool), "CAV")))
  oi_cycle      <- unname(oi[c("blue", "vermillion", "green", "orange",
                                "sky_blue", "pink", "yellow", "black")])
  tool_colors   <- setNames(oi_cycle[seq_along(tools_ordered)], tools_ordered)

  p_ec_hist <- ec_per_term |>
    mutate(tool = factor(tool, levels = tools_ordered)) |>
    ggplot(aes(x = recall, fill = tool, color = tool)) +
    geom_histogram(
      position = "identity", alpha = 0.45,
      bins = 25, boundary = 0
    ) +
    scale_fill_manual(values  = tool_colors) +
    scale_color_manual(values = tool_colors) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    base_theme() +
    labs(
      x     = "Recall per EC term  (fraction of proteins correctly predicted)",
      y     = "EC terms",
      fill  = NULL,
      color = NULL
    )
} else {
  message("Skipping EC histogram: figure_data/ec_per_term_recall.csv not found")
}

# --- EC summary panels ---
ec_summary_path <- file.path(DATA, "ec_tool_comparison_summary.csv")
ec_llr_path     <- file.path(DATA, "ec_recall_vs_llr.csv")

if (file.exists(ec_summary_path)) {
  ec_summary <- read_csv(ec_summary_path, show_col_types = FALSE)

  p_ec_recall <- ec_summary |>
    mutate(
      tool_label = fct_reorder(tool, recall_exact),
      is_cav     = tool == "CAV"
    ) |>
    ggplot(aes(x = recall_exact, y = tool_label, fill = is_cav)) +
    geom_col(width = 0.65) +
    scale_fill_manual(values = c(`FALSE` = unname(oi["sky_blue"]),
                                 `TRUE`  = TOOL_COLOR)) +
    base_theme() +
    labs(x = "Recall (exact match)", y = NULL) +
    theme(legend.position = "none") +
    xlim(0, 1)

  p_ec_coverage <- ec_summary |>
    mutate(is_cav = tool == "CAV") |>
    ggplot(aes(x = coverage, y = recall_exact,
               color = is_cav, size = is_cav, label = tool)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray60", linewidth = 0.7) +
    geom_point(alpha = 0.9) +
    geom_text_repel(size = 8 / .pt, show.legend = FALSE, max.overlaps = 20) +
    scale_color_manual(values = c(`FALSE` = unname(oi["sky_blue"]),
                                  `TRUE`  = TOOL_COLOR),
                       labels = c("Other tools", "CAV")) +
    scale_size_manual(values = c(`FALSE` = 2, `TRUE` = 3.5), guide = "none") +
    coord_equal(xlim = c(0, 1.05), ylim = c(0, 1.05)) +
    base_theme() +
    labs(x = "Coverage", y = "Recall (exact match)", color = NULL)

  if (file.exists(ec_llr_path)) {
    ec_llr <- read_csv(ec_llr_path, show_col_types = FALSE)

    p_ec_llr <- ec_llr |>
      ggplot(aes(x = llr_threshold, y = recall)) +
      geom_line(color = CAV_COLOR, linewidth = 1) +
      geom_vline(xintercept = 0, linetype = "dashed",
                 color = "gray60", linewidth = 0.7) +
      base_theme() +
      labs(x = "LLR threshold", y = "Recall (fraction of val pairs)")
  } else {
    p_ec_llr <- NULL
  }

  ec_summary_panels <- Filter(Negate(is.null),
                               list(p_ec_hist, p_ec_recall, p_ec_coverage, p_ec_llr))
  if (length(ec_summary_panels) > 0) {
    fig_ec <- plot_grid(
      plotlist   = ec_summary_panels,
      nrow       = 1,
      labels     = "AUTO",
      label_size = 8,
      label_fontfamily = FIG_FONT,
      align      = "hv",
      axis       = "tblr"
    )
    # 3.6in/panel (was 3.0): p_ec_coverage is coord_equal, so under
    # align="hv" its fixed-aspect panel is what absorbs any width the
    # 8pt tick labels take -- at 3.0in it collapses to zero width and
    # ggrepel errors with "Viewport has zero dimension(s)".
    ggsave(file.path(OUT, "fig_ec_eval.pdf"), fig_ec,
           width = 3.6 * length(ec_summary_panels), height = 3.2, device = cairo_pdf)
    message("Saved fig_ec_eval.pdf")
  }
}

# ===========================================================================
# Figure 4 — two-column layout, (A) | (B over C). Reorged 2026-09-04 (was a
# single row A=blank / B=UMAP / C=blank / D=VAV / E=Plexin; both reserved
# blank spacer panels removed at the author's request). (A) CAV
# direction-space UMAP (3-D, single azimuth-330deg view), colored by Pfam
# clan. 7688 concept vectors from the 7692-motif CAV library; coordinates
# decoded directly from the embedded Plotly data in the precomputed
# motif_clustering/results/figures/pfam_umap_3d.html
# (figure_data/pfam_clan_tsne/render_pfam_clan_umap_3d.py), not recomputed.
# No title/legend -- labels are added manually afterward. (B/C) motif-
# localization worked examples, reproducing arxiv 2511.21614v1 Figure 2:
# (B) Q9NHV9/VAV_DROME layerwise CAV score profiles (PF00621/PF00130/
# PF00017/PF00018, all 36 layers, freshly-trained --
# figure_data/vav_motif_repro/train_all_layers.log,
# score_vav_all_layers.py), domain names labeled in a fixed right-hand
# column at the end of each highlighted line instead of a legend, matching
# the reference figure's layout; (C) second worked example, Plexin A1
# (8 domain types), layer 26 only, from the 20k-library L25 CAVs
# (score_bigreceptors.py). SET1_SCHPO/PLDZ_DICDI live in the Appendix-A
# motif supplemental (draft_fig_supp_appendixA.R). B and C are assembled as
# one four-row plot_grid (curve/track/curve/track) rather than two nested
# grids so cowplot can align their panel edges -- nested grids are opaque to
# align="v" and the two y-axis labels have different widths.
# Layer-label convention: CAV filename L{k} indexes hidden_states[k]
# directly (k=0 = embedding, k=1..36 = the 36 transformer layers); paper's
# "indexed from 1" counts the embedding as layer 1, so display_layer = k + 1.
# ===========================================================================

fig4_umap_path <- file.path(OUT, "pfam_clan_umap_3d_mosaic.png")
fig4_vav_path  <- file.path(DATA, "vav_motif_repro", "vav_layerwise_scores.json")

if (file.exists(fig4_umap_path) && file.exists(fig4_vav_path)) {
  source("draft_fig_extra_proteins_common.R")
  VAV_DATA <- file.path(DATA, "vav_motif_repro")

  # --- clan UMAP (left column) ---
  # rasterGrob with both width and height left NULL fits the image to the
  # cell at its native aspect ratio (1458x1200) instead of stretching it,
  # so the figure's overall proportions are chosen to keep the left cell
  # near 1.2:1 and avoid dead space around the render.
  umap3d_png <- readPNG(fig4_umap_path)
  p_fig4_umap <- ggdraw() + draw_grob(rasterGrob(umap3d_png, interpolate = TRUE))

  # --- VAV_DROME multilayer panel ---
  vj <- fromJSON(fig4_vav_path, simplifyVector = FALSE)
  vav_seq_len <- vj$seq_len
  # Bright, high-chroma hues in the reference figure's purple/green/red/blue
  # order. The previous set (#B83280/#2E7D32/#C0392B/#1F4E96) read muddy: the
  # hues were dark and desaturated, and it hard-failed both computable color
  # checks on the all-pairs list -- worst CVD separation OKLab dE 4.2 (deutan,
  # SH2 vs C1) against a target of 8, and a worst normal-vision pair of 12.4
  # (SH2 vs RhoGEF) against a floor of 15, i.e. full-color readers could not
  # reliably separate them either. This set scores 8.9 / 21.0 and passes the
  # lightness-band, chroma, CVD, normal-vision and contrast checks (see
  # dataviz scripts/validate_palette.py, --mode light --pairs all). All-pairs
  # rather than adjacent-pairs is the right list here: the four curves overlap
  # freely across the panel, so any two can end up side by side.
  VAV_DOMAIN_COLORS <- c(PF00621 = "#A020F0", PF00130 = "#00A86B",
                        PF00017 = "#D62728", PF00018 = "#1F77B4")
  VAV_DOMAIN_LABELS <- c(PF00621 = "PF00621 (RhoGEF Domain)", PF00130 = "PF00130 (C1 Domain)",
                        PF00017 = "PF00017 (SH2 Domain)",    PF00018 = "PF00018 (SH3 Domain)")
  VAV_DOMAIN_ORDER  <- c("PF00621", "PF00130", "PF00017", "PF00018")
  VAV_HIGHLIGHT_LAYER <- 26L

  vav_curve_rows <- list(); vav_gt_rows <- list()
  for (motif_id in VAV_DOMAIN_ORDER) {
    dom <- vj$domains[[motif_id]]
    if (is.null(dom)) next
    gt <- vj$ground_truth[[motif_id]]
    vav_gt_rows[[motif_id]] <- tibble(domain = motif_id, start = gt$start, end = gt$end, name = gt$name)
    for (layer_str in names(dom$per_layer_curve)) {
      scores <- unlist(dom$per_layer_curve[[layer_str]])
      vav_curve_rows[[paste(motif_id, layer_str)]] <- tibble(
        domain = motif_id, layer = as.integer(layer_str), position = seq_along(scores), score = scores
      )
    }
  }
  vav_curves <- bind_rows(vav_curve_rows) %>% filter(!is.na(score))
  vav_gt_df  <- bind_rows(vav_gt_rows)
  vav_n_layers <- max(vav_curves$layer)

  # The track carries the colored ground-truth segments and the position axis
  # only. Its per-domain text labels are gone as of the 2/3-size version: at
  # 8pt in a track row this short, PF00130 / PF00017 / PF00018 (adjacent at
  # the C-terminus) collide with each other and with the 793 tick no matter
  # how they are staggered. Nothing is lost -- the right-hand label column
  # names all four domains in exactly these colors, and panel C's track is
  # unlabeled for the same reason, so the two now read consistently.
  build_vav_domain_track <- function() {
    ggplot() +
      geom_rect(aes(xmin = 0, xmax = vav_seq_len, ymin = 0, ymax = 1), fill = "#e8c087") +
      geom_rect(data = vav_gt_df, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = domain)) +
      scale_fill_manual(values = VAV_DOMAIN_COLORS, guide = "none") +
      scale_color_manual(values = VAV_DOMAIN_COLORS, guide = "none") +
      scale_x_continuous(limits = c(0, vav_seq_len), expand = c(0, 0), breaks = c(0, vav_seq_len)) +
      coord_cartesian(ylim = c(0, 1), clip = "off") +
      theme_void(base_family = FIG_FONT) +
      theme(axis.text.x = element_text(size = 8),
            plot.margin = margin(t = 0, b = 2, l = 2, r = VAV_LABEL_MARGIN))
  }

  # per-domain light -> dark color ramp across layers (early layers get a
  # light tint of the domain color, late layers get the full-strength
  # color) -- a plain alpha gradient on one fixed hue was too subtle to
  # tell layers apart, this is a much more extreme light-to-dark contrast.
  # The dark end used to blend 45% toward black, which is the other half of
  # why this panel read muddy: with 36 lines per domain the sweep is most of
  # the ink, and the late layers turned into near-neutral maroon/brown that
  # no longer matched the domain's hue or the track segment beneath it. The
  # blend is now 18%, so the ramp stays inside its hue and only the value
  # changes. The highlighted layer always gets the domain's pure/canonical
  # color (not whatever shade the ramp happens to land on at that layer
  # index), bold and dotted so it reads clearly against the sweep.
  vav_layer_ramps <- lapply(VAV_DOMAIN_COLORS, function(base_hex) {
    light_tint <- colorRampPalette(c("white", base_hex))(100)[18]
    dark_shade <- colorRampPalette(c(base_hex, "black"))(100)[18]
    colorRampPalette(c(light_tint, dark_shade))(vav_n_layers)
  })
  vav_curves <- vav_curves %>%
    mutate(hex = ifelse(layer == VAV_HIGHLIGHT_LAYER,
                        VAV_DOMAIN_COLORS[domain],
                        mapply(function(d, l) vav_layer_ramps[[d]][l], domain, layer)))

  # Domain names sit in a fixed column to the right of the panel, at the end
  # of each highlighted-layer line, instead of in a legend -- this is the
  # layout of the reference figure (arxiv 2511.21614v1 Fig 2). geom_text_repel
  # was used here previously but had nowhere to escape to inside a 30pt right
  # margin, so it dropped the labels back into the middle of the plot on top
  # of the curves. Placement is now computed: start from each line's terminal
  # score, sort descending, and push any label that would collide down by a
  # fixed minimum separation. scale_color_identity(guide = "none") keeps the
  # per-layer color ramp without a second color scale (avoids ggnewscale).
  VAV_LABEL_MARGIN <- 84   # pt of right margin reserved for the label column

  vav_end_labels <- vav_curves %>%
    filter(layer == VAV_HIGHLIGHT_LAYER) %>%
    group_by(domain) %>%
    slice_max(position, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(score))

  # Spread the four labels evenly over the y range, top to bottom, in order of
  # each line's terminal score. The earlier scheme started from the terminal
  # score and pushed colliding labels down by a fixed fraction of the range;
  # that fraction was calibrated against a 4in-tall panel at 6pt, and at 8pt
  # in a panel a third the height it no longer cleared two lines of text, so
  # the labels landed on top of each other. Even spacing is height-independent
  # and the ordering still matches the curves at the right edge.
  vav_y_lim   <- range(vav_curves$score, na.rm = TRUE)
  vav_y_inset <- 0.08 * diff(vav_y_lim)
  vav_end_labels$label_y <- seq(vav_y_lim[2] - vav_y_inset,
                                vav_y_lim[1] + vav_y_inset,
                                length.out = nrow(vav_end_labels))
  vav_end_labels$label_text <- sub(" \\(", "\n(", VAV_DOMAIN_LABELS[vav_end_labels$domain])

  p_vav <- ggplot() +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_line(data = vav_curves %>% filter(layer != VAV_HIGHLIGHT_LAYER),
              aes(position, score, group = interaction(domain, layer), color = hex),
              linewidth = 0.3) +
    geom_line(data = vav_curves %>% filter(layer == VAV_HIGHLIGHT_LAYER),
              aes(position, score, group = domain, color = hex),
              linewidth = 1.1, linetype = "dotted") +
    geom_text(data = vav_end_labels,
              aes(x = vav_seq_len * 1.03, y = label_y, label = label_text, color = hex),
              hjust = 0, vjust = 0.5, lineheight = 0.95,
              size = 8 / .pt, fontface = "bold", show.legend = FALSE) +
    scale_color_identity(guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, vav_seq_len), clip = "off") +
    base_theme() +
    labs(x = NULL, y = "CAV Score",
         title = "Q9NHV9 - VAV_DROME") +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.title = element_text(size = 8, face = "bold"),
          plot.margin = margin(t = 4, r = VAV_LABEL_MARGIN, b = 2, l = 5.5))

  # --- Plexin A1 (second worked example) ---
  # Composite encoding rather than eight unrelated hues: Plexin A1 carries
  # three domains from one family -- PF17960 (TIG), PF18020 (TIG, plexin-
  # specific) and PF01833 (IPT/TIG) -- and the default flat palette scattered
  # them across unrelated colors, hiding the relationship. They now share a
  # single blue hue as a light->mid->dark ramp, ordered by first sequence
  # position (569 / 714 / 864) to match the color-assignment convention used
  # everywhere else in this panel; the other five domains keep distinct hues.
  #
  # Validated with dataviz scripts/validate_palette.py in two parts, since the
  # ramp and the categorical slots answer to different checks:
  #   the five categorical slots + the ramp's mid step, --pairs all:
  #     ALL PASS (worst CVD dE 8.8 deutan, worst normal-vision pair 15.6,
  #     all inside the L band, all chroma >= 0.10, all contrast >= 3:1)
  #   the ramp itself, --ordinal: ALL PASS (monotone L, adjacent dL >= 0.06,
  #     light end 2.88:1 vs surface, hue spread 10deg)
  # The ramp's end steps sit outside the categorical lightness band and its
  # light step below the categorical chroma floor by design -- that is what a
  # ramp is, and --ordinal is the applicable gate for those three slots.
  PLXNA1_TIG_RAMP <- c(PF17960 = "#5B9BD5", PF18020 = "#1F6FB5", PF01833 = "#0B2E5C")
  PLXNA1_COLORS <- c(PF01403 = "#A81E1E", PF01437 = "#A020F0", PF24479 = "#C2298A",
                     PLXNA1_TIG_RAMP,
                     PF08337 = "#E8710A", PF20170 = "#00A86B")
  # legend order puts the three TIG-family slots in one contiguous block so the
  # shared hue reads as a family; the rest follow first sequence position.
  PLXNA1_ORDER <- c("PF01403", "PF01437", "PF24479",
                    "PF17960", "PF18020", "PF01833",
                    "PF08337", "PF20170")

  bigrec_scores <- fromJSON(file.path(VAV_DATA, "extra_proteins", "candidates_bigreceptors_scores.json"), simplifyVector = FALSE)
  # y label shortened to match panel B's, with the layer moved into the title.
  # At the 2/3 canvas the long "CAV Score (Layer 26)" reached the top-left
  # corner of its cell and the "C" panel label landed on its closing paren;
  # the short form also lines the two panels' y titles up with each other.
  plxna1_parts <- build_extra_protein_panel("Q9UIW2", bigrec_scores[["Q9UIW2"]],
                                            title = "Q9UIW2 - Plexin A1 (Layer 26)",
                                            y_lab = "CAV Score",
                                            base_theme_fn = base_theme, parts = TRUE,
                                            domain_colors = PLXNA1_COLORS,
                                            domain_order = PLXNA1_ORDER,
                                            # 8pt everywhere it fits, but not here: at a 6.33in
                                            # canvas an 8pt "PF20170 (RhoGTPase-binding domain)"
                                            # is ~2.3in of legend, more than a third of the whole
                                            # figure, and it starved the plot to about an inch.
                                            # 6pt keeps this legend near 1.6in. The track's
                                            # position axis has short labels and does fit at 8.
                                            legend_text_size = 6, axis_text_size = 8)

  # --- Right column: B (VAV curve + track) over C (Plexin curve + track) ---
  # Built as one four-row grid so align = "v" / axis = "lr" lines up all four
  # panels' left and right edges; B and C carry y-axis labels of different
  # widths ("CAV Score" vs "CAV Score (Layer 26)") and would otherwise sit at
  # different x offsets. Labels go on the two curve rows only.
  fig4_right <- plot_grid(
    p_vav, build_vav_domain_track(), plxna1_parts$curve, plxna1_parts$track,
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1, 0.16, 1, 0.16),
    labels = c("B", "", "C", ""), label_size = 8, label_fontfamily = FIG_FONT
  )

  fig4 <- plot_grid(p_fig4_umap, fig4_right, nrow = 1,
                    labels = c("A", ""), label_size = 8,
                    label_fontfamily = FIG_FONT,
                    # Re-weighted for the 2/3 canvas: point sizes are absolute, so
                    # shrinking the figure while holding text at 8pt makes the text a
                    # far larger share of the width. The right column carries all of it
                    # (two y-axis labels, panel B's label column, panel C's legend), so
                    # it now takes the larger share and panel A gives width back.
                    rel_widths = c(1, 1.35))

  # 6.33 x 2.67in -- two thirds of the previous 9.5 x 4.0, at the author's
  # request. Font sizes are deliberately NOT scaled with it: ggplot point
  # sizes are absolute, so saving at the size the figure will actually be
  # printed keeps 8pt text as true 8pt on the page. (Saving large and letting
  # the manuscript scale the graphic down is what turns 8pt into ~5pt.)
  # Everything in this figure is 8pt except panel C's legend -- see the note
  # at the build_extra_protein_panel call above.
  ggsave(file.path(OUT, "cav_fig4.pdf"), fig4, width = 6.33, height = 2.67, bg = "white", device = cairo_pdf)
  ggsave(file.path(OUT, "cav_fig4.png"), fig4, width = 6.33, height = 2.67, bg = "white", dpi = 600)
  message("Saved cav_fig4.pdf / cav_fig4.png")
} else {
  message("Skipping Figure 4: ", fig4_umap_path, " and/or ", fig4_vav_path, " not found ",
          "(run figure_data/pfam_clan_tsne/render_pfam_clan_umap_3d.py and ",
          "figure_data/vav_motif_repro/score_vav_all_layers.py first)")
}

# ===========================================================================
# Figure 6 — single-cell case studies (DE-vs-CAV scatter + transcriptional
# continuum, one column per case). Data comes from
# single_cell/scripts/export_fig5_case_studies.py (run via
# single_cell/08_export_fig5_data.sh), which reuses cav_continuum_viz.py's
# cell/gene selection, restricted to donors contributing both conditions.
# The earlier exploratory matplotlib versions (export_fig5_data.py,
# 07_paper_cases.sh) are in single_cell/archive/ -- see its README.
# ===========================================================================

sc_meta_path <- file.path(DATA, "sc_pair_meta.csv")

if (file.exists(sc_meta_path)) {

  sc_meta       <- read_csv(sc_meta_path, show_col_types = FALSE)
  sc_cont_cells <- read_csv(file.path(DATA, "sc_continuum_cells.csv"), show_col_types = FALSE)
  sc_cont_genes <- read_csv(file.path(DATA, "sc_continuum_genes.csv"), show_col_types = FALSE)
  sc_de_vs_cav  <- read_csv(file.path(DATA, "sc_de_vs_cav.csv"),       show_col_types = FALSE)

  # Skin/melanoma and lung pairs were dropped: every skin_epidermis/melanoma
  # pair in this atlas is 100% assay-confounded (normal=10x 3' v2,
  # melanoma=Smart-seq2, no overlap). The original lung pairs carried a
  # milder version of the same confound (normal ~91% 10x, cancer ~50% 10x);
  # epithelial_cell__lung10x__* is a from-scratch retrain restricted to
  # assay=="10x 3' v2" cells only (see single_cell/scripts/
  # train_10x_lung_cavs.py), so it's clean like the colorectum pairs.
  SC_PAIR_ORDER <- c(
    "neutrophil__breast__normal_vs_breast_cancer",
    "epithelial_cell__lung10x__normal_vs_lung_cancer",
    "fibroblast__colorectum__normal_vs_colorectal_cancer"
  )
  SC_PAIR_TITLES <- c(
    "neutrophil__breast__normal_vs_breast_cancer"               = "Neutrophil — breast",
    # The assay restriction (10x 3' v2 only) is stated in the Methods rather
    # than in the panel title.
    "epithelial_cell__lung10x__normal_vs_lung_cancer"           = "Epithelial cell — lung",
    "fibroblast__colorectum__normal_vs_colorectal_cancer"       = "Fibroblast — colorectum"
  )

  SC_NORMAL_COLOR     <- "#3a6fad"
  SC_CANCER_COLOR     <- "#c0392b"
  SC_BACKGROUND_COLOR <- "#2277bb"
  SC_CONTINUUM_FILL   <- "#f0a500"
  SC_CONTINUUM_EDGE   <- "#7a5300"

  # ---------------------------------------------------------------------
  # Overview row: the full top bar is left blank for panel A (a schematic
  # dropped in externally). The three hierarchy-level t-SNEs that used to
  # occupy B/C/D here were removed -- they showed near-complete overlap of
  # normal/cancer at every level, so they carried no information the case
  # columns don't carry better.
  # ---------------------------------------------------------------------
  sc_overview_row <- plot_grid(NULL, nrow = 1,
                               labels = c("A"), label_size = 8,
                               label_fontfamily = FIG_FONT)

  # ---------------------------------------------------------------------
  # DE-vs-CAV scatter, one per case. Every gene is plotted; the reported
  # Spearman r is computed over all of them. The |log2FC| >= 1.5 threshold
  # only controls appearance -- genes below it are drawn as a faint grey
  # background cloud, and labels are drawn for the "meaningful" subset
  # (passes the DE magnitude threshold, OR is one of that case's top
  # continuum genes).
  # ---------------------------------------------------------------------
  sc_build_scatter <- function(pair_id) {
    df <- sc_de_vs_cav %>% filter(pair == pair_id)
    continuum_ids <- sc_cont_genes %>% filter(pair == pair_id) %>% pull(gene) %>% unique()

    df <- df %>%
      mutate(
        is_continuum   = gene %in% continuum_ids,
        above_threshold = abs(log2fc) >= 1.5,
        category = case_when(
          is_continuum              ~ "continuum",
          above_threshold           ~ "background",
          TRUE                      ~ "below_threshold"
        )
      )

    # Panel statistic: Spearman over EVERY gene that has both a DE estimate and
    # a CAV correlation (the inner join written by export_fig5_case_studies.py).
    # (An earlier draft annotated rho over the thresholded subset only; for the
    # colorectum pair that subset was the 10 continuum genes alone -- an n = 10,
    # p = 0.10 statistic reported without its n.)
    df_all <- df %>% filter(is.finite(log2fc), is.finite(cav_r))
    rho    <- suppressWarnings(cor(df_all$log2fc, df_all$cav_r, method = "spearman"))
    n_rho  <- nrow(df_all)

    df_stat <- df %>% filter(category != "below_threshold")   # labelling only

    # Label every continuum gene, plus the top-8 by distance from the
    # origin among the meaningful subset (captures strong agreement in
    # either direction without picking up threshold noise).
    nx <- df_stat$log2fc / (max(abs(df_stat$log2fc)) + 1e-9)
    ny <- df_stat$cav_r   / (max(abs(df_stat$cav_r))   + 1e-9)
    df_stat$dist <- sqrt(nx^2 + ny^2)
    auto_label <- df_stat %>% slice_max(dist, n = 8) %>% pull(gene)
    label_ids  <- union(auto_label, df_stat %>% filter(category != "background") %>% pull(gene))
    df$show_label <- df$gene %in% label_ids

    ggplot(df, aes(log2fc, cav_r)) +
      geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
      geom_point_rast(data = filter(df, category == "below_threshold"),
                 color = "grey70", size = 0.4, alpha = 0.15, raster.dpi = RASTER_DPI) +
      geom_point_rast(data = filter(df, category == "background"),
                 color = SC_BACKGROUND_COLOR, size = 0.7, alpha = 0.4, raster.dpi = RASTER_DPI) +
      geom_point(data = filter(df, category == "continuum"),
                 shape = 21, fill = SC_CONTINUUM_FILL, color = SC_CONTINUUM_EDGE,
                 size = 1.6, stroke = 0.3, alpha = 0.9) +
      geom_text_repel(
        data = filter(df, show_label),
        aes(label = gene_name, color = category),
        size = 2.0, max.overlaps = 30, segment.size = 0.25,
        segment.color = "grey60", show.legend = FALSE,
        box.padding = 0.35, point.padding = 0.15, force = 2, force_pull = 0.5,
        min.segment.length = 0.1, seed = 42
      ) +
      scale_color_manual(values = c(background = "grey30",
                                    continuum  = SC_CONTINUUM_EDGE)) +
      base_theme() +
      labs(x = "log2FC (DE)", y = "Pearson r (CAV)",
           title = SC_PAIR_TITLES[pair_id],
           subtitle = sprintf("Spearman r = %.2f (n = %s genes)", rho,
                              format(n_rho, big.mark = ","))) +
      theme(plot.title    = element_text(size = 8, face = "bold"),
            plot.subtitle = element_text(size = 7))
  }

  # ---------------------------------------------------------------------
  # Continuum panel: top strip (cells on the L2 axis) + per-gene
  # expression strips below, each split into a normal sub-band (top
  # half) and cancer sub-band (bottom half) rather than jittering both
  # conditions together -- normal and cancer populations overlap
  # substantially along this axis in these pairs (see e.g. the
  # ADAMDEC1/fibroblast/colorectum case: within-condition L2 correlation
  # is actually stronger for normal cells (r=0.28) than cancer cells
  # (r=0.14), and normal cells in the high-L2 zone show higher expression
  # than cancer cells there), so a single interleaved jitter band
  # visually overstates the separation. Point colour is by condition
  # (cancer = warm YlOrRd ramp, normal = cool YlGnBu ramp; both from
  # cav_continuum_viz.py's 0.2-0.95 ramp), with intensity within each
  # ramp encoding expr_scaled -- via scale_color_identity() so the two
  # colour families coexist without a second colour scale. Each cell
  # gets one jitter value (keyed by within-pair row order, shared
  # between sc_continuum_cells.csv and sc_continuum_genes.csv since both
  # are written from the same per-cell array in export_fig5_case_studies.py)
  # reused across the strip and every gene row, so a given cell sits at
  # the same relative height throughout the column.
  # ---------------------------------------------------------------------
  # One ramp for both conditions -- condition is already encoded by the
  # normal/cancer sub-band position (and by the blue/red cell strip above),
  # so the gene rows use a single light-grey-to-black intensity ramp for
  # expression instead of two competing colour families.
  # Floor is a visible light grey, not near-white: ~62% of expr_scaled values
  # are exactly 0 (dropout), and those points still have to read as band
  # structure at print size.
  sc_ramp_expr <- colorRampPalette(c("#d2d2d2", "#000000"))(101)

  sc_expr_to_hex <- function(expr_scaled, is_cancer) {
    idx <- pmin(pmax(round(expr_scaled * 100) + 1, 1), 101)
    sc_ramp_expr[idx]
  }

  sc_band_ticks <- function(row_y, x_tick) {
    tibble(x   = x_tick,
           y   = c(row_y + 0.22, row_y - 0.22),
           lab = rep(c("N", "C"), each = length(row_y)))
  }

  sc_build_continuum <- function(pair_id, seed = 1, brackets = FALSE) {
    cells <- sc_cont_cells %>% filter(pair == pair_id) %>% mutate(row_id = row_number())

    set.seed(seed)
    cells <- cells %>%
      mutate(y_base   = if_else(is_baseline, 0.5, -0.5),
             jit_cell = runif(n(), -0.16, 0.16),
             y        = y_base + jit_cell)

    genes <- sc_cont_genes %>% filter(pair == pair_id) %>%
      group_by(gene) %>% mutate(row_id = row_number()) %>% ungroup() %>%
      left_join(cells %>% select(row_id, is_baseline, jit_cell), by = "row_id")

    x_lo <- min(c(cells$l2_score, genes$l2_score))
    x_hi <- max(c(cells$l2_score, genes$l2_score))
    pad  <- 0.05 * (x_hi - x_lo)
    # Extra empty gutter on the left of every panel in the column (strip
    # included, so the x scales stay identical and the panels still align)
    # to hold the N/C sub-band ticks. The top-5 / bottom-5 brackets live
    # further out still, outside the panel and left of the gene names --
    # see sc_build_bracket_col() -- and only on the first column.
    gutter   <- 0.07 * (x_hi - x_lo)
    xlim     <- c(x_lo - pad - gutter, x_hi + pad)
    x_tick   <- x_lo - pad - gutter * 0.45

    p_strip <- ggplot(cells, aes(l2_score, y, color = is_baseline)) +
      geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
      geom_point_rast(size = 0.4, alpha = 0.54, raster.dpi = RASTER_DPI) +
      # Same N/C ticks as the gene rows below, in the same gutter, so the
      # sub-band convention is stated once at the top of the column too.
      geom_text(data = tibble(x = x_tick, y = c(0.5, -0.5), lab = c("N", "C")),
                aes(x, y, label = lab), inherit.aes = FALSE,
                size = 1.7, color = "grey35", family = FIG_FONT) +
      scale_color_manual(values = c(`TRUE` = SC_NORMAL_COLOR, `FALSE` = SC_CANCER_COLOR)) +
      coord_cartesian(xlim = xlim, ylim = c(-1, 1)) +
      theme_void() +
      theme(legend.position = "none",
            plot.margin = margin(t = 2, b = 0, l = 2, r = 2))

    gene_order <- genes %>% distinct(gene_name, rank) %>% arrange(rank) %>% pull(gene_name)
    n_genes <- length(gene_order)
    genes <- genes %>%
      mutate(gene_name = factor(gene_name, levels = rev(gene_order)),
             y_ctr      = n_genes - rank)

    # The two gene blocks are the top 5 and bottom 5 genes ranked purely by
    # correlation with the L2 score (not by DE). They're pushed apart by GAP
    # and bracketed in the gutter -- an outer bracket spanning both says what
    # the ranking is on, so the split can't be read as a DE call.
    GAP <- 0.8
    genes <- genes %>% mutate(r_up = r > 0,
                              y_ctr = y_ctr + if_else(r_up, GAP, 0))
    y_top <- n_genes + 0.5 + GAP

    genes$y_off <- if_else(genes$is_baseline, 0.22, -0.22)
    genes$hex <- sc_expr_to_hex(genes$expr_scaled, !genes$is_baseline)

    gene_labels <- genes %>% distinct(y_ctr, gene_name, r, r_up) %>%
      mutate(lab = sprintf("%s (r=%.2f)", gene_name, r)) %>% arrange(y_ctr)

    blocks <- gene_labels %>% group_by(r_up) %>%
      summarise(ymin = min(y_ctr) - 0.5, ymax = max(y_ctr) + 0.5, .groups = "drop") %>%
      # Deliberately neutral: blue/red are already spoken for by the cell
      # strip (= a cell's condition), and colouring these blocks with the
      # same swatches would read as "normal genes / cancer genes" rather
      # than "ranked by correlation with the L2 score".
      mutate(color = "grey25",
             lab   = if_else(r_up, "top 5", "bottom 5"),
             ymid  = (ymin + ymax) / 2)

    outer_bracket <- tibble(ymin = min(blocks$ymin), ymax = max(blocks$ymax)) %>%
      mutate(ymid = (ymin + ymax) / 2,
             lab  = "Correlation with cancer\ndirection on L2 axis")

    p_genes <- ggplot(genes, aes(l2_score, y_ctr + y_off + jit_cell)) +
      geom_hline(yintercept = sort(unique(c(gene_labels$y_ctr - 0.5, gene_labels$y_ctr + 0.5))),
                 color = "grey90", linewidth = 0.3) +
      geom_hline(yintercept = gene_labels$y_ctr, color = "grey95", linewidth = 0.2) +
      {if (x_lo < 0 && x_hi > 0) geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4, linetype = "dashed")} +
      geom_point_rast(aes(color = hex), size = 0.35, alpha = 0.54, raster.dpi = RASTER_DPI) +
      # N/C ticks: with expression now on a single grey ramp, nothing inside
      # a gene row says which sub-band is which, so mark them explicitly.
      geom_text(data = sc_band_ticks(gene_labels$y_ctr, x_tick),
                aes(x, y, label = lab), inherit.aes = FALSE,
                size = 1.7, color = "grey35", family = FIG_FONT) +
      scale_color_identity() +
      scale_y_continuous(breaks = gene_labels$y_ctr, labels = gene_labels$lab,
                          limits = c(0.5, y_top), expand = c(0, 0)) +
      coord_cartesian(xlim = xlim) +
      base_theme() +
      labs(x = "L2 score  (← normal    cancer →)", y = NULL) +
      theme(axis.text.y  = element_text(size = 6),
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank(),
            plot.margin  = margin(t = 0, b = 2, l = 2, r = 2))

    inner <- plot_grid(p_strip, p_genes, ncol = 1, align = "v", axis = "lr",
                       rel_heights = c(0.9, n_genes + GAP))

    if (!brackets) return(inner)

    # Bracket column: its own plot so the labels sit outside the panel, to
    # the left of the gene names. It carries an invisible copy of p_genes'
    # x axis so the two panels end up the same height, and is stacked under
    # a spacer matching the cell strip's share of the column.
    p_ann <- ggplot() +
      geom_segment(data = blocks, aes(x = 0.70, xend = 0.70,
                                      y = ymin + 0.08, yend = ymax - 0.08),
                   color = "grey25", linewidth = 0.7, lineend = "round") +
      geom_text(data = blocks, aes(x = 0.88, y = ymid, label = lab),
                color = "grey25", angle = 90, size = 2.4, family = FIG_FONT) +
      geom_segment(data = outer_bracket, aes(x = 0.48, xend = 0.48,
                                             y = ymin, yend = ymax),
                   color = "grey45", linewidth = 0.5) +
      geom_text(data = outer_bracket, aes(x = 0.26, y = ymid, label = lab),
                color = "grey25", angle = 90, size = 2.4, family = FIG_FONT,
                lineheight = 0.9) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.5, y_top), expand = c(0, 0)) +
      base_theme() +
      labs(x = "L2 score", y = NULL) +
      theme(axis.title.x = element_text(color = "white"),
            axis.text.x  = element_text(color = "white"),
            axis.ticks.x = element_line(color = "white"),
            axis.line    = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin  = margin(t = 0, b = 2, l = 2, r = 0))

    left_col <- plot_grid(NULL, p_ann, ncol = 1, rel_heights = c(0.9, n_genes + GAP))
    plot_grid(left_col, inner, nrow = 1, rel_widths = c(0.20, 1))
  }

  sc_columns <- map2(SC_PAIR_ORDER, seq_along(SC_PAIR_ORDER), function(pid, i) {
    plot_grid(sc_build_scatter(pid), sc_build_continuum(pid, brackets = (i == 1)),
              ncol = 1, align = "v", axis = "lr",
              rel_heights = c(0.62, 1))
  })

  case_labels <- c("B", "C", "D")
  fig6_body <- plot_grid(plotlist = sc_columns, nrow = 1,
                         labels = case_labels, label_size = 8, label_fontfamily = FIG_FONT)

  sc_legend_strip_df <- tibble(x = 1, y = 1,
                                cond = factor(c("normal", "cancer"), levels = c("normal", "cancer")))
  sc_legend_strip <- get_legend(
    ggplot(sc_legend_strip_df, aes(x, y, color = cond)) +
      geom_point(size = 2.5) +
      scale_color_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR),
                         name = "Cell strip") +
      base_theme() + theme(legend.position = "right", legend.key.size = unit(9, "pt"))
  )

  sc_legend_scatter_df <- tibble(
    x = 1:2, y = 1,
    cat = factor(c("other gene (|log2FC| >= threshold)",
                   "top continuum gene"),
                 levels = c("other gene (|log2FC| >= threshold)",
                            "top continuum gene"))
  )
  # Keys drawn as shape 21 with the same fill/outline the scatter points use,
  # so the continuum key reads as the points' amber fill rather than their
  # dark outline colour.
  sc_legend_scatter <- get_legend(
    ggplot(sc_legend_scatter_df, aes(x, y, fill = cat, color = cat)) +
      geom_point(size = 2.5, shape = 21, stroke = 0.3) +
      scale_fill_manual(values = setNames(c(SC_BACKGROUND_COLOR, SC_CONTINUUM_FILL),
                                          levels(sc_legend_scatter_df$cat)),
                        name = NULL) +
      scale_color_manual(values = setNames(c(SC_BACKGROUND_COLOR, SC_CONTINUUM_EDGE),
                                           levels(sc_legend_scatter_df$cat)),
                         name = NULL) +
      base_theme() + theme(legend.position = "right", legend.key.size = unit(9, "pt"),
                            legend.text = element_text(size = 7))
  )

  # Colourbar key for the grey expression ramp used in the gene rows, plus a
  # reminder of what the N/C sub-band ticks mean.
  sc_legend_expr <- get_legend(
    ggplot(tibble(x = 1:2, y = 1, e = c(0, 1)), aes(x, y, color = e)) +
      geom_point(size = 2.5) +
      scale_color_gradient(low = "#d2d2d2", high = "#000000",
                           name = "Scaled expression",
                           breaks = c(0, 1), labels = c("low", "high"),
                           guide = guide_colourbar(title.position = "top",
                                                   barwidth = unit(30, "pt"),
                                                   barheight = unit(5, "pt"),
                                                   direction = "horizontal")) +
      base_theme() + theme(legend.position = "right",
                            legend.title = element_text(size = 6),
                            legend.text  = element_text(size = 6))
  )

  sc_legend_row <- plot_grid(NULL, sc_legend_scatter, sc_legend_strip, sc_legend_expr, NULL,
                             nrow = 1, rel_widths = c(0.2, 1.4, 0.8, 1.2, 0.2))

  fig6_cases <- plot_grid(fig6_body, sc_legend_row, ncol = 1, rel_heights = c(1, 0.13))

  # 7.5in ~ 190mm, a standard full double-column journal figure width;
  # heights rescaled (x 7.5/12) to preserve the original aspect ratios.
  # Panel A's blank bar is 30% shorter than the old t-SNE row (0.6 -> 0.42),
  # with the overall height rescaled to match so the case columns keep their
  # original aspect ratio.
  fig6 <- plot_grid(sc_overview_row, fig6_cases, ncol = 1, rel_heights = c(0.42, 1))
  fig6_height <- 5.82

  fig6_width <- 7.5
  ggsave(file.path(OUT, "cav_fig6.pdf"), fig6, width = fig6_width, height = fig6_height,
         bg = "white", device = cairo_pdf)
  ggsave(file.path(OUT, "cav_fig6.png"), fig6, width = fig6_width, height = fig6_height,
         bg = "white", dpi = 300)
  message("Saved cav_fig6.pdf / cav_fig6.png")

  # -------------------------------------------------------------------
  # Candidate supplemental — does orthogonalizing the condition CAV
  # against cell-type/tissue baseline structure change the DE-vs-CAV
  # relationship, compared to using the raw (unprocessed) condition CAV?
  # Three rows = three stages of the SAME subtraction chain that produces
  # L2: raw -> minus cell-type baseline (L0) -> minus cell-type + tissue
  # (L0+L1 = L2, identical to the L2 axis used everywhere else in Figure
  # 6). There is no fourth stage to subtract in this two-level hierarchy.
  # Same three case-study pairs and paired-donor cell populations as
  # panels E-G. Data from single_cell/scripts/
  # export_projection_comparison_supp_data.py.
  # -------------------------------------------------------------------
  proj_path <- file.path(DATA, "proj_de_vs_cav.csv")

  if (file.exists(proj_path)) {

    proj_df <- read_csv(proj_path, show_col_types = FALSE)

    PROJ_LEVEL_ORDER <- c("raw", "L0_removed", "L2")
    PROJ_LEVEL_TITLES <- c(
      raw        = "Raw condition CAV\n(no orthogonalization)",
      L0_removed = "Minus cell-type baseline (L0)",
      L2         = "Minus cell-type + tissue (L0+L1 = L2)"
    )

    proj_build_scatter <- function(pair_id, level_id) {
      df <- proj_df %>% filter(pair == pair_id, level == level_id) %>%
        filter(is.finite(log2fc), is.finite(cav_r))

      # Same convention as Figure 6: rho over all genes. No gene labels and no
      # above-threshold recolouring here -- neither is referred to in the text,
      # and this figure exists for the correlations alone.
      rho   <- suppressWarnings(cor(df$log2fc, df$cav_r, method = "spearman"))
      n_rho <- nrow(df)

      ggplot(df, aes(log2fc, cav_r)) +
        geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_point_rast(color = "grey65", size = 0.4, alpha = 0.2,
                        raster.dpi = RASTER_DPI) +
        base_theme() +
        labs(x = "log2FC (DE)", y = "Pearson r (CAV)",
             subtitle = sprintf("Spearman r = %.2f (n = %s genes)", rho,
                                format(n_rho, big.mark = ","))) +
        theme(axis.text     = element_text(size = 6),
              axis.title    = element_text(size = 6.5),
              plot.subtitle = element_text(size = 6.5),
              plot.margin   = margin(t = 2, r = 4, b = 2, l = 2))
    }

    proj_col_headers <- plot_grid(plotlist = map(SC_PAIR_ORDER, function(pid) {
      ggdraw() + draw_label(SC_PAIR_TITLES[pid], fontface = "bold", size = 7, fontfamily = FIG_FONT)
    }), nrow = 1)

    proj_row_label <- function(lvl) {
      ggdraw() + draw_label(PROJ_LEVEL_TITLES[lvl], fontface = "bold", size = 7, angle = 90, fontfamily = FIG_FONT)
    }

    proj_body_rows <- map(seq_along(PROJ_LEVEL_ORDER), function(i) {
      lvl <- PROJ_LEVEL_ORDER[i]
      row_plots <- map(SC_PAIR_ORDER, function(pid) proj_build_scatter(pid, lvl))
      plot_grid(proj_row_label(lvl), plot_grid(plotlist = row_plots, nrow = 1),
                nrow = 1, rel_widths = c(0.07, 1))
    })

    proj_body <- plot_grid(plotlist = proj_body_rows, ncol = 1,
                           labels = LETTERS[1:3], label_size = 8, label_fontfamily = FIG_FONT)

    fig_supp_projection <- plot_grid(
      plot_grid(NULL, proj_col_headers, nrow = 1, rel_widths = c(0.07, 1)),
      proj_body, ncol = 1, rel_heights = c(0.05, 1)
    )

    # Two-thirds of the original 11 x 9.5 in; type trimmed to suit rather than
    # scaled down with the panels.
    proj_w <- 7.35
    proj_h <- 6.35
    ggsave(file.path(OUT, "fig_supp_projection_comparison.pdf"), fig_supp_projection,
           width = proj_w, height = proj_h, bg = "white", device = cairo_pdf)
    ggsave(file.path(OUT, "fig_supp_projection_comparison.png"), fig_supp_projection,
           width = proj_w, height = proj_h, bg = "white", dpi = 300)
    message("Saved fig_supp_projection_comparison.pdf / .png")

  } else {
    message("Skipping projection-comparison supplemental: figure_data/proj_de_vs_cav.csv not found ",
            "(run single_cell/10_export_projection_comparison_supp_data.sh first)")
  }

} else {
  message("Skipping Figure 6: figure_data/sc_pair_meta.csv not found ",
          "(run single_cell/08_export_fig5_data.sh first)")
}

# ===========================================================================
# Supplemental — CAV detects subpopulation-abundance shifts that
# population-average DE misses. Worked example: fibroblast/colorectum's
# quiescent-fibroblast module (MGP, DCN, OGN, C3, CCDC80) is DE-null
# (padj > 0.5 for all five) but a strong, highly significant CAV correlate
# (r=-0.33 to -0.37, p ~1e-53 to 1e-65) -- consistent with a subpopulation
# whose relative abundance shifts with disease, rather than uniform
# per-cell downregulation, which a group-mean DE test is powered to detect
# but a compositional shift is not. Data from single_cell/scripts/
# export_subpop_supp_data.py.
# ===========================================================================

subpop_scatter_path <- file.path(DATA, "subpop_gene_scatter.csv")

if (file.exists(subpop_scatter_path)) {

  subpop_bg     <- read_csv(file.path(DATA, "subpop_gene_scatter_background.csv"), show_col_types = FALSE)
  subpop_genes  <- read_csv(subpop_scatter_path, show_col_types = FALSE) |> rename(cav_r = r)
  subpop_corr   <- read_csv(file.path(DATA, "subpop_gene_corr.csv"), show_col_types = FALSE)
  subpop_cells  <- read_csv(file.path(DATA, "subpop_cell_scores.csv"), show_col_types = FALSE)

  subpop_theme <- function() {
    base_theme() +
      theme(axis.text    = element_text(size = 6),
            axis.title   = element_text(size = 6.5),
            legend.text  = element_text(size = 5.5),
            legend.title = element_text(size = 6),
            plot.title   = element_text(size = 7, face = "bold"),
            legend.key.size = unit(6, "pt"),
            legend.margin   = margin(t = 0, b = 0),
            plot.margin     = margin(t = 2, r = 3, b = 2, l = 2))
  }

  MODULE_COLORS <- c(quiescent = "#3a6fad", activated = "#c0392b", independent = "#7a5300")
  MODULE_LABELS <- c(quiescent = "quiescent-fibroblast module (MGP/DCN/OGN/C3/CCDC80)",
                     activated = "CXCR4 (anti-correlated with the module)",
                     independent = "ADAMDEC1 / CXCL14 (independent of the module)")

  # --- Panel A: DE log2FC vs CAV r, background cloud + highlighted genes ---
  p_subpop_a <- ggplot(subpop_bg, aes(log2fc, cav_r)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_point_rast(color = "grey70", size = 0.5, alpha = 0.25, raster.dpi = RASTER_DPI) +
    geom_point(data = subpop_genes, aes(color = module), size = 1.8) +
    ggrepel::geom_text_repel(data = subpop_genes, aes(label = gene_name, color = module),
                             size = 2.0, fontface = "bold", show.legend = FALSE,
                             box.padding = 0.5, point.padding = 0.25,
                             force = 4, force_pull = 0.4,
                             min.segment.length = 0, segment.size = 0.25,
                             segment.color = "grey55", max.overlaps = 20, seed = 42) +
    scale_color_manual(values = MODULE_COLORS, labels = MODULE_LABELS, name = NULL) +
    subpop_theme() +
    labs(x = "log2FC (DE, mixedlm)", y = "Pearson r (CAV)",
         title = "DE-null, CAV-strong genes") +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 3))

  # --- Panel B: co-expression heatmap among the 7 genes ---
  gene_order <- c("MGP", "DCN", "OGN", "C3", "CCDC80", "CXCR4", "ADAMDEC1", "CXCL14")
  p_subpop_b <- subpop_corr |>
    mutate(gene1 = factor(gene1, levels = gene_order),
           gene2 = factor(gene2, levels = rev(gene_order))) |>
    ggplot(aes(gene1, gene2, fill = r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 1.7) +
    scale_fill_gradient2(low = unname(oi["vermillion"]), mid = "white", high = unname(oi["blue"]),
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    subpop_theme() +
    labs(x = NULL, y = NULL, title = "Pairwise co-expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 5.5),
          legend.key.width = unit(5, "pt"), legend.key.height = unit(14, "pt"))

  # --- Panel C: compositional shift -- fraction of cells in each
  # module-defined bucket, by disease. Real but modest (smaller than the
  # raw CAV r's alone might suggest): module-negative cells go from 44%
  # (normal) to 51% (cancer); top-quartile-by-score cells go from 27%
  # (normal) to 23% (cancer).
  composition_df <- subpop_cells |>
    mutate(disease = if_else(str_detect(tolower(disease), "normal"), "normal", "cancer"),
           disease = factor(disease, levels = c("normal", "cancer"))) |>
    group_by(disease) |>
    summarise(
      `Module-negative\n(all 5 = 0)`      = mean(module_negative),
      `High module score\n(top quartile)` = mean(high_quiescent),
      .groups = "drop"
    ) |>
    pivot_longer(-disease, names_to = "metric", values_to = "frac")

  p_subpop_c <- composition_df |>
    ggplot(aes(metric, frac, fill = disease)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = scales::percent(frac, accuracy = 1)),
              position = position_dodge(width = 0.7), vjust = -0.4, size = 2.0) +
    scale_fill_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR), name = NULL) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
    subpop_theme() +
    labs(x = NULL, y = "Share of fibroblasts") +
    theme(legend.position = "top",
          axis.title.y = element_text(size = 6.5, margin = margin(r = 4)),
          plot.margin  = margin(t = 2, r = 3, b = 2, l = 5))

  # --- Panel D: quiescent-module score vs. the continuous CAV L2 axis ---
  p_subpop_d <- subpop_cells |>
    mutate(disease = if_else(str_detect(tolower(disease), "normal"), "normal", "cancer"),
           disease = factor(disease, levels = c("normal", "cancer"))) |>
    ggplot(aes(l2_score, quiescent_score, color = disease)) +
    geom_point_rast(size = 0.6, alpha = 0.4, raster.dpi = RASTER_DPI) +
    geom_smooth(aes(group = 1), method = "loess", color = "black", linewidth = 0.6, se = TRUE) +
    scale_color_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR), name = NULL) +
    subpop_theme() +
    labs(x = "L2 score  (← normal    cancer →)", y = "Quiescent-module score") +
    theme(legend.position = "top")

  fig_subpop <- plot_grid(p_subpop_a, p_subpop_b, p_subpop_c, p_subpop_d,
                          nrow = 2, labels = "AUTO", label_size = 8, label_fontfamily = FIG_FONT,
                          align = "hv", axis = "tblr")

  # Two-thirds of the original 8.5 x 8 in; font sizes are left as-is so the
  # text stays legible at print size rather than shrinking with the panels.
  subpop_w <- 5.7
  subpop_h <- 5.35
  ggsave(file.path(OUT, "fig_supp_fibroblast_subpop.pdf"), fig_subpop,
         width = subpop_w, height = subpop_h, bg = "white", device = cairo_pdf)
  ggsave(file.path(OUT, "fig_supp_fibroblast_subpop.png"), fig_subpop,
         width = subpop_w, height = subpop_h, bg = "white", dpi = 300)
  message("Saved fig_supp_fibroblast_subpop.pdf / .png")

} else {
  message("Skipping fibroblast-subpopulation supplemental: figure_data/subpop_gene_scatter.csv not found ",
          "(run single_cell/scripts/export_subpop_supp_data.py first)")
}

# ===========================================================================
# Figure 5 — SERK vs. CIK LRR-RK co-receptor specificity
# ===========================================================================
# Data source: interface_search/fastas/pairs/lrr_serkcik/, exported to
# figure_data/fig6_heatmap.csv, fig6_densities.csv, fig6_auroc.csv by the
# python analysis in that directory (see SESSION_FINDINGS.md there).
# Layout: 3 rows -- [A B] placeholders / [C (0.6) placeholder, D (0.4)
# heatmap] / [E, full width, two densities]. A/B/C are Illustrator or
# not-yet-built panels; D/E are built here.

fig6_heat_path  <- file.path(DATA, "fig6_heatmap.csv")
fig6_dens_path  <- file.path(DATA, "fig6_densities.csv")
fig6_auroc_path <- file.path(DATA, "fig6_auroc.csv")

if (file.exists(fig6_heat_path) && file.exists(fig6_dens_path) && file.exists(fig6_auroc_path)) {

  SERK_COLOR <- CAV_COLOR
  CIK_COLOR  <- unname(oi["vermillion"])

  placeholder_panel <- function(label_text) {
    ggplot() +
      annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.02, ymax = 0.98,
               fill = "#fafaf9", color = "#b9b8b3", linewidth = 0.4, linetype = "42") +
      annotate("text", x = 0.5, y = 0.5, label = label_text, size = 6 / .pt,
               color = "#a3a29d", fontface = "italic") +
      xlim(0, 1) + ylim(0, 1) +
      theme_void() +
      theme(plot.margin = margin(1, 1, 1, 1))
  }

  p_5a <- placeholder_panel("structure overview\n(Illustrator)")
  p_5b <- placeholder_panel("SERK-binding interface\n(Illustrator)")
  p_5c <- placeholder_panel("(TBD)")

  # --- Panel D: BLAST heatmap, no-close-paralog genes vs. top SERK/CIK hits ---
  heat <- read_csv(fig6_heat_path, show_col_types = FALSE) |>
    mutate(col_gene = factor(col_gene, levels = unique(col_gene)),
           row_gene = factor(row_gene, levels = rev(unique(row_gene))))

  col_label_colors <- heat |>
    distinct(col_gene, group) |>
    arrange(col_gene) |>
    mutate(color = if_else(group == "SERK", SERK_COLOR, CIK_COLOR)) |>
    pull(color)

  p_5d <- ggplot(heat, aes(col_gene, row_gene, fill = bitscore)) +
    geom_tile() +
    scale_fill_viridis_c(name = "BLAST\nbitscore",
                         guide = guide_colorbar(barwidth = unit(6, "pt"), barheight = unit(28, "pt"))) +
    base_theme() +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6, color = col_label_colors),
          axis.text.y = element_text(size = 6, face = "bold", color = SERK_COLOR),
          legend.title = element_text(size = 6), legend.text = element_text(size = 6),
          # extra bottom margin left blank for a hand-added SERK-binding /
          # non-SERK-binding (or more specific) label under the x-axis
          plot.margin = margin(2, 8, 26, 2))

  # --- Panel E: cross-species densities, cosine 1-NN margin then CAV score ---
  # Both facets are 1-D densities of a signed score, so they read in parallel.
  # cosine 1-NN margin = cosine to the nearest SERK-partner training gene
  # minus cosine to the nearest CIK-partner training gene. Stated as a
  # similarity margin rather than a distance margin so that, like the CAV
  # score, positive = more SERK-like (a distance margin would invert the sign).
  # Both use the same in-context pooled embeddings. Regenerate the CSVs with
  # interface_search/fastas/pairs/lrr_serkcik/export_fig5E_data.py.
  # The scatter decomposition of this margin (which shows that all 36 genes sit
  # at cosine > 0.985 to both classes) is kept in fig_supp_noparalog_loo.
  dens <- read_csv(fig6_dens_path, show_col_types = FALSE) |>
    mutate(label = factor(label, levels = c("SERK-partner", "CIK-partner")))
  auroc <- read_csv(fig6_auroc_path, show_col_types = FALSE)
  auroc_cos <- auroc$auroc[auroc$metric == "cosine_margin"]
  auroc_cav <- auroc$auroc[auroc$metric == "cav_score"]

  pal <- c("SERK-partner" = SERK_COLOR, "CIK-partner" = CIK_COLOR)

  dens_long <- dens |>
    pivot_longer(cols = c(cosine_margin, cav_score), names_to = "metric", values_to = "value") |>
    mutate(metric = if_else(metric == "cosine_margin", "cosine 1-NN margin", "CAV score"),
           metric = factor(metric, levels = c("cosine 1-NN margin", "CAV score")))

  auroc_df <- tibble(
    metric = factor(c("cosine 1-NN margin", "CAV score"), levels = c("cosine 1-NN margin", "CAV score")),
    auroc_text = paste0("AUROC=", sprintf("%.3f", c(auroc_cos, auroc_cav)))
  )

  # geom_density(trim = TRUE): stop each group's curve at its own data range
  # instead of extrapolating a long near-zero tail out to the expanded limits.
  p_5e <- ggplot(dens_long, aes(value, color = label, fill = label)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "#b9b8b3", linetype = "22") +
    geom_density(alpha = 0.18, linewidth = 0.6, trim = TRUE) +
    geom_rug(linewidth = 0.4, alpha = 0.9, length = unit(0.03, "npc")) +
    scale_color_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.08)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    facet_wrap(~metric, scales = "free", strip.position = "bottom") +
    geom_text(data = auroc_df, aes(x = Inf, y = Inf, label = auroc_text),
              inherit.aes = FALSE, hjust = 1.1, vjust = 1.8, size = 6 / .pt, fontface = "bold") +
    labs(x = NULL, y = NULL) +
    base_theme() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 6),
          strip.background = element_blank(), strip.placement = "outside",
          strip.text = element_text(size = 7),
          panel.spacing = unit(14, "pt"),
          legend.position = "none",
          plot.margin = margin(10, 14, 10, 14))


  legend_5 <- get_legend(
    p_5e + theme(legend.position = "bottom", legend.key.size = unit(7, "pt"),
                 legend.text = element_text(size = 7))
  )

  row1 <- plot_grid(p_5a, p_5b, nrow = 1, labels = c("A", "B"), label_size = 8,
                    label_fontfamily = FIG_FONT)
  row2 <- plot_grid(p_5c, p_5d, nrow = 1, rel_widths = c(0.6, 0.4),
                     labels = c("C", "D"), label_size = 8, label_fontfamily = FIG_FONT)
  row3_plots <- plot_grid(p_5e, labels = "E", label_size = 8, label_fontfamily = FIG_FONT)
  row3 <- plot_grid(row3_plots, legend_5, ncol = 1, rel_heights = c(1, 0.1))

  fig5 <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(1, 1.15, 1.15))

  # 5.5in ~ 140mm, a standard "1.5-column" journal figure width; height
  # rescaled to preserve the original 5:7 aspect ratio.
  ggsave(file.path(OUT, "cav_fig5.pdf"), fig5, width = 5.5, height = 7.7, bg = "white", device = cairo_pdf)
  ggsave(file.path(OUT, "cav_fig5.png"), fig5, width = 5.5, height = 7.7, bg = "white", dpi = 300)
  message("Saved cav_fig5.pdf / cav_fig5.png")

  # -------------------------------------------------------------------------
  # Supplemental: within-Arabidopsis leave-one-out, cosine 1-NN vs. matched CAV
  # -------------------------------------------------------------------------
  # A DIFFERENT experiment from Figure 5E (which is Arabidopsis-train ->
  # non-Arabidopsis-test). Here both scores are leave-one-out over the same 35
  # Arabidopsis genes and the same in-context embeddings, so the two AUROCs are
  # directly comparable: cosine 0.876 vs. CAV 0.850. Panel D's four
  # no-close-paralog genes are labelled; all four sit on the y = x boundary,
  # i.e. their nearest-neighbour class calls are ties.
  # Regenerate with export_fig5E_noparalog.py.
  np_path <- file.path(DATA, "fig6_noparalog_scatter.csv")

  if (file.exists(np_path)) {
    npd <- read_csv(np_path, show_col_types = FALSE) |>
      mutate(label = factor(label, levels = c("SERK-partner", "CIK-partner")))
    np_lim <- range(c(npd$cos_nearest_serk, npd$cos_nearest_cik))
    np_lim <- np_lim + c(-1, 1) * 0.04 * diff(np_lim)

    auroc_np_cos <- 0.876
    auroc_np_cav <- 0.850

    p_np_cos <- ggplot(npd, aes(cos_nearest_cik, cos_nearest_serk, color = label)) +
      geom_abline(slope = 1, intercept = 0, linewidth = 0.35,
                  color = "#8a8985", linetype = "22") +
      geom_point(data = ~ dplyr::filter(.x, !panel_d), size = 1.1, alpha = 0.45) +
      geom_point(data = ~ dplyr::filter(.x, panel_d), size = 2.1, alpha = 1) +
      ggrepel::geom_text_repel(data = ~ dplyr::filter(.x, panel_d),
                               aes(label = gene), size = 7 / .pt,
                               segment.linewidth = 0.25, min.segment.length = 0,
                               box.padding = 0.4, show.legend = FALSE) +
      scale_color_manual(values = pal, name = NULL) +
      coord_fixed(xlim = np_lim, ylim = np_lim) +
      annotate("text", x = Inf, y = -Inf,
               label = paste0("cosine 1-NN  AUROC=", sprintf("%.3f", auroc_np_cos)),
               hjust = 1.05, vjust = -1.0, size = 7 / .pt, fontface = "bold") +
      labs(x = "cosine to nearest CIK-partner", y = "cosine to nearest SERK-partner") +
      base_theme() +
      theme(legend.position = "none", plot.margin = margin(8, 8, 4, 8))

    # Full range, no clipping: LOO CAV scores span -26.9 to +31.2 (RLP23 and
    # CLV2, the two kinase-domain-lacking receptor-like proteins, are the
    # documented LOO outliers). Clipping was tried and rejected -- it hid 22 of
    # 35 points. The wide spread and class overlap are the honest picture at
    # AUROC 0.850.
    p_np_cav <- ggplot(npd, aes(cav_loo, label, color = label)) +
      geom_vline(xintercept = 0, linewidth = 0.35, color = "#8a8985", linetype = "22") +
      geom_point(size = 1.9, alpha = 0.85,
                 position = position_jitter(height = 0.16, width = 0, seed = 1)) +
      scale_color_manual(values = pal, name = NULL) +
      scale_x_continuous(expand = expansion(mult = 0.06)) +
      annotate("text", x = Inf, y = Inf,
               label = paste0("CAV (refit LOO)  AUROC=", sprintf("%.3f", auroc_np_cav)),
               hjust = 1.05, vjust = 1.8, size = 7 / .pt, fontface = "bold") +
      labs(x = "CAV score (leave-one-out)", y = NULL) +
      base_theme() +
      theme(legend.position = "none", plot.margin = margin(8, 8, 4, 8))

    legend_np <- get_legend(
      p_np_cos + theme(legend.position = "bottom", legend.key.size = unit(8, "pt"),
                       legend.text = element_text(size = 8))
    )

    fig_np <- plot_grid(
      plot_grid(p_np_cos, p_np_cav, nrow = 1, rel_widths = c(1, 1),
                labels = c("A", "B"), label_size = 9, label_fontfamily = FIG_FONT),
      legend_np, ncol = 1, rel_heights = c(1, 0.09))

    ggsave(file.path(OUT, "fig_supp_noparalog_loo.pdf"), fig_np,
           width = 7.2, height = 3.6, bg = "white", device = cairo_pdf)
    ggsave(file.path(OUT, "fig_supp_noparalog_loo.png"), fig_np,
           width = 7.2, height = 3.6, bg = "white", dpi = 300)
    message("Saved fig_supp_noparalog_loo.pdf / .png")
  }

} else {
  message("Skipping Figure 5: run the export step in ",
          "interface_search/fastas/pairs/lrr_serkcik/ to populate ",
          "figure_data/fig6_heatmap.csv, fig6_densities.csv, fig6_auroc.csv")
}

message("\nAll figures written to ", OUT, "/")
