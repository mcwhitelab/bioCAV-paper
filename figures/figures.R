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
})

RASTER_DPI <- 150

DATA <- "figure_data"
OUT  <- "figures"
dir.create(OUT, showWarnings = FALSE)

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
base_theme <- function(...) theme_cowplot(font_size = 8, ...)

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
# Panel A: intro panel -- CAV score (z-scored against each term's own
# negative/background distribution) for held-out positive vs. negative
# examples, pooled across ontologies. This is the "does the method work at
# all" panel: negatives sit centered near z=0 (they define the reference
# frame), positives shift well above it. Source: figure_data/
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
    labs(x = "CAV score (z, vs. own-term negative background)", y = "Density") +
    theme(
      legend.position = "top",
      legend.key.size = unit(7, "pt"),
      legend.text     = element_text(size = 6),
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
        legend.key.size   = unit(7, "pt"),
        legend.text       = element_text(size = 6),
        legend.margin     = margin(t = 0, b = 0),
        strip.background  = element_blank(),
        strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0),
        panel.spacing.y   = unit(3, "pt"),
        plot.margin       = margin(t = 2, r = 10, b = 2, l = 2)
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
      guides(fill = guide_legend(ncol = 5, byrow = TRUE, keywidth = unit(5, "pt"))) +
      theme(
        legend.position   = "top",
        legend.key.size   = unit(5, "pt"),
        legend.text       = element_text(size = 5),
        legend.spacing.x  = unit(1.5, "pt"),
        legend.margin     = margin(t = 0, b = 0),
        strip.background  = element_blank(),
        strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0, size = 7),
        panel.spacing.y   = unit(4, "pt"),
        plot.margin       = margin(t = 2, r = 10, b = 2, l = 2)
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
      legend.text     = element_text(size = 6),
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
               size = 2.3, hjust = 1.05, vjust = 1.3) +
      scale_fill_manual(values = c("Annotated match" = GE_MATCH_COLOR), name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      base_theme() +
      labs(x = paste0("Rank of annotated EC partner\n(out of ", n_ec_cavs, " trained EC CAVs)"),
           y = "Number of GO terms") +
      GE_LEGEND_THEME +
      theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5))

    # align = "h", axis = "tb": match the top/bottom plot-panel edges of E
    # and F so their x-axes sit at the same height despite F's rotated tick
    # labels taking more vertical space than E's horizontal ones.
    p_2e <- plot_grid(p_2e_solo, p_2f_solo, nrow = 1, align = "h", axis = "tb",
                      labels = c("E", "F"), label_size = 8,
                      rel_widths = c(1, 1))
  } else {
    message("Skipping GO-EC cosine validation: figure_data/go_ec_cosine_*.csv not found")
  }

  # ---------------------------------------------------------------------------
  # Figure 3, panel B: CAV projection score along the EC hierarchy (level
  # 1->4) for a handful of validation proteins with a full lineage. Not a
  # cosine-similarity metric -- shows the score sharpening as the CAV gets
  # more specific down the hierarchy. Source: hierarchy_decay/results/
  # ec_hierarchy_decay.csv (copied to figure_data/).
  # ---------------------------------------------------------------------------
  ec_decay_path <- file.path(DATA, "ec_hierarchy_decay.csv")

  p_3b <- NULL
  if (file.exists(ec_decay_path)) {
    ec_decay <- read_csv(ec_decay_path, show_col_types = FALSE)

    p_3b <- ec_decay |>
      mutate(level = factor(level, levels = 1:4,
                            labels = c("Level 1", "Level 2", "Level 3", "Level 4\n(fully specific)"))) |>
      ggplot(aes(x = level, y = score, group = protein_id, color = protein_id)) +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      geom_point(size = 1.2) +
      scale_color_manual(values = scales::hue_pal()(n_distinct(ec_decay$protein_id)), name = "Protein") +
      base_theme() +
      labs(x = "EC hierarchy depth", y = "CAV projection score") +
      theme(legend.position = "right", legend.key.size = unit(7, "pt"),
            legend.text = element_text(size = 6), legend.title = element_text(size = 7))
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
      guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
      base_theme() +
      labs(x = "UMAP 1", y = "UMAP 2",
           title = paste0("Combined CAV space (n=", nrow(combined_umap), ")")) +
      theme(plot.title = element_text(size = 8, face = "bold"),
            legend.position = "right", legend.key.size = unit(8, "pt"),
            legend.text = element_text(size = 6), axis.text = element_blank(),
            axis.ticks = element_blank())
  } else {
    message("Skipping Figure 3 combined-embedding panel: figure_data/go_ec_combined_umap.csv not found ",
            "(run motif_clustering/combined_cav_umap_3d.py with --coords-csv-out first)")
  }

  # ---------------------------------------------------------------------------
  # Assemble Figure 2 — GO evaluation. Stacked rows: A & B (intro); C
  # (rank-tier composition); D & E (AUC/AUPR violins).
  # ---------------------------------------------------------------------------
  intro_row <- plot_grid(
    if (!is.null(p_2a_zscore)) p_2a_zscore else ggplot() + theme_void(),
    NULL,
    nrow       = 1,
    labels     = c("A", "B"),
    label_size = 8
  )

  legend_method <- get_legend(
    p_2b + theme(legend.position = "bottom", legend.key.size = unit(9, "pt"))
  )

  de_row <- plot_grid(
    p_2b, p_2c,
    nrow       = 1,
    labels     = c("D", "E"),
    label_size = 8
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
    rel_heights <- c(rel_heights, 1.05)
  }
  rows        <- c(rows, list(violin_block))
  row_labels  <- c(row_labels, "")
  rel_heights <- c(rel_heights, 1)

  fig2 <- plot_grid(
    plotlist    = rows,
    ncol        = 1,
    labels      = row_labels,
    label_size  = 8,
    rel_heights = rel_heights
  )

  ggsave(file.path(OUT, "fig2.pdf"), fig2, width = 4.2, height = 7.2)
  message("Saved fig2.pdf")

  # ---------------------------------------------------------------------------
  # Assemble Figure 3 — EC evaluation & cross-ontology validation. Stacked
  # rows: A (EC rank composition); B (EC hierarchy-depth score); C & D
  # (GO-EC cosine validation); E (combined 2D CAV overview).
  # ---------------------------------------------------------------------------
  fig3_rows        <- list()
  fig3_row_labels  <- character(0)
  fig3_rel_heights <- numeric(0)

  if (!is.null(p_2d)) {
    fig3_rows        <- c(fig3_rows, list(p_2d))
    fig3_row_labels  <- c(fig3_row_labels, "A")
    fig3_rel_heights <- c(fig3_rel_heights, 0.85)
  }
  if (!is.null(p_3b)) {
    fig3_rows        <- c(fig3_rows, list(p_3b))
    fig3_row_labels  <- c(fig3_row_labels, "B")
    fig3_rel_heights <- c(fig3_rel_heights, 0.85)
  }
  if (!is.null(p_2e)) {
    p_2e_relabeled <- plot_grid(p_2e_solo, p_2f_solo, nrow = 1, align = "h", axis = "tb",
                                labels = c("C", "D"), label_size = 8,
                                rel_widths = c(1, 1))
    fig3_rows        <- c(fig3_rows, list(p_2e_relabeled))
    fig3_row_labels  <- c(fig3_row_labels, "")
    fig3_rel_heights <- c(fig3_rel_heights, 0.85)
  }
  if (!is.null(p_3e)) {
    fig3_rows        <- c(fig3_rows, list(p_3e))
    fig3_row_labels  <- c(fig3_row_labels, "E")
    fig3_rel_heights <- c(fig3_rel_heights, 1.1)
  }

  if (length(fig3_rows) > 0) {
    fig3 <- plot_grid(
      plotlist    = fig3_rows,
      ncol        = 1,
      labels      = fig3_row_labels,
      label_size  = 8,
      rel_heights = fig3_rel_heights
    )

    ggsave(file.path(OUT, "fig3.pdf"), fig3, width = 4.6,
           height = sum(fig3_rel_heights) * 2.1, bg = "white")
    message("Saved fig3.pdf")
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
    geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 20) +
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
      align      = "hv",
      axis       = "tblr"
    )
    ggsave(file.path(OUT, "fig_ec_eval.pdf"), fig_ec,
           width = 3.0 * length(ec_summary_panels), height = 3.0)
    message("Saved fig_ec_eval.pdf")
  }
}

# ===========================================================================
# Figure 5 — single-cell case studies (DE-vs-CAV scatter + transcriptional
# continuum, one column per case). Data comes from
# single_cell/scripts/export_fig5_data.py (run via
# single_cell/08_export_fig5_data.sh), which reuses cav_continuum_viz.py's
# cell/gene selection so this matches the exploratory PNGs it replaces
# (single_cell/07_paper_cases.sh) exactly, just rendered as native vector
# ggplot panels instead of embedded matplotlib rasters.
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
    "epithelial_cell__lung10x__normal_vs_lung_cancer"           = "Epithelial cell — lung (10x-matched)",
    "fibroblast__colorectum__normal_vs_colorectal_cancer"       = "Fibroblast — colorectum"
  )

  SC_NORMAL_COLOR     <- "#3a6fad"
  SC_CANCER_COLOR     <- "#c0392b"
  SC_BACKGROUND_COLOR <- "#2277bb"
  SC_CONTINUUM_FILL   <- "#f0a500"
  SC_CONTINUUM_EDGE   <- "#7a5300"

  # ---------------------------------------------------------------------
  # Overview row: t-SNE of the CAV hierarchy's three levels -- L0 (cell_type
  # axes), L1 (tissue axes), L2 (condition-residual axes), per
  # cav_hierarchy.py's default group-context-condition level order. Each
  # panel is colored by whatever that level is actually built to separate.
  # Data comes from single_cell/scripts/export_fig5_umap_data.py (run via
  # single_cell/09_export_fig5_umap_data.sh).
  # ---------------------------------------------------------------------
  sc_umap_path <- file.path(DATA, "sc_umap_overview.csv")

  sc_overview_row <- NULL
  if (file.exists(sc_umap_path)) {
    sc_umap <- read_csv(sc_umap_path, show_col_types = FALSE)

    sc_umap_theme <- function() {
      base_theme() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              plot.title = element_text(size = 8, face = "bold"),
              plot.subtitle = element_text(size = 7))
    }

    sc_umap <- sc_umap %>% mutate(condition = if_else(is_baseline, "normal", "cancer"))

    sc_condition_scale <- scale_color_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR),
                                             name = NULL)
    sc_condition_guide <- guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

    sc_raw_df <- sc_umap %>% filter(level == "raw")
    p_umap_raw <- ggplot(sc_raw_df, aes(x, y, color = condition)) +
      geom_point_rast(size = 0.15, alpha = 0.45, raster.dpi = RASTER_DPI) +
      sc_condition_scale + sc_condition_guide +
      sc_umap_theme() +
      theme(legend.text = element_text(size = 6), legend.key.size = unit(7, "pt")) +
      labs(x = "t-SNE 1", y = "t-SNE 2", title = "Raw — individual baseline CAVs",
           subtitle = "colored by condition")

    sc_l1_df <- sc_umap %>% filter(level == "L1")
    p_umap_l1 <- ggplot(sc_l1_df, aes(x, y, color = condition)) +
      geom_point_rast(size = 0.15, alpha = 0.45, raster.dpi = RASTER_DPI) +
      sc_condition_scale + sc_condition_guide +
      sc_umap_theme() +
      theme(legend.text = element_text(size = 6), legend.key.size = unit(7, "pt")) +
      labs(x = "t-SNE 1", y = "t-SNE 2", title = "L1 — cell type axis (tissue projected away)",
           subtitle = "colored by condition")

    sc_l2_df <- sc_umap %>% filter(level == "L2")
    p_umap_l2 <- ggplot(sc_l2_df, aes(x, y, color = condition)) +
      geom_point_rast(size = 0.15, alpha = 0.45, raster.dpi = RASTER_DPI) +
      sc_condition_scale + sc_condition_guide +
      sc_umap_theme() +
      theme(legend.text = element_text(size = 6), legend.key.size = unit(7, "pt")) +
      labs(x = "t-SNE 1", y = "t-SNE 2",
           title = "L2 — condition axis (tissue + cell type away)",
           subtitle = "colored by condition")

    sc_overview_row <- plot_grid(NULL, p_umap_raw, p_umap_l1, p_umap_l2,
                                 nrow = 1, rel_widths = c(1, 1, 1, 1),
                                 labels = c("A", "B", "C", "D"), label_size = 9)
  } else {
    message("Figure 5 overview row skipped: figure_data/sc_umap_overview.csv not found ",
            "(run single_cell/09_export_fig5_umap_data.sh first)")
  }

  # ---------------------------------------------------------------------
  # DE-vs-CAV scatter, one per case. Every gene is plotted -- genes below
  # the |log2FC| threshold are shown as a faint background cloud for
  # context, but the Spearman r annotation (and the auto-labelled genes)
  # is still computed on the same "meaningful" subset as
  # 07_paper_cases.sh's de_vs_cav_scatter_paper_cases.png (passes the DE
  # magnitude threshold, OR is one of that case's top continuum genes) so
  # the stat doesn't get diluted by the thousands of near-zero-log2FC genes.
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

    df_stat <- df %>% filter(category != "below_threshold")
    rho <- suppressWarnings(cor(df_stat$log2fc, df_stat$cav_r, method = "spearman"))

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
           subtitle = sprintf("Spearman r = %.2f", rho)) +
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
  sc_ramp_pos <- colorRampPalette(brewer.pal(9, "YlOrRd"))(101)
  sc_ramp_neg <- colorRampPalette(brewer.pal(9, "YlGnBu"))(101)

  sc_expr_to_hex <- function(expr_scaled, is_cancer) {
    idx <- pmin(pmax(round((0.2 + 0.75 * expr_scaled) * 100) + 1, 1), 101)
    ifelse(is_cancer, sc_ramp_pos[idx], sc_ramp_neg[idx])
  }

  sc_build_continuum <- function(pair_id, seed = 1) {
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
    xlim <- c(x_lo - pad, x_hi + pad)

    p_strip <- ggplot(cells, aes(l2_score, y, color = is_baseline)) +
      geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
      geom_point_rast(size = 0.4, alpha = 0.6, raster.dpi = RASTER_DPI) +
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

    genes$y_off <- if_else(genes$is_baseline, 0.22, -0.22)
    genes$hex <- sc_expr_to_hex(genes$expr_scaled, !genes$is_baseline)

    gene_labels <- genes %>% distinct(y_ctr, gene_name, r) %>%
      mutate(lab = sprintf("%s (r=%.2f)", gene_name, r)) %>% arrange(y_ctr)

    p_genes <- ggplot(genes, aes(l2_score, y_ctr + y_off + jit_cell)) +
      geom_hline(yintercept = seq(0.5, n_genes + 0.5, 1), color = "grey90", linewidth = 0.3) +
      geom_hline(yintercept = seq(1, n_genes, 1), color = "grey95", linewidth = 0.2) +
      {if (x_lo < 0 && x_hi > 0) geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4, linetype = "dashed")} +
      geom_point_rast(aes(color = hex), size = 0.35, alpha = 0.6, raster.dpi = RASTER_DPI) +
      scale_color_identity() +
      scale_y_continuous(breaks = gene_labels$y_ctr, labels = gene_labels$lab,
                          limits = c(0.5, n_genes + 0.5), expand = c(0, 0)) +
      coord_cartesian(xlim = xlim) +
      base_theme() +
      labs(x = "L2 score  (← normal    cancer →)", y = NULL) +
      theme(axis.text.y  = element_text(size = 6),
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank(),
            plot.margin  = margin(t = 0, b = 2, l = 2, r = 2))

    plot_grid(p_strip, p_genes, ncol = 1, align = "v", axis = "lr",
              rel_heights = c(0.9, n_genes))
  }

  sc_columns <- map(SC_PAIR_ORDER, function(pid) {
    plot_grid(sc_build_scatter(pid), sc_build_continuum(pid),
              ncol = 1, align = "v", axis = "lr",
              rel_heights = c(0.62, 1))
  })

  case_labels <- if (is.null(sc_overview_row)) c("A", "B", "C") else c("E", "F", "G")
  fig5_body <- plot_grid(plotlist = sc_columns, nrow = 1,
                         labels = case_labels, label_size = 9)

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
  sc_legend_scatter <- get_legend(
    ggplot(sc_legend_scatter_df, aes(x, y, color = cat)) +
      geom_point(size = 2.5) +
      scale_color_manual(values = setNames(c(SC_BACKGROUND_COLOR, SC_CONTINUUM_EDGE),
                                           levels(sc_legend_scatter_df$cat)),
                         name = NULL) +
      base_theme() + theme(legend.position = "right", legend.key.size = unit(9, "pt"),
                            legend.text = element_text(size = 7))
  )

  sc_legend_row <- plot_grid(NULL, sc_legend_scatter, sc_legend_strip, NULL,
                             nrow = 1, rel_widths = c(0.3, 1.4, 1, 0.3))

  fig5_cases <- plot_grid(fig5_body, sc_legend_row, ncol = 1, rel_heights = c(1, 0.09))

  if (!is.null(sc_overview_row)) {
    fig5 <- plot_grid(sc_overview_row, fig5_cases, ncol = 1, rel_heights = c(0.6, 1))
    fig5_height <- 10.5
  } else {
    fig5 <- fig5_cases
    fig5_height <- 7
  }

  fig5_width <- 12
  ggsave(file.path(OUT, "fig5.pdf"), fig5, width = fig5_width, height = fig5_height,
         bg = "white", device = cairo_pdf)
  ggsave(file.path(OUT, "fig5.png"), fig5, width = fig5_width, height = fig5_height,
         bg = "white", dpi = 300)
  message("Saved fig5.pdf / fig5.png")

  # -------------------------------------------------------------------
  # Candidate supplemental — does orthogonalizing the condition CAV
  # against cell-type/tissue baseline structure change the DE-vs-CAV
  # relationship, compared to using the raw (unprocessed) condition CAV?
  # Three rows = three stages of the SAME subtraction chain that produces
  # L2: raw -> minus cell-type baseline (L0) -> minus cell-type + tissue
  # (L0+L1 = L2, identical to the L2 axis used everywhere else in Figure
  # 5). There is no fourth stage to subtract in this two-level hierarchy.
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
        mutate(above_threshold = abs(log2fc) >= 1.5)

      df_stat <- df %>% filter(above_threshold)
      # fibroblast/colorectum has no genes clearing |log2FC| >= 1.5 at all
      # (max |log2FC| ~1.27) -- fall back to all genes rather than an
      # undefined stat for that pair.
      if (nrow(df_stat) < 5) df_stat <- df
      rho <- suppressWarnings(cor(df_stat$log2fc, df_stat$cav_r, method = "spearman"))

      nx <- df_stat$log2fc / (max(abs(df_stat$log2fc)) + 1e-9)
      ny <- df_stat$cav_r   / (max(abs(df_stat$cav_r))   + 1e-9)
      df_stat$dist <- sqrt(nx^2 + ny^2)
      label_ids <- df_stat %>% slice_max(dist, n = 6) %>% pull(gene)
      df$show_label <- df$gene %in% label_ids

      ggplot(df, aes(log2fc, cav_r)) +
        geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
        geom_point_rast(data = filter(df, !above_threshold),
                   color = "grey70", size = 0.4, alpha = 0.15, raster.dpi = RASTER_DPI) +
        geom_point_rast(data = filter(df, above_threshold),
                   color = SC_BACKGROUND_COLOR, size = 0.7, alpha = 0.4, raster.dpi = RASTER_DPI) +
        geom_text_repel(
          data = filter(df, show_label),
          aes(label = gene_name), color = "grey20",
          size = 1.9, max.overlaps = 30, segment.size = 0.25,
          segment.color = "grey60",
          box.padding = 0.3, point.padding = 0.15, force = 2, force_pull = 0.5,
          min.segment.length = 0.1, seed = 42
        ) +
        base_theme() +
        labs(x = "log2FC (DE)", y = "Pearson r (CAV)",
             subtitle = sprintf("Spearman r = %.2f", rho)) +
        theme(plot.subtitle = element_text(size = 7))
    }

    proj_col_headers <- plot_grid(plotlist = map(SC_PAIR_ORDER, function(pid) {
      ggdraw() + draw_label(SC_PAIR_TITLES[pid], fontface = "bold", size = 8)
    }), nrow = 1)

    proj_row_label <- function(lvl) {
      ggdraw() + draw_label(PROJ_LEVEL_TITLES[lvl], fontface = "bold", size = 7.5, angle = 90)
    }

    proj_body_rows <- map(seq_along(PROJ_LEVEL_ORDER), function(i) {
      lvl <- PROJ_LEVEL_ORDER[i]
      row_plots <- map(SC_PAIR_ORDER, function(pid) proj_build_scatter(pid, lvl))
      plot_grid(proj_row_label(lvl), plot_grid(plotlist = row_plots, nrow = 1),
                nrow = 1, rel_widths = c(0.07, 1))
    })

    proj_body <- plot_grid(plotlist = proj_body_rows, ncol = 1,
                           labels = LETTERS[1:3], label_size = 9)

    fig_supp_projection <- plot_grid(
      plot_grid(NULL, proj_col_headers, nrow = 1, rel_widths = c(0.07, 1)),
      proj_body, ncol = 1, rel_heights = c(0.05, 1)
    )

    ggsave(file.path(OUT, "fig_supp_projection_comparison.pdf"), fig_supp_projection,
           width = 11, height = 9.5, bg = "white", device = cairo_pdf)
    ggsave(file.path(OUT, "fig_supp_projection_comparison.png"), fig_supp_projection,
           width = 11, height = 9.5, bg = "white", dpi = 300)
    message("Saved fig_supp_projection_comparison.pdf / .png")

  } else {
    message("Skipping projection-comparison supplemental: figure_data/proj_de_vs_cav.csv not found ",
            "(run single_cell/10_export_projection_comparison_supp_data.sh first)")
  }

} else {
  message("Skipping Figure 5: figure_data/sc_pair_meta.csv not found ",
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

  MODULE_COLORS <- c(quiescent = "#3a6fad", activated = "#c0392b", independent = "#7a5300")
  MODULE_LABELS <- c(quiescent = "quiescent-fibroblast module (MGP/DCN/OGN/C3/CCDC80)",
                     activated = "CXCR4 (anti-correlated with the module)",
                     independent = "ADAMDEC1 (independent of the module)")

  # --- Panel A: DE log2FC vs CAV r, background cloud + highlighted genes ---
  p_subpop_a <- ggplot(subpop_bg, aes(log2fc, cav_r)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_point_rast(color = "grey70", size = 0.5, alpha = 0.25, raster.dpi = RASTER_DPI) +
    geom_point(data = subpop_genes, aes(color = module), size = 2.2) +
    ggrepel::geom_text_repel(data = subpop_genes, aes(label = gene_name, color = module),
                             size = 2.3, fontface = "bold", show.legend = FALSE,
                             box.padding = 0.4, seed = 42) +
    scale_color_manual(values = MODULE_COLORS, labels = MODULE_LABELS, name = NULL) +
    base_theme() +
    labs(x = "log2FC (DE, mixedlm)", y = "Pearson r (CAV)",
         title = "Fibroblast — colorectum: DE-null, CAV-strong genes") +
    theme(legend.position = "bottom", legend.text = element_text(size = 6),
          legend.key.size = unit(7, "pt"),
          plot.title = element_text(size = 8, face = "bold")) +
    guides(color = guide_legend(nrow = 3))

  # --- Panel B: co-expression heatmap among the 7 genes ---
  gene_order <- c("MGP", "DCN", "OGN", "C3", "CCDC80", "CXCR4", "ADAMDEC1", "CXCL14")
  p_subpop_b <- subpop_corr |>
    mutate(gene1 = factor(gene1, levels = gene_order),
           gene2 = factor(gene2, levels = rev(gene_order))) |>
    ggplot(aes(gene1, gene2, fill = r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.2) +
    scale_fill_gradient2(low = unname(oi["vermillion"]), mid = "white", high = unname(oi["blue"]),
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    base_theme() +
    labs(x = NULL, y = NULL, title = "Pairwise co-expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size = 8, face = "bold"),
          legend.key.size = unit(9, "pt"))

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
      `Module-negative\n(all 5 markers = 0)` = mean(module_negative),
      `High module score\n(top quartile)`     = mean(high_quiescent),
      .groups = "drop"
    ) |>
    pivot_longer(-disease, names_to = "metric", values_to = "frac")

  p_subpop_c <- composition_df |>
    ggplot(aes(metric, frac, fill = disease)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = scales::percent(frac, accuracy = 1)),
              position = position_dodge(width = 0.7), vjust = -0.4, size = 2.6) +
    scale_fill_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR), name = NULL) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
    base_theme() +
    labs(x = NULL, y = "Share of fibroblasts",
         title = "Compositional shift: real, but modest") +
    theme(legend.position = "top", plot.title = element_text(size = 8, face = "bold"))

  # --- Panel D: quiescent-module score vs. the continuous CAV L2 axis ---
  p_subpop_d <- subpop_cells |>
    mutate(disease = if_else(str_detect(tolower(disease), "normal"), "normal", "cancer"),
           disease = factor(disease, levels = c("normal", "cancer"))) |>
    ggplot(aes(l2_score, quiescent_score, color = disease)) +
    geom_point_rast(size = 0.6, alpha = 0.4, raster.dpi = RASTER_DPI) +
    geom_smooth(aes(group = 1), method = "loess", color = "black", linewidth = 0.6, se = TRUE) +
    scale_color_manual(values = c(normal = SC_NORMAL_COLOR, cancer = SC_CANCER_COLOR), name = NULL) +
    base_theme() +
    labs(x = "L2 score  (← normal    cancer →)", y = "Quiescent-module score",
         title = "Module score tracks the continuous CAV axis") +
    theme(legend.position = "top", plot.title = element_text(size = 8, face = "bold"))

  fig_subpop <- plot_grid(p_subpop_a, p_subpop_b, p_subpop_c, p_subpop_d,
                          nrow = 2, labels = "AUTO", label_size = 9,
                          align = "hv", axis = "tblr")

  ggsave(file.path(OUT, "fig_supp_fibroblast_subpop.pdf"), fig_subpop,
         width = 8.5, height = 8, bg = "white", device = cairo_pdf)
  message("Saved fig_supp_fibroblast_subpop.pdf")

} else {
  message("Skipping fibroblast-subpopulation supplemental: figure_data/subpop_gene_scatter.csv not found ",
          "(run single_cell/scripts/export_subpop_supp_data.py first)")
}

message("\nAll figures written to ", OUT, "/")
