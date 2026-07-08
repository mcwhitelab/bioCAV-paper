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
})

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
# Figure 2 — GO and EC evaluation
# ===========================================================================

# Load and combine all three ontologies (used by panels A, B, C)
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

  # ---------------------------------------------------------------------------
  # Panel A: GO term specificity rank histograms, faceted by ontology
  # For each (protein, GO term) validation pair, both methods rank the true
  # GO term among all trained CAVs.  x = rank / n_go_terms (0 = top, 1 = last).
  # ggridges stat = "binline" — same histogram-ridge style as the per-prediction
  # split plots in compare_tool_temporal.py.
  # ---------------------------------------------------------------------------
  ranks_all <- bind_rows(
    load_ont("go_specificity_ranks", "mf"),
    load_ont("go_specificity_ranks", "bp"),
    load_ont("go_specificity_ranks", "cc")
  )

  if (nrow(ranks_all) > 0) {
    bins_n     <- 40
    bar_w      <- 1 / bins_n   # 0.025
    EXCL_X     <- 1.10         # center of the excluded bar
    EXCL_SEP   <- 1.025        # dotted separator line position

    # Denominator = all pairs per ontology (matches Python's n_total)
    n_total_df <- ranks_all |> count(ontology, name = "n_total")

    # Rank histogram bins — CAV uses llr > 0; DeepGoSE uses tool_predicted == TRUE
    cav_bins <- ranks_all |>
      filter(llr > 0) |>
      mutate(bin = floor(cav_rank / n_go_terms * bins_n) / bins_n,
             method = "CAV") |>
      count(method, ontology, bin)

    tool_bins <- ranks_all |>
      filter(tool_predicted == TRUE) |>
      mutate(bin = floor(tool_rank / n_go_terms * bins_n) / bins_n,
             method = "DeepGoSE") |>
      count(method, ontology, bin)

    # Excluded bars (lighter alpha, placed after gap)
    cav_excl <- ranks_all |>
      filter(llr <= 0) |>
      count(ontology) |>
      mutate(method = "CAV",      bin = EXCL_X - 0.5 * bar_w)

    tool_excl <- ranks_all |>
      filter(!tool_predicted) |>
      count(ontology) |>
      mutate(method = "DeepGoSE", bin = EXCL_X - 0.5 * bar_w)

    rank_binned <- bind_rows(cav_bins, tool_bins, cav_excl, tool_excl) |>
      left_join(n_total_df, by = "ontology") |>
      mutate(
        prop      = n / n_total,
        bar_alpha = 0.75,
        method    = factor(method,   levels = c("CAV", "DeepGoSE")),
        ontology  = factor(ontology, levels = c("MF", "BP", "CC"))
      )

    p_2a <- rank_binned |>
      ggplot(aes(x = bin + 0.5 * bar_w, y = prop,
                 fill = method, color = method, alpha = I(bar_alpha))) +
      geom_col(width = bar_w, linewidth = 0.15) +
      geom_vline(xintercept = EXCL_SEP, linetype = "dotted",
                 color = "gray50", linewidth = 0.5) +
      facet_grid(method ~ ontology) +
      scale_fill_manual(values  = c(CAV = CAV_COLOR, DeepGoSE = TOOL_COLOR), name = NULL) +
      scale_color_manual(values = c(CAV = CAV_COLOR, DeepGoSE = TOOL_COLOR), name = NULL) +
      scale_x_continuous(
        limits = c(0, EXCL_X + bar_w * 2),
        breaks = c(0, 0.25, 0.50, 0.75, 1.00, EXCL_X),
        labels = c("0%", "25%", "50%", "75%", "100%", "Excl."),
        expand = c(0, 0)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      base_theme() +
      panel_border() +
      labs(x = "Rank percentile  (0% = true GO term ranked first)", y = "Proportion") +
      theme(
        legend.position  = "none",
        strip.background = element_blank(),
        strip.text.y     = element_text(angle = 0, hjust = 0)
      )
  } else {
    message("Skipping rank histogram: no go_specificity_ranks_*.csv files found")
    p_2a <- NULL
  }

  # ---------------------------------------------------------------------------
  # Panel B: ggridges for AUC — one row per ontology, both methods overlaid
  # ---------------------------------------------------------------------------
  p_2b <- ont_comp |>
    select(ontology, go_term,
           CAV      = auc_val_vs_test_neg,
           DeepGoSE = tool_auc) |>
    pivot_longer(c(CAV, DeepGoSE), names_to = "method", values_to = "auc") |>
    drop_na(auc) |>
    mutate(ontology = factor(ontology, levels = ONT_LEVELS)) |>
    ggplot(aes(x = auc, y = ontology, fill = method, color = method)) +
    do.call(geom_density_ridges, ridges_opts) +
    scale_fill_manual(values  = c(CAV = CAV_COLOR, DeepGoSE = TOOL_COLOR),
                      name = NULL) +
    scale_color_manual(values = c(CAV = CAV_COLOR, DeepGoSE = TOOL_COLOR),
                       name = NULL) +
    scale_x_continuous(limits = c(0.4, 1.0), expand = c(0, 0)) +
    base_theme() +
    labs(x = "AUC", y = NULL) +
    theme(legend.position = "none")

  # ---------------------------------------------------------------------------
  # Panel C: ggridges for AUPR — same layout as B
  # ---------------------------------------------------------------------------
  p_2c <- ont_comp |>
    select(ontology, go_term,
           CAV      = aupr_val_vs_test_neg,
           DeepGoSE = tool_aupr) |>
    pivot_longer(c(CAV, DeepGoSE), names_to = "method", values_to = "aupr") |>
    drop_na(aupr) |>
    mutate(ontology = factor(ontology, levels = ONT_LEVELS)) |>
    ggplot(aes(x = aupr, y = ontology, fill = method, color = method)) +
    do.call(geom_density_ridges, ridges_opts) +
    scale_fill_manual(values  = c(CAV = CAV_COLOR, DeepGoSE = TOOL_COLOR),
                      name = NULL) +
    scale_color_manual(values = c(CAV = CAV_COLOR, DeepGoSE = TOOL_COLOR),
                       name = NULL) +
    scale_x_continuous(limits = c(0, 1.0), expand = c(0, 0)) +
    base_theme() +
    labs(x = "AUPR", y = NULL) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )

  # ---------------------------------------------------------------------------
  # Assemble Figure 2
  # A on top (full width); B and C side by side on bottom
  # Legend extracted from B (always present) and placed right of B+C
  # ---------------------------------------------------------------------------
  legend_fig2 <- get_legend(
    p_2b + theme(legend.position = "right")
  )

  bottom_row <- plot_grid(
    p_2b, p_2c, legend_fig2,
    nrow       = 1,
    labels     = c("B", "C", ""),
    label_size = 8,
    rel_widths = c(1, 1, 0.2)
  )

  if (!is.null(p_2a)) {
    p_2a_noleg <- p_2a + theme(legend.position = "none")
    fig2 <- plot_grid(
      p_2a_noleg,
      bottom_row,
      ncol        = 1,
      labels      = c("A", ""),
      label_size  = 8,
      rel_heights = c(1, 1.3)
    )
  } else {
    fig2 <- plot_grid(bottom_row, ncol = 1, labels = c("B"), label_size = 8)
  }

  ggsave(file.path(OUT, "fig2.pdf"), fig2, width = 4.7, height = 3.7)
  message("Saved fig2.pdf")

} else {
  message("Skipping Figure 2: no temporal_tool_comparison_*.csv files found")
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

message("\nAll figures written to ", OUT, "/")
