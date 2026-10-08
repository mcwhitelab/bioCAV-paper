#!/usr/bin/env Rscript
# draft_fig2C_allterms.R — standalone alternative to Fig 2C (rank-tier composition)
#
# Shipped Fig 2C ranks each (validation protein, true GO term) pair only among the
# terms evaluated in that ontology (MF 226 / BP 395 / CC 258). This draft ranks the
# same pairs against wider candidate sets; data from
# go_eval/rank_all_terms/compute_ranks_all_terms.py ->
# figure_data/go_specificity_ranks_allterms_{mf,bp,cc}.csv.
#
# "No prediction" gates are unchanged from the shipped panel (CAV: LLR <= 0;
# DeepGo-SE: term not emitted). Not invoked by 00_make_figures.sh.
#
# Outputs (figures/):
#   draft_fig2C_allterms.{pdf,png}             — all 6,677 GO CAVs, filtered rank,
#                                                drop-in size for panel C
#   draft_fig2C_allterms_variants.{pdf,png}    — side-by-side of every candidate set

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

DATA <- "figure_data"
OUT  <- "figures"
FIG_FONT <- "Arial"

base_theme <- function(...) {
  theme_cowplot(font_size = 8, font_family = FIG_FONT, ...) +
    theme(
      axis.text   = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text  = element_text(size = 8)
    )
}

METHOD_LEVELS <- c("CAV", "DeepGoSE")
ONT_DISPLAY   <- c("MF", "BP", "CC")
RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
RANK_COLORS <- c(
  "Rank 1"        = "#005A32",
  "Rank 2-3"      = "#238B45",
  "Rank 4-10"     = "#74C476",
  "Rank >10"      = "#C4B7DC",
  "No prediction" = "#8073AC"
)

tier_of <- function(rank, predicted) {
  case_when(
    !predicted ~ "No prediction",
    rank == 1  ~ "Rank 1",
    rank <= 3  ~ "Rank 2-3",
    rank <= 10 ~ "Rank 4-10",
    TRUE       ~ "Rank >10"
  )
}

ranks <- map_dfr(c("mf", "bp", "cc"), \(o)
  read_csv(file.path(DATA, paste0("go_specificity_ranks_allterms_", o, ".csv")),
           show_col_types = FALSE) |> mutate(ontology = toupper(o)))

VARIANTS <- c(
  eval_shipped = "Shipped: evaluated terms\n(226 / 395 / 258)",
  ont          = "Same-ontology CAVs\n(1,161 / 4,687 / 707)",
  ont_filt     = "Same-ontology CAVs,\nfiltered",
  all          = "All 6,677 GO CAVs",
  all_filt     = "All 6,677 GO CAVs,\nfiltered"
)

tiers <- map_dfr(names(VARIANTS), \(v)
  ranks |>
    transmute(
      ontology, variant = v,
      CAV      = tier_of(.data[[paste0("cav_rank_",  v)]], llr > 0),
      DeepGoSE = tier_of(.data[[paste0("tool_rank_", v)]], tool_predicted == TRUE)
    )) |>
  pivot_longer(c(CAV, DeepGoSE), names_to = "method", values_to = "tier") |>
  count(variant, ontology, method, tier) |>
  group_by(variant, ontology, method) |>
  mutate(prop = n / sum(n)) |>
  ungroup() |>
  mutate(
    tier     = factor(tier,     levels = RANK_LEVELS),
    method   = factor(method,   levels = METHOD_LEVELS),
    ontology = factor(ontology, levels = rev(ONT_DISPLAY)),
    variant  = factor(variant,  levels = names(VARIANTS), labels = VARIANTS)
  )

rank_bars <- function(d, legend_rows = 1) {
  ggplot(d, aes(x = prop, y = ontology, fill = tier)) +
    geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = RANK_COLORS, breaks = RANK_LEVELS, name = NULL,
                      drop = FALSE) +
    scale_x_continuous(labels = scales::percent,
                       expand = expansion(mult = c(0, 0.02))) +
    base_theme() +
    labs(x = "Share of validation protein-GO pairs", y = NULL) +
    guides(fill = guide_legend(nrow = legend_rows)) +
    theme(
      legend.position    = "top",
      legend.location    = "plot",
      legend.justification = "left",
      legend.key.size    = unit(6, "pt"),
      legend.margin      = margin(t = 0, b = 2, r = 0, l = 16),
      legend.spacing.x   = unit(0, "pt"),
      legend.box.spacing = unit(2, "pt"),
      strip.background   = element_blank(),
      strip.placement    = "outside",
      strip.text.y.left  = element_text(angle = 0, size = 8, margin = margin(l = 4, r = 6)),
      axis.text.y        = element_text(margin = margin(r = 4)),
      panel.spacing.y    = unit(10, "pt"),
      plot.margin        = margin(t = 2, r = 10, b = 2, l = 2)
    )
}

# Drop-in panel: all 6,677 CAVs, filtered rank
p_main <- tiers |>
  filter(variant == VARIANTS[["all_filt"]]) |>
  rank_bars(legend_rows = 2) +
  facet_grid(method ~ ., switch = "y")

ggsave(file.path(OUT, "draft_fig2C_allterms.pdf"), p_main,
       width = 3.6, height = 2.05, device = cairo_pdf)
ggsave(file.path(OUT, "draft_fig2C_allterms.png"), p_main,
       width = 3.6, height = 2.05, dpi = 300, bg = "white")

# Comparison sheet: every candidate set side by side
p_var <- rank_bars(tiers) +
  facet_grid(method ~ variant, switch = "y") +
  theme(strip.text.x = element_text(size = 7, lineheight = 0.9),
        panel.spacing.x = unit(8, "pt"),
        axis.text.x = element_text(size = 6.5))

ggsave(file.path(OUT, "draft_fig2C_allterms_variants.pdf"), p_var,
       width = 9, height = 2.4, device = cairo_pdf)
ggsave(file.path(OUT, "draft_fig2C_allterms_variants.png"), p_var,
       width = 9, height = 2.4, dpi = 300, bg = "white")
