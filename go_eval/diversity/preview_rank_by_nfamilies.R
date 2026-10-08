# Exploratory preview (2026-09-29): Fig 2C-style rank composition, binned by the
# absolute number of sequence families among each GO term's CAV training
# positives (term_diversity.py: DIAMOND, >=30% identity over >=50% of the shorter
# sequence, connected components). Known-positive-excluded data; gated ranks as
# in panel C. Bar labels: terms / pairs. Run from bioCAV-paper/go_eval/diversity/.
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

STAGE <- "../../figures/figure_data_staging/noknown"
FIG_FONT <- "Arial"
RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
RANK_COLORS <- c("Rank 1" = "#005A32", "Rank 2-3" = "#238B45", "Rank 4-10" = "#74C476",
                 "Rank >10" = "#C4B7DC", "No prediction" = "#8073AC")
FAM_BREAKS <- c(0, 5, 10, 15, 30, 60, Inf)
FAM_LEVELS <- c("1-4", "5-9", "10-14", "15-29", "30-59", ">=60")
ONTS <- c("mf", "bp", "cc")
tier <- function(rank, nopred) {
  case_when(nopred ~ "No prediction", rank == 1 ~ "Rank 1", rank <= 3 ~ "Rank 2-3",
            rank <= 10 ~ "Rank 4-10", TRUE ~ "Rank >10")
}

pairs <- map_dfr(ONTS, \(o) {
  read_csv(sprintf("%s/go_specificity_ranks_allterms_%s.csv", STAGE, o), show_col_types = FALSE) |>
    inner_join(read_csv(sprintf("term_diversity_%s.csv", o), show_col_types = FALSE),
               by = "go_term") |>
    mutate(ontology = toupper(o))
}) |>
  mutate(fam_bin = cut(n_families, FAM_BREAKS, right = FALSE, labels = FAM_LEVELS))
pairs <- bind_rows(pairs, pairs |> mutate(ontology = "All")) |>
  mutate(ontology = factor(ontology, levels = c(toupper(ONTS), "All")))

long <- bind_rows(
  pairs |> transmute(ontology, go_term, fam_bin, method = "CAV",
                     tier = tier(cav_rank_ont_filt, llr <= 0)),
  pairs |> transmute(ontology, go_term, fam_bin, method = "DeepGoSE",
                     tier = tier(tool_rank_ont_filt, !as.logical(tool_predicted)))
) |> mutate(tier = factor(tier, levels = rev(RANK_LEVELS)),
            method = factor(method, levels = c("CAV", "DeepGoSE")))

bin_n <- pairs |> group_by(ontology, fam_bin) |>
  summarise(lab = sprintf("%d\n%d", n_distinct(go_term), n()), .groups = "drop")

fig <- ggplot(long, aes(fam_bin, fill = tier)) +
  geom_bar(position = "fill", width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(data = bin_n, aes(fam_bin, 1.09, label = lab), inherit.aes = FALSE,
            size = 2, lineheight = 0.85, family = FIG_FONT, color = "gray30") +
  facet_grid(method ~ ontology) +
  scale_fill_manual(values = RANK_COLORS, breaks = RANK_LEVELS, name = NULL) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1),
                     limits = c(0, 1.16), expand = c(0, 0)) +
  theme_cowplot(font_size = 8, font_family = FIG_FONT) +
  theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8),
        strip.text = element_text(size = 8), strip.background = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top") +
  labs(x = "Sequence families among the term's CAV training positives",
       y = "Share of validation protein-GO pairs",
       caption = "Bar labels: GO terms / pairs")

ggsave("preview_rank_by_nfamilies.png", fig, width = 7.2, height = 4.2, dpi = 300, bg = "white")
ggsave("preview_rank_by_nfamilies.pdf", fig, width = 7.2, height = 4.2, device = cairo_pdf)

long |> group_by(ontology, fam_bin, method) |>
  summarise(terms = n_distinct(go_term), pairs = n(),
            rank1 = round(100 * mean(tier == "Rank 1"), 1),
            top10 = round(100 * mean(tier %in% c("Rank 1", "Rank 2-3", "Rank 4-10")), 1),
            nopred = round(100 * mean(tier == "No prediction"), 1), .groups = "drop") |>
  filter(ontology == "All") |> print(n = 20)
