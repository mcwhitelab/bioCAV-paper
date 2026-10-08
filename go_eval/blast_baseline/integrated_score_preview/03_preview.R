# Preview: Fig 2 rank panel (live Fig 2B) with two CAV x BLAST integrated scores added.
suppressPackageStartupMessages({ library(tidyverse); library(cowplot) })
here <- "/groups/clairemcwhite/cav_workspace/bioCAV-paper/go_eval/blast_baseline/integrated_score_preview"
FIG_FONT <- "Arial"
RANK_LEVELS <- c("Rank 1", "Rank 2-3", "Rank 4-10", "Rank >10", "No prediction")
RANK_COLORS <- c("Rank 1" = "#005A32", "Rank 2-3" = "#238B45", "Rank 4-10" = "#74C476",
                 "Rank >10" = "#C4B7DC", "No prediction" = "#8073AC")
METHODS <- c("CAV", "DeepGoSE", "BLAST", "CAV x BLAST (rank product)", "CAV x BLAST (soft, eps 0.1)")
LAB <- c("CAV", "DeepGoSE", "BLAST", "CAV×BLAST\nrank product", "CAV×BLAST\nsoft (ε 0.1)")
FAM <- c("<20 families", "20-59 families", ">=60 families")
comp <- read_csv(file.path(here, "integrated_tier_composition.csv"), show_col_types = FALSE) |>
  filter(method %in% METHODS) |>
  mutate(method = factor(method, levels = rev(METHODS), labels = rev(LAB)),
         tier = factor(tier, levels = RANK_LEVELS),
         ontology = factor(ontology, levels = c("MF", "BP", "CC")),
         fam = factor(fam, levels = FAM, labels = c("<20 families", "20-59 families", "≥60 families")))
fam_n <- comp |> filter(method == "CAV") |> group_by(ontology, fam, method) |> summarise(n = sum(n), .groups = "drop")
p <- ggplot(comp, aes(x = prop, y = method, fill = tier)) +
  geom_col(width = 0.72, position = position_stack(reverse = TRUE)) +
  geom_text(data = fam_n, aes(x = 1.01, y = method, label = sprintf("n=%d", n)), inherit.aes = FALSE,
            hjust = 0, size = 2, family = FIG_FONT, color = "gray30") +
  facet_grid(ontology ~ fam) +
  scale_fill_manual(values = RANK_COLORS, breaks = RANK_LEVELS, name = NULL) +
  scale_x_continuous(labels = scales::percent, breaks = c(0, .5, 1), limits = c(0, 1.17), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_cowplot(font_size = 8, font_family = FIG_FONT) +
  theme(axis.text = element_text(size = 7), legend.position = "top", legend.justification = "left",
        legend.key.size = unit(6, "pt"), strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"), panel.spacing.x = unit(10, "pt")) +
  labs(x = "Share of validation protein-GO pairs", y = NULL,
       title = "PREVIEW: Fig 2 rank panel + integrated CAV×BLAST scores",
       caption = "Integrated 'No prediction' = CAV LLR ≤ 0 and no BLAST hit to the true term.")
ggsave(file.path(here, "preview_integrated_ranks.png"), p, width = 7.2, height = 5.6, dpi = 300, bg = "white")
ggsave(file.path(here, "preview_integrated_ranks.pdf"), p, width = 7.2, height = 5.6, device = cairo_pdf)
