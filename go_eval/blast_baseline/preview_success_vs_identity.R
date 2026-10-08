# Exploratory preview (2026-09-29): Fig 2C success rate vs homology to the
# true term's training positives. x = best DIAMOND identity of the validation
# protein to any CAV training positive of its true term (blast_ranks.py,
# self-hits removed; "No hit" = nothing at e <= 1e-3). y = share of pairs with
# the true term at Rank 1 (top row) or in the top 10 (bottom row), gated as in
# preview_ranks_pf.R (CAV llr > 0, DeepGo-SE predicted, BLAST hit). BLAST solid
# = tie-averaged rank, dotted = optimistic rank. Candidates exclude obsolete GO
# terms (noobsolete/, from run_noobsolete.py). Error bars: Wilson 95% CI.
# Not part of figures.R. Run from bioCAV-paper/go_eval/blast_baseline/.
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

DATA <- "noobsolete"
FIG_FONT <- "Arial"
ONTS <- c("mf", "bp", "cc")
ID_BREAKS <- c(-Inf, 30, 40, 50, 70, 90, Inf)
ID_LEVELS <- c("No hit", "<30", "30-40", "40-50", "50-70", "70-90", ">=90")
# Okabe-Ito, matching figures.R (CAV blue, DeepGo-SE vermillion); BLAST neutral.
METHOD_COLORS <- c("CAV" = "#0072B2", "DeepGo-SE" = "#D55E00", "BLAST" = "gray35")
METHOD_SHAPES <- c("CAV" = 16, "DeepGo-SE" = 17, "BLAST" = 15)

pairs <- map_dfr(ONTS, \(o) {
  # one (protein, term) pair per ontology appears twice as identical rows; keep one
  read_csv(sprintf("%s/go_specificity_ranks_allterms_%s.csv", DATA, o), show_col_types = FALSE) |>
    distinct() |>
    inner_join(read_csv(sprintf("%s/blast_ranks_%s.csv", DATA, o), show_col_types = FALSE) |> distinct(),
               by = c("protein_id", "go_term")) |>
    mutate(ontology = toupper(o))
}) |>
  mutate(id_bin = if_else(blast_predicted,
                          as.character(cut(blast_true_pident, ID_BREAKS, right = FALSE,
                                           labels = ID_LEVELS[-1])),
                          "No hit") |> factor(levels = ID_LEVELS),
         ontology = factor(ontology, levels = c("MF", "BP", "CC")))

long <- pairs |>
  transmute(ontology, id_bin,
            `CAV|1`        = llr > 0 & cav_rank_ont_filt_pf < 1.5,
            `CAV|10`       = llr > 0 & cav_rank_ont_filt_pf <= 10,
            `DeepGo-SE|1`  = tool_predicted == TRUE & tool_rank_ont_filt_pf < 1.5,
            `DeepGo-SE|10` = tool_predicted == TRUE & tool_rank_ont_filt_pf <= 10,
            `BLAST|1`      = blast_predicted & blast_rank_mid < 1.5,
            `BLAST|10`     = blast_predicted & blast_rank_mid <= 10,
            `BLAST opt|1`  = blast_predicted & blast_rank_opt < 1.5,
            `BLAST opt|10` = blast_predicted & blast_rank_opt <= 10) |>
  pivot_longer(-c(ontology, id_bin), names_to = c("method", "k"), names_sep = "\\|",
               values_to = "hit")

wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n
  c <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  tibble(lo = pmax(0, c - h), hi = pmin(1, c + h))
}

rates <- long |>
  group_by(ontology, id_bin, method, k) |>
  summarise(x = sum(hit), n = n(), .groups = "drop") |>
  mutate(rate = x / n, wilson(x, n),
         opt = method == "BLAST opt",
         method = factor(if_else(opt, "BLAST", method), levels = names(METHOD_COLORS)),
         k = factor(if_else(k == "1", "Rank 1", "Top 10"), levels = c("Rank 1", "Top 10")))
write_csv(rates, "preview_success_vs_identity.csv")

n_lab <- pairs |> count(ontology, id_bin)
dodge <- position_dodge(width = 0.5)
main <- filter(rates, !opt)

p <- ggplot(main, aes(x = id_bin, y = rate, color = method, group = method)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "gray94") +
  geom_line(data = filter(main, id_bin != "No hit"), position = dodge, linewidth = 0.6) +
  geom_line(data = filter(rates, opt, id_bin != "No hit"), position = dodge,
            linewidth = 0.5, linetype = "11") +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = dodge, width = 0, linewidth = 0.4) +
  geom_point(aes(shape = method), position = dodge, size = 1.4) +
  geom_text(data = n_lab, aes(x = id_bin, y = Inf, label = n), inherit.aes = FALSE,
            vjust = 1.4, size = 1.9, family = FIG_FONT, color = "gray40") +
  facet_grid(k ~ ontology, scales = "free_y", switch = "y") +
  scale_x_discrete(limits = ID_LEVELS) +
  scale_color_manual(values = METHOD_COLORS, name = NULL) +
  scale_shape_manual(values = METHOD_SHAPES, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.02, 0.14))) +
  theme_cowplot(font_size = 8, font_family = FIG_FONT) +
  labs(x = "Best identity (%) to the true term's training positives (DIAMOND)", y = NULL,
       caption = "Numbers: pairs per bin. BLAST dotted: optimistic tie rank. Bars: Wilson 95% CI.") +
  theme(axis.text = element_text(size = 7), axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 8), strip.background = element_blank(),
        strip.placement = "outside", strip.text = element_text(size = 8, face = "bold"),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
        plot.caption = element_text(size = 6.5, color = "gray40", hjust = 0),
        panel.spacing = unit(8, "pt"))
ggsave("preview_success_vs_identity.png", p, width = 7.2, height = 4.2, dpi = 300, bg = "white")
ggsave("preview_success_vs_identity.pdf", p, width = 7.2, height = 4.2, device = cairo_pdf)

print(main |> mutate(pct = round(100 * rate, 1)) |>
        select(ontology, k, id_bin, method, pct) |>
        pivot_wider(names_from = method, values_from = pct) |>
        arrange(k, ontology, id_bin), n = 50)
