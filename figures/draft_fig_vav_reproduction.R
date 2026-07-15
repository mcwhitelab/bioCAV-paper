#!/usr/bin/env Rscript
# DRAFT / exploratory -- reproduction of arxiv 2511.21614v1 Figure 2
# (Q9NHV9/VAV_DROME layerwise CAV score profiles across PF00621/PF00130/
# PF00017/PF00018), using freshly-trained 36-layer CAVs (see
# figure_data/vav_motif_repro/train_all_layers.log and
# score_vav_all_layers.py) plus the already-existing layer-25 CAVs.
# Slated to become Figure 4 (fills the fig2/fig3/[gap]/fig5 numbering gap
# in figures.R) once validated -- kept standalone for now, same
# draft-then-promote pattern used for the Fig 5 split-strip panels (see
# memory fig5-continuum-split-strips). Not yet part of the figures.R
# pipeline; run standalone with:
#   Rscript draft_fig_vav_reproduction.R
#
# Layer-label convention: the CAV filename number L{k} indexes
# hidden_states[k] directly (k=0 is the raw embedding output, k=1..36 are
# the 36 transformer layers). The paper's "indexed from 1" caption counts
# the embedding output as layer 1, so display_layer = k + 1 -- this makes
# the pre-existing L25 CAV (used for the paper's panel B) display as
# "Layer 26", matching the paper exactly. score_vav_all_layers.py already
# applies this offset when writing vav_layerwise_scores.json.

source("draft_fig_extra_proteins_common.R")

DATA <- "figure_data/vav_motif_repro"
OUT  <- "figures"

base_theme <- function(...) theme_cowplot(font_size = 8, ...)

j <- fromJSON(file.path(DATA, "vav_layerwise_scores.json"), simplifyVector = FALSE)
seq_len <- j$seq_len

DOMAIN_COLORS <- c(PF00621 = "#B83280", PF00130 = "#2E7D32",
                   PF00017 = "#C0392B", PF00018 = "#1F4E96")
DOMAIN_LABELS <- c(PF00621 = "PF00621 (RhoGEF Domain)", PF00130 = "PF00130 (C1 Domain)",
                   PF00017 = "PF00017 (SH2 Domain)",    PF00018 = "PF00018 (SH3 Domain)")
DOMAIN_ORDER  <- c("PF00621", "PF00130", "PF00017", "PF00018")
HIGHLIGHT_LAYER <- 26L  # display layer (1-indexed, embedding counted as layer 1)

# --- long-format curve data: domain, layer, position, score ---
curve_rows <- list()
gt_rows <- list()
for (motif_id in DOMAIN_ORDER) {
  dom <- j$domains[[motif_id]]
  if (is.null(dom)) next
  gt <- j$ground_truth[[motif_id]]
  gt_rows[[motif_id]] <- tibble(domain = motif_id, start = gt$start, end = gt$end, name = gt$name)

  plc <- dom$per_layer_curve
  for (layer_str in names(plc)) {
    scores <- unlist(plc[[layer_str]])
    curve_rows[[paste(motif_id, layer_str)]] <- tibble(
      domain = motif_id, layer = as.integer(layer_str),
      position = seq_along(scores), score = scores
    )
  }
}
curves <- bind_rows(curve_rows) %>% filter(!is.na(score))
gt_df  <- bind_rows(gt_rows)

n_layers <- max(curves$layer)
message(sprintf("Loaded %d domains, %d layers, %d residues", length(unique(curves$domain)), n_layers, seq_len))

# --- domain track (shared by both panels) ---
build_domain_track <- function() {
  ggplot() +
    geom_rect(aes(xmin = 0, xmax = seq_len, ymin = 0, ymax = 1), fill = "#e8c087") +
    geom_rect(data = gt_df, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = domain)) +
    geom_text(data = gt_df, aes(x = (start + end) / 2, y = -0.6, label = domain, color = domain),
              size = 2.3, fontface = "bold") +
    scale_fill_manual(values = DOMAIN_COLORS, guide = "none") +
    scale_color_manual(values = DOMAIN_COLORS, guide = "none") +
    scale_x_continuous(limits = c(0, seq_len), expand = c(0, 0),
                        breaks = c(0, seq_len)) +
    coord_cartesian(ylim = c(-1.4, 1), clip = "off") +
    theme_void() +
    theme(axis.text.x = element_text(size = 6),
          plot.margin = margin(t = 0, b = 10, l = 2, r = 2))
}

# --- Panel A: all layers, light (early) -> dark (late), bold highlight ---
p_a <- ggplot() +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_line(data = curves %>% filter(layer != HIGHLIGHT_LAYER),
            aes(position, score, group = interaction(domain, layer), color = domain,
                alpha = layer),
            linewidth = 0.25) +
  geom_line(data = curves %>% filter(layer == HIGHLIGHT_LAYER),
            aes(position, score, group = domain, color = domain),
            linewidth = 0.9) +
  scale_color_manual(values = DOMAIN_COLORS, labels = DOMAIN_LABELS, name = NULL) +
  scale_alpha_continuous(range = c(0.08, 0.55), guide = "none") +
  scale_x_continuous(limits = c(0, seq_len), expand = c(0, 0)) +
  base_theme() +
  labs(x = NULL, y = "CAV Score",
       title = sprintf("Q9NHV9 - VAV_DROME  (all %d layers, bold = layer %d)", n_layers, HIGHLIGHT_LAYER)) +
  theme(legend.position = "right", legend.text = element_text(size = 6.5),
        legend.key.size = unit(9, "pt"), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 8, face = "bold"))

panel_a <- plot_grid(p_a, build_domain_track(), ncol = 1, align = "v", axis = "lr",
                     rel_heights = c(1, 0.16))

# --- Panel B: layer 26 only, with dashed lines from each domain's peak
# position down to its ground-truth interval ---
curves_hl <- curves %>% filter(layer == HIGHLIGHT_LAYER)
peak_df <- curves_hl %>% group_by(domain) %>% slice_max(score, n = 1, with_ties = FALSE) %>% ungroup()

p_b <- ggplot() +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_segment(data = peak_df, aes(x = position, xend = position, y = score, yend = -16),
               color = "grey50", linewidth = 0.3, linetype = "dashed") +
  geom_line(data = curves_hl, aes(position, score, color = domain), linewidth = 0.9) +
  scale_color_manual(values = DOMAIN_COLORS, labels = DOMAIN_LABELS, name = NULL) +
  scale_x_continuous(limits = c(0, seq_len), expand = c(0, 0)) +
  coord_cartesian(ylim = c(min(curves_hl$score, na.rm = TRUE) - 1,
                           max(curves_hl$score, na.rm = TRUE) + 1), clip = "off") +
  base_theme() +
  labs(x = NULL, y = sprintf("CAV Score (Layer %d)", HIGHLIGHT_LAYER),
       title = "Q9NHV9 - VAV_DROME") +
  theme(legend.position = "right", legend.text = element_text(size = 6.5),
        legend.key.size = unit(9, "pt"), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 8, face = "bold"))

panel_b <- plot_grid(p_b, build_domain_track(), ncol = 1, align = "v", axis = "lr",
                     rel_heights = c(1, 0.16))

# --- Panel C: two more example proteins showing repeated-domain
# localization (Q9Y7R4/SET1_SCHPO: PF00076 RRM has an unannotated second
# occurrence; Q54SA1/PLDZ_DICDI: PF00614 similarly). Layer 26 only, from
# the 20k-library L25 CAVs (see score_extra_proteins.py) -- these organically
# reproduce the paper's "second peak missed by UniProt/InterPro-N" finding
# since the extra peak just falls out of real find_peaks() calls, not a
# hardcoded position. ---
extra_scores <- fromJSON(file.path(DATA, "extra_proteins", "extra_proteins_scores.json"), simplifyVector = FALSE)
p_set1  <- build_extra_protein_panel("Q9Y7R4", extra_scores[["Q9Y7R4"]], title = "Q9Y7R4 - SET1_SCHPO", base_theme_fn = base_theme)
p_pldz  <- build_extra_protein_panel("Q54SA1", extra_scores[["Q54SA1"]], title = "Q54SA1 - PLDZ_DICDI", base_theme_fn = base_theme)
panel_c <- plot_grid(p_set1, p_pldz, nrow = 1, labels = c("C", "D"), label_size = 10)

fig <- plot_grid(panel_a, panel_b, panel_c, ncol = 1, labels = c("A", "B", NULL), label_size = 10,
                 rel_heights = c(1, 1, 0.75))

ggsave(file.path(OUT, "draft_fig4_vav_reproduction.pdf"), fig,
       width = 7.5, height = 9, bg = "white", device = cairo_pdf)
ggsave(file.path(OUT, "draft_fig4_vav_reproduction.png"), fig,
       width = 7.5, height = 9, bg = "white", dpi = 220)
message("Saved draft_fig4_vav_reproduction.pdf/.png")
