#!/usr/bin/env Rscript
# Shared plotting helpers for the two draft reproductions built from
# extra_proteins_scores.json (arxiv 2511.21614v1 Appendix A grid +
# Figure 3's repeat-domain showcase proteins). Sourced by
# draft_fig_supp_appendixA.R and draft_fig_vav_reproduction.R (for the
# SET1_SCHPO / PLDZ_DICDI panels). Not part of the figures.R pipeline.

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(jsonlite)
})

# figures.R sets FIG_FONT globally before sourcing this file; the standalone
# draft_* scripts may not, so fall back to grid's default family.
EXTRA_FONT <- function() if (exists("FIG_FONT")) get("FIG_FONT") else ""

EXTRA_DOMAIN_PALETTE <- c("#C0392B", "#1F4E96", "#2E7D32", "#B8860B", "#7B3F9E",
                          "#008080", "#D2691E", "#5D3FD3", "#8B8B00", "#C71585")

# Build one protein's panel: curve(s) + domain track below, with dashed
# lines from each peak down to the track. Peaks with annotated == FALSE
# (the "second, unannotated occurrence" cases) get a hatched/outline-only
# track segment instead of a solid fill, and a lighter dashed peak line.
# With parts = TRUE the two grobs are returned unstacked, so a caller can
# interleave them into a larger plot_grid and get real cowplot alignment
# against other panels (a nested plot_grid cannot be aligned).
# domain_colors / domain_order are optional per-call overrides. The default
# is the flat EXTRA_DOMAIN_PALETTE assigned in sequence order, which treats
# every Pfam family as unrelated; pass a named color vector when some of the
# domains belong to one family and should read as a ramp of a single hue
# (composite encoding), and a domain_order vector to group those slots
# together in the legend. Without domain_order the legend follows ggplot's
# discrete-scale default (alphabetical by Pfam accession), which is not the
# same as the sequence order used to assign colors.
build_extra_protein_panel <- function(acc, protein_data, title = acc, base_theme_fn, legend_position = "right",
                                      parts = FALSE, domain_colors = NULL, domain_order = NULL,
                                      legend_text_size = 6, axis_text_size = 6,
                                      y_lab = "CAV Score (Layer 26)") {
  seq_len <- protein_data$seq_len

  # order domains (and thus their color/legend order) by earliest
  # sequence position they appear at, not JSON insertion order
  domain_ids_raw <- names(protein_data$domains)
  domain_first_pos <- map_dbl(domain_ids_raw, function(m) {
    gt <- protein_data$domains[[m]]$ground_truth
    if (length(gt) > 0) min(map_dbl(gt, ~ .x$start)) else Inf
  })
  domain_ids <- domain_ids_raw[order(domain_first_pos)]

  if (is.null(domain_colors)) {
    domain_colors <- setNames(EXTRA_DOMAIN_PALETTE[seq_along(domain_ids)], domain_ids)
  } else {
    missing_cols <- setdiff(domain_ids, names(domain_colors))
    if (length(missing_cols) > 0) {
      stop("build_extra_protein_panel(", acc, "): domain_colors is missing ",
           paste(missing_cols, collapse = ", "))
    }
    domain_colors <- domain_colors[domain_ids]
  }

  # legend shows "PFxxxxx (Pfam family name)" -- name comes from the first
  # ground-truth annotation for that domain (same across repeated instances)
  domain_labels <- setNames(map_chr(domain_ids, function(m) {
    gt <- protein_data$domains[[m]]$ground_truth
    nm <- if (length(gt) > 0) gt[[1]]$name else NULL
    if (is.null(nm) || !nzchar(nm)) m else sprintf("%s (%s)", m, nm)
  }), domain_ids)

  curve_rows <- list()
  peak_rows <- list()
  track_rows <- list()
  for (motif_id in domain_ids) {
    dom <- protein_data$domains[[motif_id]]
    scores <- unlist(dom$curve)
    curve_rows[[motif_id]] <- tibble(domain = motif_id, position = seq_along(scores), score = scores)

    peaks <- dom$peaks
    if (length(peaks) > 0) {
      peak_rows[[motif_id]] <- map_dfr(peaks, function(p) {
        tibble(domain = motif_id, position = p$position, score = p$score, annotated = isTRUE(p$annotated))
      })
    }

    gt <- dom$ground_truth
    if (length(gt) > 0) {
      track_rows[[paste0(motif_id, "_gt")]] <- map_dfr(gt, function(g) {
        tibble(domain = motif_id, start = g$start, end = g$end, kind = "annotated")
      })
    }
    # candidate (unannotated) peaks get their own track segment, width = window_length
    if (length(peaks) > 0) {
      cand <- keep(peaks, ~ !isTRUE(.x$annotated))
      if (length(cand) > 0) {
        half_w <- dom$window_length %/% 2
        track_rows[[paste0(motif_id, "_cand")]] <- map_dfr(cand, function(p) {
          tibble(domain = motif_id, start = p$position - half_w, end = p$position + half_w, kind = "candidate")
        })
      }
    }
  }
  curves <- bind_rows(curve_rows) %>% filter(!is.na(score))
  peaks_df <- bind_rows(peak_rows)
  track_df <- bind_rows(track_rows)

  # legend order: ggplot sorts a character mapping alphabetically, so an
  # explicit order has to be baked in as factor levels on every frame that
  # carries `domain`.
  if (!is.null(domain_order)) {
    lv <- c(intersect(domain_order, domain_ids), setdiff(domain_ids, domain_order))
    curves$domain <- factor(curves$domain, levels = lv)
    if (nrow(peaks_df) > 0) peaks_df$domain <- factor(peaks_df$domain, levels = lv)
    if (nrow(track_df) > 0) track_df$domain <- factor(track_df$domain, levels = lv)
  }

  y_lo <- min(curves$score, na.rm = TRUE)
  y_hi <- max(curves$score, na.rm = TRUE)

  p_curve <- ggplot() +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3)
  if (nrow(peaks_df) > 0) {
    p_curve <- p_curve +
      geom_segment(data = peaks_df,
                   aes(x = position, xend = position, y = score, yend = y_lo - 1, color = domain,
                       linetype = annotated),
                   linewidth = 0.3, alpha = 0.7, show.legend = FALSE)
  }
  p_curve <- p_curve +
    geom_line(data = curves, aes(position, score, color = domain), linewidth = 0.7) +
    scale_color_manual(values = domain_colors, labels = domain_labels, name = NULL) +
    scale_linetype_manual(values = c(`TRUE` = "22", `FALSE` = "13")) +
    scale_x_continuous(limits = c(0, seq_len), expand = c(0, 0)) +
    coord_cartesian(ylim = c(y_lo - 1, y_hi + 1), clip = "off") +
    base_theme_fn() +
    labs(x = NULL, y = y_lab, title = title) +
    guides(color = if (legend_position == "bottom") guide_legend(nrow = 2) else "legend") +
    theme(legend.position = legend_position, legend.text = element_text(size = legend_text_size),
          legend.key.size = unit(8, "pt"), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 8, face = "bold"))

  p_track <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = seq_len, ymin = 0, ymax = 1), fill = "#e8c087")
  ann <- track_df %>% filter(kind == "annotated")
  cand <- track_df %>% filter(kind == "candidate")
  if (nrow(ann) > 0) {
    p_track <- p_track + geom_rect(data = ann, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = domain))
  }
  if (nrow(cand) > 0) {
    p_track <- p_track + geom_rect(data = cand, aes(xmin = start, xmax = end, ymin = 0, ymax = 1, color = domain),
                                    fill = NA, linewidth = 0.5, linetype = "dashed")
  }
  p_track <- p_track +
    scale_fill_manual(values = domain_colors, guide = "none") +
    scale_color_manual(values = domain_colors, guide = "none") +
    scale_x_continuous(limits = c(0, seq_len), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    theme_void(base_family = EXTRA_FONT()) +
    theme(axis.text.x = element_text(size = axis_text_size),
          plot.margin = margin(t = 0, b = 2, l = 2, r = 2))

  if (parts) return(list(curve = p_curve, track = p_track))
  plot_grid(p_curve, p_track, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 0.14))
}
