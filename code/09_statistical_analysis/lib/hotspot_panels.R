# The trans-hotspot panel: the trans-protein count per locus across the genome
# (01_pqtl_architecture.png c). Plotting only, from the table 01b_hotspot_data.R writes,
# so 01c_plot.R draws it without re-running 01b (which needs the server's VEP
# annotation). Needs ggplot2, lib/common.R (CHR_LEN_HG19, HOTSPOT_MIN) and
# lib/fig_theme.R (FLUID_COL).
#   loci     - 01b_hotspot_loci.csv

hotspot_autosomes <- as.character(1:22)
hotspot_chr_len <- CHR_LEN_HG19[hotspot_autosomes]   # autosomes only here
hotspot_chr_offset <- setNames(c(0, cumsum(hotspot_chr_len[-length(hotspot_chr_len)])),
                               hotspot_autosomes)

# Panel a as a function so the same code serves any figure size: narrower copies need
# smaller hotspot labels with more vertical headroom to avoid overlap.
make_hotspot_map <- function(loci, label_size = 3, headroom = 0.12,
                             base_size = 14, title_size = 13, legend_size = 11,
                             short_labels = FALSE, y_title = "N of trans-target proteins",
                             near_frac = 0.04, pair_push = 0.3) {
  to_linear <- function(chr, pos) hotspot_chr_offset[as.character(chr)] + as.numeric(pos) / 1e6
  autosome_end <- sum(hotspot_chr_len)

  plot_d <- copy(loci)
  plot_d[is.na(regulator), regulator := ""]
  plot_d[, x := to_linear(chr, lead_pos)]
  # Long form: one row per locus per fluid. The hotspot definition stays on the pooled
  # target count (a locus is one locus in both fluids); the facets then show how that
  # locus's targets are distributed, which is where the CSF/plasma split becomes visible.
  plot_long <- melt(
    plot_d, id.vars = c("x", "chr", "n_proteins"),
    measure.vars = c("n_csf", "n_plasma"),
    variable.name = "fluid", value.name = "n_fluid"
  )
  plot_long[, fluid := factor(fifelse(fluid == "n_csf", "CSF", "Plasma"),
                              levels = c("CSF", "Plasma"))]
  # Re-derived per facet rather than carried over from the locus: this is the same rule
  # as loci$is_hotspot_csf / is_hotspot_plasma, applied to whichever fluid the row is.
  # A locus with no target in this fluid is not drawn in its facet: a row of points at
  # zero only repeats the other facet's loci.
  plot_long <- plot_long[n_fluid > 0]
  plot_long[, is_hotspot := n_fluid >= HOTSPOT_MIN]
  # A hotspot takes its facet's fluid as its class, so it is coloured like panel a.
  plot_long[, class := factor(fifelse(is_hotspot, as.character(fluid), "other"),
                              levels = c("other", "CSF", "Plasma"))]
  # plotmath expressions so only the gene name is italic (chrN: prefix stays upright)
  # Built for every locus: which ones get drawn is decided per facet by is_hotspot.
  plot_d[, lab := fifelse(regulator != "",
                          sprintf("'chr%s: '*italic('%s')", chr, regulator),
                          sprintf("'chr%s'", chr))]
  # Short label (gene only) for the compact copy: the chromosome is already on the
  # x-axis, and dropping the prefix roughly halves the label width, which is what
  # lets the dense chr15-22 hotspots be repelled in y alone without crossing leaders.
  plot_d[, lab_short := fifelse(regulator != "", sprintf("italic('%s')", regulator), lab)]
  n_before <- nrow(plot_long)
  plot_long <- merge(plot_long, plot_d[, .(x, lab, lab_short)], by = "x", all.x = TRUE)
  stopifnot(nrow(plot_long) == n_before)

  x_breaks <- unname(hotspot_chr_offset + hotspot_chr_len / 2)
  odd <- seq(1, length(hotspot_autosomes), 2)
  stripe_df <- data.frame(xmin = hotspot_chr_offset[odd],
                          xmax = hotspot_chr_offset[odd] + hotspot_chr_len[odd])

  hs <- plot_long[is_hotspot == TRUE]
  # Two hotspots this close in one facet (DCBLD2 and KNG1 on chr3) print their names
  # into each other at a centred hjust; pushing the pair apart sideways keeps the
  # fixed-vjust design for every label instead of bringing in a repel layout. A
  # narrower panel widens the labels relative to the genome, so it asks for a larger
  # near_frac (gap, as a fraction of the genome) and pair_push (hjust shift from 0.5).
  near <- near_frac * autosome_end
  setorder(hs, fluid, x)
  hs[, `:=`(gap_next = shift(x, type = "lead") - x, gap_prev = x - shift(x)), by = fluid]
  hs[, hjust := 0.5]
  hs[(gap_next < near) %in% TRUE, hjust := 0.5 + pair_push]
  hs[(gap_prev < near) %in% TRUE, hjust := 0.5 - pair_push]
  ggplot(plot_long, aes(x, n_fluid)) +
    geom_rect(
      data = stripe_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "grey94", inherit.aes = FALSE
    ) +
    # The hotspot threshold, so the coloured points read as the ones above the line.
    geom_hline(yintercept = HOTSPOT_MIN, linetype = "dashed", colour = "grey50",
               linewidth = 0.3) +
    geom_segment(aes(xend = x, yend = 0), colour = "grey70", linewidth = 0.3) +
    # Grey = any other trans locus; a hotspot takes its fluid's colour (blue CSF,
    # orange plasma), the same encoding as panel a, so the two panels read alike.
    # One layer, so the classes can be keyed by a legend the way panel a keys fluid.
    # Most loci act on one or two proteins and sit side by side at the bottom;
    # transparency lets the overlapping ones show as darker points.
    geom_point(aes(fill = class, size = class, alpha = class), shape = 21, stroke = 0.25,
               colour = "grey30") +
    # The sparse hotspot labels do not need collision avoidance. A fixed vjust keeps
    # every name the same visual distance above its point despite the free y scales.
    geom_text(
      data = hs, aes(label = if (short_labels) lab_short else lab, hjust = hjust),
      parse = TRUE, size = label_size, colour = "grey15", vjust = -0.8
    ) +
    # The hotspot key doubles as the definition, so no free-floating annotation has
    # to be fitted between the labelled loci -- earlier versions of that text collided
    # with the loci just above the threshold no matter where it was placed.
    scale_fill_manual(
      values = c(other = "grey55", FLUID_COL),
      labels = c(other = "Trans locus", CSF = "CSF hotspot", Plasma = "Plasma hotspot"),
      name = NULL,
      guide = guide_legend(override.aes = list(size = 2.4))
    ) +
    scale_size_manual(values = c(other = 1.65, CSF = 2.1, Plasma = 2.1), guide = "none") +
    scale_alpha_manual(values = c(other = 0.7, CSF = 1, Plasma = 1), guide = "none") +
    # Free y: on a shared scale the CSF facet is flattened by the plasma ABO locus
    # (75 targets vs 17), and the point of the facets is to read each fluid's own
    # hotspots. The axis numbers carry the between-fluid asymmetry.
    # Strips on the right: on the left they would sit between the y-axis title and
    # the axis numbers.
    # facet_wrap, not facet_grid with space = "free_y": the y here is a count, and
    # sizing the facets by their data range would shrink CSF (0-17) to a sliver next
    # to plasma (0-75), which is exactly what free scales are here to avoid.
    facet_wrap(~ fluid, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_x_continuous(
      breaks = x_breaks, labels = hotspot_autosomes, limits = c(0, autosome_end),
      expand = expansion(mult = 0, add = 0), name = "Chromosome"
    ) +
    scale_y_continuous(
      # Linear on purpose: a square-root axis lifted the 5-9 hotspots off the plasma
      # baseline but was harder to read (tried and dropped 2026-10-29).
      # Zero sits on the lower border, so the stems stand on the axis instead of the
      # CSF facet looking lifted off it; clip = "off" keeps the low-count points whole.
      expand = expansion(mult = c(0, headroom)),
      name = y_title
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
      axis.title = element_text(size = title_size),
      # No box, no bold. The facet label is a quiet index, not a heading competing
      # with the panel tag.
      strip.background = element_blank(),
      strip.text = element_text(size = title_size - 1),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = legend_size),
      legend.key.spacing.x = unit(6, "pt"),
      legend.margin = margin(0, 0, -2, 0),
      plot.margin = margin(6, 8, 6, 6)
      )
}
