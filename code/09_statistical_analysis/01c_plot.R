#!/usr/bin/env Rscript
# The pQTL map, hotspot and replication figures listed below (the study overview is
# the hand-made output/overview_plot/overview.png, used as is). Plotting only: every number here
# comes from the tables written by 01a_pqtl_data.R, 01b_hotspot_data.R and 13, so a
# layout change never re-runs the analysis, and this script runs on the Mac. Re-run 01a / 01b when the underlying data
# change, this script when the figures do.
#   01_pqtl_architecture.png
#     a  credible-set genome map, with the per-chromosome CSF/plasma proportion strip
#     b  plasma Olink--SomaScan concordance (13_somascan_olink_plasma.R): SomaScan match
#        of each Olink lead (top), Olink vs SomaScan effects of the exact leads (bottom)
#     c  trans-pQTL hotspot loci
# The BF2 vs BF1 replication panel (lib/repl_panel.R) is drawn by 09f_plot_novelty.R.
# The hotspot panels are drawn by lib/hotspot_panels.R from 01b's tables.
# Usage:
#   Rscript code/09_statistical_analysis/01c_plot.R [out_dir]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("code/09_statistical_analysis/lib/common.R")   # assay_platform, assay_fluid, CHR_LEN_HG19
source("code/09_statistical_analysis/lib/hotspot_panels.R")
source("code/09_statistical_analysis/lib/fig_theme.R")      # theme_pgd, save_fig, pt_size
DEFAULT_OUT <- "output"

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[1] else DEFAULT_OUT
finemap_png <- file.path(OUT_DIR, "figures", "01_pqtl_architecture.png")

read_panel <- function(name) {
  f <- file.path(OUT_DIR, name)
  if (!file.exists(f)) {
    stop("Missing ", f, ": run 01a_pqtl_data.R / 01b_hotspot_data.R / 01b_hotspot_enrichment.R first.")
  }
  fread(f)
}

map_d <- read_panel("01a_fig1_map.csv")
prop_d <- read_panel("01a_fig1_chr_prop.csv")
hs_loci <- read_panel("01b_hotspot_loci.csv")
# Plasma Olink leads and their SomaScan match (13_somascan_olink_plasma.R).
plat_d <- read_panel("13_somascan_olink_signals.csv.gz")
stopifnot(nrow(plat_d) > 0, all(plat_d$comparison %in% names(PLATFORM_COL)))
plat_d[, comparison := factor(comparison, levels = names(PLATFORM_COL))]
plat_d[, signal_type := factor(signal_type, levels = c("cis", "trans"))]
stopifnot(hs_loci[is_hotspot_csf == TRUE, all(n_csf >= HOTSPOT_MIN)],
          hs_loci[is_hotspot_plasma == TRUE, all(n_plasma >= HOTSPOT_MIN)])
for (d in list(map_d, prop_d)) {
  d[, fluid := factor(fluid, levels = c("CSF", "Plasma"))]
}
map_d[, cis_trans := factor(cis_trans, levels = c("cis", "trans"))]

# --- Genome-wide plotting scales ------------------------------------------------
chr_len <- CHR_LEN_HG19
autosomes <- as.character(1:22)
chr_offset <- c(0, cumsum(chr_len[-length(chr_len)]))
names(chr_offset) <- names(chr_len)

autosome_end <- sum(chr_len[autosomes])
genome_end <- autosome_end
x_breaks <- unname(chr_offset[autosomes] + chr_len[autosomes] / 2)
x_labels <- autosomes
y_breaks <- x_breaks
y_labels <- autosomes

stripe_df <- data.frame(
  xmin = chr_offset[seq(1, length(autosomes), 2)],
  xmax = chr_offset[seq(1, length(autosomes), 2)] +
    chr_len[seq(1, length(autosomes), 2)]
)
diag_x <- seq(0, autosome_end, length.out = 500)
diag_df <- data.frame(
  x = diag_x,
  ymin = pmax(diag_x - 1, 0),
  ymax = pmin(diag_x + 1, autosome_end)
)

make_map <- function(d = map_d) {

  plot_d <- copy(d)
  plot_d[, color_group := fluid]
  color_values <- FLUID_COL

  p <- ggplot(plot_d, aes(
    variant_linear, gene_linear, color = color_group, shape = cis_trans
  )) +
    geom_rect(
      data = stripe_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "grey94", inherit.aes = FALSE
    ) +
    geom_ribbon(
      data = diag_df, aes(x = x, ymin = ymin, ymax = ymax),
      fill = "grey60", alpha = 0.25, inherit.aes = FALSE
    ) +
    geom_point(size = 1.2, alpha = 0.55) +
    # No coord_fixed here: a fixed aspect ratio locks the panel's width, so patchwork
    # cannot stretch it to the layout cell and centres it instead, leaving panel a
    # indented relative to panels b and c. Both axes are genomic position over the
    # same total span, so the cis diagonal stays a straight corner-to-corner line.
    coord_cartesian() +
    scale_color_manual(
      values = color_values,
      name = "Fluid type",
      drop = FALSE,
      guide = guide_legend(
        nrow = 1, title.position = "left",
        override.aes = list(size = 3, alpha = 1, shape = 16)
      )
    ) +
    scale_shape_manual(
      values = c(cis = 17, trans = 16),
      labels = c(cis = expression(italic("cis")), trans = expression(italic("trans"))),
      name = "pQTL type",
      guide = guide_legend(
        nrow = 1, title.position = "left",
        override.aes = list(size = 3, alpha = 1, color = "grey25")
      )
    ) +
    scale_x_continuous(
      breaks = x_breaks, labels = x_labels, limits = c(0, autosome_end),
      expand = expansion(mult = 0, add = 0),
      name = "Genomic position of variant"
    ) +
    scale_y_continuous(
      breaks = y_breaks, labels = y_labels, limits = c(0, genome_end),
      expand = expansion(mult = 0, add = 0),
      name = "Genomic position of protein-coding gene"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.spacing = unit(0.08, "lines"),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.x = element_text(size = 10, margin = margin(t = 9, b = 5)),
      axis.title.y = element_text(size = 10),
      strip.background = element_rect(fill = "white", colour = "grey70"),

      # Above the panel, as in panel b, so neither legend sits on top of data and the
      # two panels are keyed the same way. Inside the map it had to be nudged with a
      # large negative legend.margin, which broke whenever the panel was resized.
      legend.position = "top",
      # Nudged right of centre so the key sits over the map rather than over the
      # y-axis title; 0.5 would centre it on the whole plot, axis furniture included.
      legend.justification = "right",
      legend.margin = margin(t = 0, r = 0, b = -2, l = 0),
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9),

      plot.margin = margin(t = 4, r = 0, b = 4, l = 2)
    )

  p
}

make_fluid_prop <- function(d = prop_d) {
  prop_d <- copy(d)
  prop_d[, variant_chr := factor(as.character(variant_chr), levels = autosomes)]
  # Both the count and the share: the percentage alone hides how few signals some
  # chromosomes carry.
  prop_d[, label := sprintf("%.0f%%(%d)", 100 * prop, n)]
  # Labels are anchored to the outer end of their own segment (plasma at the left
  # edge, CSF at the right) rather than centred: a centred label on a narrow segment
  # spills across the colour boundary, or off the panel entirely.
  prop_d[, `:=`(
    label_x = fifelse(fluid == "CSF", 0.985, 0.015),
    hjust = fifelse(fluid == "CSF", 1, 0),
    # A segment narrower than its own label still reaches into the other colour, so
    # the few narrow ones (chr 9 CSF) get a smaller label rather than a misplaced one.
    label_size = fifelse(prop < 0.28, pt_size(5), pt_size(5.5))
  )]

  prop_bg <- data.frame(
    ymin = seq(0.5, length(autosomes) - 0.5, by = 2),
    ymax = seq(1.5, length(autosomes) + 0.5, by = 2)
  )

  ggplot(prop_d, aes(x = prop, y = variant_chr, fill = fluid)) +
    geom_rect(
      data = prop_bg,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "grey94", inherit.aes = FALSE
    ) +
    geom_col(width = 0.7, color = "white", linewidth = 0.2) +
    geom_text(
      aes(x = label_x, label = label, hjust = hjust, size = label_size),
      color = "white"
    ) +
    scale_size_identity() +
    scale_fill_manual(
      values = FLUID_COL,
      drop = FALSE,
      guide = "none"
    ) +
    scale_x_continuous(
      labels = function(x) round(100 * x),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      expand = expansion(mult = 0, add = 0)
    ) +
    # Clip rather than filter: the proportions are read back from CSV, so a
    # chromosome's two segments can stack to 1 + 1e-16 and scale limits would drop
    # the top segment entirely.
    coord_cartesian(xlim = c(0, 1), expand = FALSE) +
    # Axis on the left, between the map and the bars. On the right it was pushed off
    # its own panel by ~0.25 in: panel b's right-hand facet strips claim a column in
    # the shared patchwork gtable, and because patchwork aligns columns across rows
    # this row had to reserve the same width, stranding the axis beyond it. With the
    # axis on the left that column is simply empty here, and the bars get the space.
    # The y title is dropped as well -- it repeated the map's x-axis title verbatim.
    scale_y_discrete(drop = FALSE, limits = autosomes, position = "left", expand = c(0, 0.5)) +
    labs(x = "Proportion, % (n)", y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 5)),
      # Matched to the map panel so the two x-axis titles of panel a sit on one
      # baseline: the map has a taller gap (axis.text margin + title margin t = 9),
      # and the strip's rotated "100" label is taller than the map's "22", so the
      # strip needs a larger title margin, not the same one.
      axis.title.x = element_text(size = 10, margin = margin(t = 8, b = 5)),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 4, r = 4, b = 4, l = 2)
    )
}




# --- Assembly -------------------------------------------------------------------
# Per-panel margins rather than one `&` theme: `&` overrides plot.margin everywhere.
gutter <- function(l = 12, t = 5, b = 5, r = 5) {
  theme(plot.margin = margin(t = t, r = r, b = b, l = l))
}

# The hotspot panels default to a 12.6 in wide figure; these are 7.3 in, so they are
# asked for smaller text and gene-only labels.
# Axis titles at 7, as on the genome map, so the panels of one figure match.
hotspot_shrink <- theme(
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 7),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 7),
  legend.key.size = unit(0.7, "lines")
)

# 01_pqtl_architecture.png c: trans-pQTL hotspot loci.
make_trans_hotspot_fig <- function() {
  make_hotspot_map(hs_loci,
    label_size = pt_size(6), headroom = 0.45, short_labels = TRUE,
    # c is 9/13 of the page width (b at its left), so close pairs are pushed apart
    # further: chr3 DCBLD2/KNG1 and chr12/14 GNPTAB/SERPINA1 in plasma.
    near_frac = 0.1, pair_push = 0.5,
    base_size = 11, title_size = 8, legend_size = 7.5,
    # t = 14: at t = 0 the key sat right under a's x-axis title.
    y_title = "Trans-target proteins") + hotspot_shrink + gutter(8, t = 14) +
    # Key above the plot (make_hotspot_map's default), with no white boxes behind it.
    theme(legend.background = element_blank(), legend.key = element_blank())
}

# 01_pqtl_architecture.png b, top: per Olink plasma lead, the closest same-type SomaScan
# signal for the same target (13's comparison), as one stacked bar per cis / trans.
make_platform_overlap <- function(d = plat_d) {
  n <- d[, .N, by = .(signal_type, comparison)]
  tot <- d[, .N, by = signal_type]
  ggplot(n, aes(N, signal_type, fill = comparison)) +
    geom_col(position = position_fill(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = PLATFORM_COL, name = NULL,
                      labels = c(`exact lead` = "Same lead", `same locus` = "Within 1 Mb",
                                 `different locus` = "Beyond 1 Mb",
                                 `no SomaScan signal` = "No signal"),
                      guide = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_continuous(labels = function(x) paste0(100 * x, "%"), breaks = c(0, 0.5, 1),
                       expand = expansion(mult = c(0, 0.06))) +
    scale_y_discrete(limits = rev(levels(d$signal_type)),
                     labels = function(x) parse(text = paste0(
                       "italic(", x, ")~'('*italic(n)*' = ", format(tot[match(x, signal_type), N], big.mark = ","), ")'"))) +
    labs(x = "Proportion of Olink leads", y = NULL) +
    theme_bw(base_size = 7) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.ticks.y = element_blank(), axis.text = element_text(size = 6),
          axis.title.x = element_text(size = 7),
          legend.position = "bottom", legend.justification = "left",
          legend.margin = margin(-2, 0, 0, 0), legend.key.size = unit(6, "pt"))
}

# 01_pqtl_architecture.png b, bottom: Olink and SomaScan effects of the exact leads with
# harmonized alleles; plasma orange, cis / trans shapes as in a.
make_platform_scatter <- function(d = plat_d) {
  s <- d[!is.na(concordant_direction)]
  lab <- s[, .(k = sum(concordant_direction), n = .N), by = signal_type][order(signal_type)]
  lab <- paste0("atop('Same direction', ", paste0("italic(", lab$signal_type, ")~'", lab$k, "/", lab$n, "'",
                                                  collapse = "*','~"), ")")
  lim <- range(c(s$olink_beta, s$soma_beta_aligned))
  ggplot(s, aes(olink_beta, soma_beta_aligned, shape = signal_type)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(colour = FLUID_COL[["Plasma"]], size = 1, alpha = 0.7) +
    annotate("text", x = lim[1], y = lim[2], label = lab, parse = TRUE,
             hjust = 0, vjust = 1, size = pt_size(6)) +
    scale_shape_manual(values = c(cis = 17, trans = 16), guide = "none") +
    coord_cartesian(xlim = lim, ylim = lim) +
    labs(x = "Olink effect", y = "SomaScan effect") +
    theme_bw(base_size = 7) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 6),
          axis.title = element_text(size = 7))
}

# 01_pqtl_architecture.png: the genome map (a) above the plasma Olink--SomaScan
# comparison (b) and the trans hotspots (c) right of b, in the order the text cites
# them; heights fit the 170-mm page.
make_section1_fig <- function() {
  small <- theme(axis.title.x = element_text(size = 7, margin = margin(t = 6, b = 3)),
                 axis.title.y = element_text(size = 7),
                 legend.title = element_text(size = 7),
                 legend.text = element_text(size = 6))
  # One grid, not two wrapped pictures: wrapped parts are never aligned, so b's frame
  # sat left of a's, and a's legend (collected outside its sub-patchwork) fell below
  # the tag while b's sat above it. Here each part keeps its legend in its own plot
  # (the proportion strip shares the map's fluid colours and drops its own), and every
  # tag sits at the top left of its plot, above the legend.
  map  <- make_map() + gutter(b = -4) + small + labs(tag = "a") +
    theme(legend.position = "top")
  prop <- make_fluid_prop() + gutter(b = -4) + small + theme(legend.position = "none")
  hot  <- make_trans_hotspot_fig() + labs(tag = "c")
  ovl  <- make_platform_overlap() + gutter(t = 14, b = 2) + labs(tag = "b")
  # One analysis, one tag: the scatter is the lower half of b.
  sct  <- make_platform_scatter() + gutter(t = 6)
  # One 13-column grid so b's left edge stays aligned with a's: a:strip 10:3 as before,
  # b:c 4:9.
  (map + prop + hot + ovl + sct) +
    plot_layout(design = "AAAAAAAAAABBB\nDDDDCCCCCCCCC\nEEEECCCCCCCCC",
                heights = c(100, 20, 50)) &
    theme_pgd() &
    theme(plot.tag = element_text(size = 8, face = "bold"),
          plot.tag.position = "topleft", plot.tag.location = "plot")
}

# --- Output ----------------------------------------------------------------------
# Drawn at the final double-column width (lib/fig_theme.R); PNG
# for the HTML manuscript.
save_fig(make_section1_fig(), file.path(OUT_DIR, "figures"), "01_pqtl_architecture", 170)
