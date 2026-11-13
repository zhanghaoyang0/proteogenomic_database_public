#!/usr/bin/env Rscript
# Locus panels of the redesigned protein-disease figure, from the tables 05f writes:
# one column per case (CD58 - multiple sclerosis, GRN - Alzheimer disease), stacking
# the disease GWAS, CSF pQTL and plasma pQTL, each as -log10 P coloured by LD r2 with
# the disease lead above a strip with the SuSiE PIP of that track's credible-set
# variants, and the protein's own gene below. Vertical lines mark the variants in both the
# disease and the CSF credible sets, the ones the colocalization rests on (CD58 five,
# all in one LD block; GRN rs5848 alone); a diamond marks such a variant in each track
# whose own credible set contains it. The view is the disease lead +/- 250 kb (the SuSiE
# window is +/- 500 kb, too wide to label the genes). CD58 is shown at
# 117.00-117.16 Mb and GRN at the disease lead +/-100 kb. Look after the GWAS / eQTL locus
# panels the author pointed to (2026-10-27; Fig. 3c,f of doi:10.1038/s41588-024-01896-3).
#
# Plotting only. Panels b-c of 05_protein_disease.png: 05c_plot.R sources this file
# with PANEL_EXPORT_ONLY set and calls case_column(). Run on its own it saves a preview,
# tmp/proteogenomic_database/05g_locus.png at the repo root.
# Usage: Rscript code/09_statistical_analysis/05g_locus_plot.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("code/09_statistical_analysis/lib/fig_theme.R")
PGD <- "."
OUT_DIR <- "tmp/proteogenomic_database"
dir.create(OUT_DIR, showWarnings = FALSE)

tr    <- fread(file.path(PGD, "output", "05f_locus_tracks.csv.gz"))
genes <- fread(file.path(PGD, "output", "05f_locus_genes.csv"))
meta  <- fread(file.path(PGD, "output", "05f_locus_meta.csv"))

CASE_ORDER <- c("CD58", "GRN")
DISEASE_SHORT <- c("Multiple sclerosis" = "MS", "Alzheimer disease" = "AD")
# Colours, LD bins and line styles from lib/fig_theme.R, shared with Fig. 5c-e (10d).
NA_COL <- "grey80"
GENE_COL <- "grey15"
CD58_VIEW <- c(117e6, 117.16e6)
GRN_VIEW_HALF <- 1e5

tr[, `:=`(mlp = -log10(pmax(p, 1e-300)), ld = cut(r2, LD_BREAKS, LD_LABELS, right = FALSE))]

ax_theme <- theme_classic(base_size = 7) +
  theme(axis.text = element_text(size = 5.6, colour = "black"),
        axis.title = element_text(size = 6),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(size = 6.5, hjust = 0.5, margin = margin(b = 1)),
        plot.title.position = "panel",
        plot.margin = margin(1, 2, 1, 2))
no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.title.x = element_blank(), axis.line.x = element_blank())

case_column <- function(g) {
  d  <- tr[gene == g]
  m  <- meta[gene == g]
  dz <- DISEASE_SHORT[[m$disease[1]]]
  lead_snp <- m$disease_lead[1]
  lead_pos <- d[track == "Disease" & snp == lead_snp, pos[1]]
  if (is.na(lead_pos)) stop(g, ": disease lead ", lead_snp, " not in the disease track")
  view_bp <- switch(g,
                    CD58 = CD58_VIEW,
                    GRN = lead_pos + c(-GRN_VIEW_HALF, GRN_VIEW_HALF),
                    stop("no display window for ", g))
  xr <- view_bp / 1e6
  x_breaks <- if (g == "CD58") c(117.00, 117.04, 117.08, 117.12) else scales::pretty_breaks(4)(xr)
  d <- d[pos / 1e6 >= xr[1] & pos / 1e6 <= xr[2]]
  shared <- intersect(d[track == "Disease" & in_cs, pos], d[track == "CSF" & in_cs, pos])
  stopifnot(length(shared) > 0)
  shared_all <- intersect(tr[gene == g & track == "Disease" & in_cs, pos],
                          tr[gene == g & track == "CSF" & in_cs, pos])
  if (!all(shared_all >= view_bp[1] & shared_all <= view_bp[2]))
    stop(g, ": display window excludes disease/CSF credible-set variants")
  vl <- shared / 1e6
  xs <- scale_x_continuous(limits = xr, breaks = x_breaks,
                           labels = function(x) sprintf("%.2f", x), expand = c(0, 0))
  # Variants in both the disease and CSF credible sets: one line for a single variant
  # (GRN); for several, one light band from the first to the last (CD58: five variants
  # in one LD block), since five near-touching lines read as five signals. The shared
  # variants themselves stay marked by the diamonds in the association tracks.
  vlines <- if (length(vl) == 1L)
    geom_vline(xintercept = vl, colour = LEAD_LINE_COL, linewidth = 0.25, alpha = 0.7) else
    annotate("rect", xmin = min(vl), xmax = max(vl), ymin = -Inf, ymax = Inf,
             fill = LEAD_LINE_COL, alpha = 0.15)
  # the shaded band is explained in the caption; a label in it collided with the rsID

  assoc <- function(t, title) {
    dt <- d[track == t]
    setorder(dt, r2, na.last = FALSE)                     # high-LD points drawn on top
    # The top 25% of the track is left for the lead's rsID, above both the data and the
    # genome-wide line.
    ytop <- max(dt$mlp, -log10(5e-8)) / 0.75
    p <- ggplot(dt, aes(pos / 1e6, mlp)) + GW_LINE + vlines +
      geom_blank(data = data.frame(x = xr[1], y = ytop), aes(x, y)) +
      geom_point(aes(colour = ld), size = 0.8, stroke = 0) +
      geom_point(data = dt[pos %in% shared & in_cs], shape = 23, size = 1.5,
                 fill = LEAD_COL, colour = "black", stroke = 0.2) +
      scale_colour_manual(values = setNames(LD_COLS, LD_LABELS), breaks = LD_LABELS,
                          na.value = NA_COL, drop = FALSE,
                          name = expression(italic(r)^2 ~ "with disease lead"),
                          guide = guide_legend(nrow = 1, override.aes = list(size = 1.6))) +
      # ticks over the data only, about three to a track (the tracks are ~8 mm tall);
      # bottom headroom keeps the 0 label clear of the PIP strip's 1 below it
      scale_y_continuous(breaks = function(l) pretty(c(0, l[2] * 0.6), 2),
                         expand = expansion(mult = c(0.1, 0))) +
      # Same wording as the Fig. 3b and 4c,d locus plots; expression(-log[10]...) set the
      # minus apart as "− log".
      # Track title centred above the panel (not over the y axis), as in Fig. 3b and 4c,d.
      xs + labs(y = "−log₁₀ P", title = parse(text = title)[[1]]) + ax_theme + no_x
    if (t == "Disease") p <- p +
      # centred above the lead, clear of the neighbouring credible-set diamonds; a white
      # box so the credible-set lines pass behind it rather than through it
      geom_label(data = dt[pos == lead_pos], aes(label = snp), hjust = 0.5, vjust = -0.6,
                 size = pt_size(5), family = FIG_FAMILY, fill = "white", linewidth = 0,
                 label.padding = unit(0.5, "pt"))
    p
  }
  pip <- function(t) {
    dt <- d[track == t & in_cs]
    p <- ggplot(dt, aes(pos / 1e6, pip)) +
      geom_hline(yintercept = 1, colour = "grey85", linewidth = 0.2) + vlines +
      # padded so its 1 and 0 labels clear the 0 of the track above and the next title
      xs + scale_y_continuous(limits = c(-0.3, 1.3), breaks = c(0, 1), expand = c(0, 0)) +
      # Vertical like the -log10 P titles (2026-11-12): horizontal, it widened the title
      # column and pushed every -log10 P title away from its axis.
      labs(y = "PIP") + ax_theme + no_x
    if (nrow(dt)) p + geom_point(shape = 23, size = 1, fill = LEAD_COL, colour = "black",
                                 stroke = 0.15)
    else p + annotate("text", x = xr[1], y = 0.5, label = "No fine-mapped signal", hjust = -0.05,
                      size = pt_size(5), colour = "grey40", family = FIG_FAMILY)
  }

  # the protein's own gene only, named below its bar; the other genes in the window
  # (20 at GRN) are too small to name and are left out
  gg <- genes[gene == g & symbol == g]
  stopifnot(nrow(gg) == 1)
  gg[, `:=`(x0 = pmax(start, xr[1] * 1e6), x1 = pmin(end, xr[2] * 1e6))]
  # label kept inside the view, so a gene running off either edge keeps its name
  gg[, lx := pmin(pmax((x0 + x1) / 2e6, xr[1] + 0.04 * diff(xr)), xr[2] - 0.04 * diff(xr))]
  gp <- ggplot(gg) +
    geom_segment(aes(x = x0 / 1e6, xend = x1 / 1e6, y = 0, yend = 0),
                 colour = GENE_COL, linewidth = 1.2) +
    geom_text(aes(x = lx, y = -0.35, label = symbol), fontface = "bold.italic",
              colour = GENE_COL, size = pt_size(5), family = FIG_FAMILY, vjust = 1) +
    scale_y_continuous(limits = c(-1, 0.4)) +
    coord_cartesian(clip = "off") +
    xs + labs(x = sprintf("Chromosome %s (Mb)", d$chr[1]), y = NULL) + ax_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

  # a gap after each track's PIP strip, so the association and PIP rows read as one
  # group per track (disease, CSF, plasma)
  gap <- plot_spacer()
  wrap_plots(
    assoc("Disease", sprintf("%s~GWAS", dz)), pip("Disease"), gap,
    assoc("CSF", sprintf("CSF~'%s'~pQTL", g)), pip("CSF"), gap,
    assoc("Plasma", sprintf("Plasma~'%s'~pQTL", g)), pip("Plasma"),
    gp, ncol = 1, heights = c(3.2, 1.3, 0.35, 3.2, 1.3, 0.35, 3.2, 1.3, 1))
}

stopifnot(setequal(unique(tr$gene), CASE_ORDER))
if (!exists("PANEL_EXPORT_ONLY")) {
# the identical r2 scales of the six association panels collect into one legend
fig <- (case_column("CD58") | case_column("GRN")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))
fig <- fig & theme(text = element_text(family = FIG_FAMILY))
save_fig(fig, OUT_DIR, "05g_locus", height_mm = 95)
message("05g_locus_plot done")
}
