#!/usr/bin/env Rscript
# Supplementary tissue-context figure from 04's exported counts (plotting only), 183 mm
# wide (2026-11-11: the former brain and peripheral figures in one):
#   04_supp_coloc.png  MetaBrain regions, GTEx brain regions, SingleBrain cell types and
#                      GTEx peripheral tissues, one panel each with its own x axis
#     a  colocalizing proteins, inclusive (counts)
#     b  colocalizing proteins, head-to-head (% of the cis-both denominator)
# The numbers at the right of each row are the counts (CSF, plasma) in both panels. The
# main-text tissue panels are in 10d.
# Usage (repo root): Rscript code/09_statistical_analysis/04c_plot_tissue_context.R
# (04_protein_tissue.R sources this file and calls plot_tissue_context()).
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("code/09_statistical_analysis/lib/fig_theme.R")          # theme_pgd, save_fig, pt_size
source("code/09_statistical_analysis/lib/eqtl_coloc_scope.R")   # CELLTYPE_LAB

plot_tissue_context <- function(output_dir) {
  counts <- fread(file.path(output_dir, "04_tissue_dumbbell.csv"))
  sets <- fread(file.path(output_dir, "03_analysis_sets.csv"))
  denom <- sets[step == "cis_both", n]   # cis pQTL in both fluids (03)
  cols <- FLUID_COL
  block_labels <- c(MetaBrain = "MetaBrain", `GTEx Brain` = "GTEx brain",
                    SingleBrain = "SingleBrain", `GTEx Peripheral` = "GTEx peripheral")
  counts[, block := factor(block_labels[fifelse(resource == "GTEx", paste(resource, tissue_class),
                                                resource)], levels = block_labels)]
  stopifnot(!anyNA(counts$block))
  base <- theme_bw(base_size = 7) + theme(panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(), strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    plot.title = element_text(hjust = 0.5),
    legend.margin = margin(t = 0), legend.box.spacing = unit(1, "pt"),
    legend.position = "bottom")
  # One column per denominator. Each resource is its own panel with its own x axis
  # (2026-11-12): on one shared axis MetaBrain cortex (180) squeezed the GTEx brain and
  # SingleBrain rows into the first fifth. Separate panels also keep a region name shared
  # by two resources (cortex in MetaBrain and GTEx) as two rows.
  dots <- function(d, title, counts_axis) {
    d <- copy(d)
    d[, label := tools::toTitleCase(gsub("_", " ", sub("^brain_", "", tissue)))]
    d[resource == "SingleBrain", label := CELLTYPE_LAB[tissue]]
    d[label == "Basalganglia", label := "Basal ganglia"]
    d[label == "Spinalcord", label := "Spinal cord"]
    d[, x := if (counts_axis) as.numeric(n) else 100 * n / denom]
    blocks <- levels(droplevels(d$block))
    ps <- lapply(seq_along(blocks), function(k) {
      db <- d[block == blocks[k]]
      ord <- db[, .(total = sum(n)), by = label][order(total)]$label
      db[, row := factor(label, levels = ord)]
      seg <- dcast(db, block + row ~ fluid, value.var = "x")
      stopifnot(!anyNA(seg), nrow(db) == 2 * nrow(seg))
      # Counts as two columns right of the dots (CSF, then plasma): printed beside the
      # points they collided with the neighbouring rows once 36 tissues share 170 mm.
      x_max <- max(db$x)
      db[, x_lab := x_max * fifelse(fluid == "CSF", 1.16, 1.31)]
      last <- k == length(blocks)
      ggplot(db, aes(x, row)) +
        geom_segment(data = seg, aes(x = CSF, xend = Plasma, y = row, yend = row),
                     inherit.aes = FALSE, colour = "grey70", linewidth = 0.4) +
        geom_point(aes(colour = fluid), size = 1.3) +
        geom_text(aes(x_lab, label = n, colour = fluid), size = pt_size(5), hjust = 1,
                  family = FIG_FAMILY) +
        facet_wrap(~ block) +   # resource name above the panel, not beside it
        scale_colour_manual(values = cols, guide = "none") +
        scale_x_continuous(limits = c(0, x_max * 1.33), expand = expansion(mult = c(.02, 0)),
                           breaks = scales::breaks_pretty(4)(c(0, x_max))) +
        labs(title = if (k == 1) title, y = NULL,
             x = if (last) { if (counts_axis) "Colocalizing proteins (count)" else
               sprintf("Colocalizing proteins (%% of %d)", denom) }) +
        base + theme_pgd() +
        theme(strip.text = element_text(size = 6.5, hjust = 0, margin = margin(b = 2)),
              axis.text = element_text(size = 6), legend.position = "none")
    })
    # Panel heights follow the number of rows (plus room for the title and x title).
    n_rows <- vapply(blocks, function(b) uniqueN(d[block == b, label]), 1L)
    wrap_plots(ps, ncol = 1, heights = n_rows)
  }
  fig <- function(group, height_mm, name, max_h_mm = FIG_H_MAX_MM) {
    d <- counts[block %in% group]
    tag <- theme(plot.tag = element_text(size = 8, face = "bold", family = FIG_FAMILY))
    # Each column is wrapped so the panel letter covers the column, not every resource.
    p <- (wrap_elements(full = dots(d[analysis == "inclusive"], "Inclusive", TRUE)) + labs(tag = "a") + tag) |
         (wrap_elements(full = dots(d[analysis == "head_to_head"], "Head-to-head", FALSE)) + labs(tag = "b") + tag)
    save_fig(p, file.path(output_dir, "figures"), name, height_mm, max_h_mm = max_h_mm)
  }
  # 61 rows: the height of a full supplementary page.
  fig(block_labels, 235, "04_supp_coloc", max_h_mm = FIG_H_MAX_SUPP_MM)
}
if (sys.nframe() == 0L) plot_tissue_context("output")
