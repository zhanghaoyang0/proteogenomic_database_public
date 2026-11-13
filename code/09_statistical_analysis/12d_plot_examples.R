#!/usr/bin/env Rscript
# CSF against plasma for one protein, after Lu et al.'s rs5848 Fig. 1 and 2: the section 12
# examples (docs/pipeline.md). Plotting only, from the tables 11a and 11c
# write in protein mode (PROTEINS=..., OUT_NAME 11a_protein_brain / 11c_protein_clinical).
# One block of three panels per gene (comma-separated), CSF left and plasma right:
#   baseline subcortical volumes; regional cortical thickness x time (std beta per SD of
#   protein level); mPACC by protein tertile (CU) and event-free probability for AD-CI.
# The supplementary figure is GRN alone (panels a-c; 2026-11-12: NPTX2 and NCAM1 dropped,
# their clinical panels repeat Fig. 6):
#   Rscript code/09_statistical_analysis/12d_plot_examples.R GRN output 12_supp_examples.png
# (the clinical examples of Fig. 7, then the AD genetically supported protein; LYPD3 left
# with FinnGen's wide AD definition, 2026-11-01; one figure since 2026-11-11). It reads only summary tables, so it runs
# on the cluster too, with the house fonts the Mac lacks; there sf needs
# PROJ_DATA=$HOME/miniforge3/envs/gwas/share/proj or the brain maps fail.
# Usage (repo root): Rscript code/09_statistical_analysis/12d_plot_examples.R <gene[,gene]> [out_dir] [fig_name]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(sf)
  library(ggseg)
})
source("code/09_statistical_analysis/lib/brain_maps.R")
source("code/09_statistical_analysis/lib/fig_theme.R")   # theme_pgd, save_fig, p_expr, pt_size
IN_DIR <- "output"
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript code/09_statistical_analysis/12d_plot_examples.R <gene[,gene]> [out_dir] [fig_name]")
GENES <- strsplit(args[1], ",")[[1]]
OUT_DIR <- if (length(args) >= 2) args[2] else IN_DIR
FIG_NAME <- sub("\\.png$", "", if (length(args) >= 3) args[3] else
  paste0("12_supp_examples_", paste(GENES, collapse = "_")))   # save_fig adds .png
rd <- function(f) { p <- file.path(IN_DIR, f); if (!file.exists(p)) stop("missing ", p); fread(p) }
BA <- rd("11a_protein_brain_assoc.csv"); CT <- rd("11c_protein_clinical_assoc.csv")
TJ <- rd("11c_protein_clinical_traj.csv"); SV <- rd("11c_protein_clinical_surv.csv")

panels <- function(GENE) {

  ba <- BA[startsWith(snp, paste0(GENE, "__"))]
  ct <- CT[startsWith(snp, paste0(GENE, "__"))]
  tj <- TJ[startsWith(snp, paste0(GENE, "__"))]
  sv <- SV[startsWith(snp, paste0(GENE, "__")) & outcome == "AD-CI"]
  fluid_of <- function(a) fifelse(grepl("csf", a, ignore.case = TRUE), "CSF", "Plasma")
  for (d in list(ba, ct, tj, sv)) d[, fluid := factor(fluid_of(snp), levels = c("CSF", "Plasma"))]
  if (uniqueN(ba$fluid) != 2) stop(GENE, ": need one CSF and one plasma assay in 11a_protein_brain_assoc.csv")
  if (ba[analysis == "longitudinal" & measure == "cortical thickness", all(fit == "skipped")])
    stop("regional longitudinal models were skipped: rerun 11a with LONG_ALL=1")

  map_row <- function(meas, ana, atlas, to_label, title, fill_lab, label_ns, bar_pt, max_labels) {
    d <- ba[measure == meas & analysis == ana]
    d[, label := to_label(region)]
    lim <- MAP_LIM[[meas]]
    # One colour scale per outcome, the same for every gene (MAP_LIM), so colours compare
    # across the gene blocks; collected to the right of the two maps: two bottom
    # colourbars side by side ran into each other.
    wrap_plots(lapply(c("CSF", "Plasma"), function(fl)
      brain(d[fluid == fl], atlas, lim, paste(GENE, "|", fl, "|", title), fill_lab, label_ns = label_ns,
            max_labels = max_labels) |>
        style_map()), nrow = 1, guides = "collect") & theme_map_row &
      theme(legend.key.height = unit(bar_pt / 5, "pt"))   # the bar is five keys long
  }
  ra <- map_row("subcortical volume", "cross-sectional", aseg, to_aseg, "Baseline subcortical volume",
                "Std. beta\nper SD", FALSE, 36, 8)
  rc <- map_row("cortical thickness", "longitudinal", dk, to_dk, "Longitudinal cortical thickness",
                "Std. beta per\nSD per year", FALSE, 36, 4)
  # (thickness: at most the 4 smallest-P names per map; more piled up in the short row)

  TER <- c(`0` = "Low", `1` = "Middle", `2` = "High")
  TCOL <- LEVEL_COL   # as in Fig. 7 (12e)
  tj[, tertile_lab := factor(TER[as.character(genotype)], levels = TER)]
  sv[, tertile_lab := factor(TER[as.character(genotype)], levels = TER)]
  # Label wording as in Fig. 7 (12e).
  lab_m <- ct[model == "mPACC, CU" & term == "years:g",
              .(fluid, lab = sub("italic(P)", "italic(P)[interaction]", p_expr(p), fixed = TRUE))]
  lab_s <- ct[model == "Cox AD-CI, CU" & term == "g",
              .(fluid, lab = sprintf('HR~"%.2f;"~%s*";"~"%d events"', est, p_expr(p), events))]
  stopifnot(nrow(lab_m) == 2, nrow(lab_s) == 2)
  pm <- ggplot(tj, aes(years, fit, colour = tertile_lab, fill = tertile_lab)) +
    geom_ribbon(aes(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.45) +
    geom_label(data = lab_m, aes(x = -Inf, y = -Inf, label = lab), inherit.aes = FALSE, parse = TRUE,
              hjust = -0.06, vjust = -0.6, size = TXT, colour = "grey20", fill = "white",
              border.colour = "black", linewidth = .25, label.padding = unit(1.5, "pt")) +
    facet_wrap(~ fluid) +
    scale_colour_manual(values = TCOL, name = "Protein level") +
    scale_fill_manual(values = TCOL, guide = "none") +
    labs(x = "Years from baseline", y = "mPACC", title = paste(GENE, "| Cognition")) + theme_line
  # Bottom left, inset from the frame: survival curves are still near 1 there.
  ps <- ggplot(sv, aes(years, surv, colour = tertile_lab, fill = tertile_lab)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.13, colour = NA) +
    geom_step(linewidth = 0.45) +
    geom_label(data = lab_s, aes(x = -Inf, y = -Inf, label = lab), inherit.aes = FALSE, parse = TRUE,
              hjust = -0.06, vjust = -0.6, size = TXT, colour = "grey20", fill = "white",
              border.colour = "black", linewidth = .25, label.padding = unit(1.5, "pt")) +
    facet_wrap(~ fluid) +
    scale_colour_manual(values = TCOL, name = "Protein level") +
    scale_fill_manual(values = TCOL, guide = "none") +
    labs(x = "Years from baseline", y = "Free of A\u03b2+\nMCI or dementia",
         title = paste(GENE, "| Progression")) + theme_line
  list(ra = ra, rc = rc, pm = font_layers(pm), ps = font_layers(ps))
}

# Drawn at the final size (183 mm wide): text 5-7 pt, panel letters 8 pt bold.
TXT <- pt_size(6)
theme_line <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        axis.text = element_text(size = 6), legend.position = "bottom",
        plot.title = element_text(size = 7, hjust = 0.5))
theme_map_row <- theme_pgd() + theme(legend.position = "right", legend.direction = "vertical",
  legend.title = element_text(size = 6, vjust = .5), legend.text = element_text(size = 5.5),
  # Title beside the bar: a collected legend gets only the panel's height, and title
  # above bar overflowed the short thickness row (its bar was cut off).
  legend.title.position = "left",
  legend.key.width = unit(5, "pt"), legend.margin = margin(0, 0, 0, 0),
  # Side margins keep the two long titles of a row apart (the fixed-aspect maps are
  # narrower than their titles).
  plot.margin = margin(2, 12, 2, 12),
  plot.title = element_text(size = 7, hjust = 0.5))
# font_layers() is in lib/fig_theme.R.
# Significant region names only (label_ns = FALSE above, as in Fig. 7): at this size the
# grey names of all ~30 other subcortical regions piled into an unreadable layer. The
# boxed names get more solver time so they do not touch in the small maps.
style_map <- function(p) {
  for (l in p$layers) if (inherits(l$geom, "GeomLabelRepel")) {
    l$aes_params$size <- pt_size(5.5); l$geom_params$box.padding <- unit(.15, "lines")
    l$geom_params$force <- 3; l$geom_params$max.time <- 3; l$geom_params$max.iter <- 1e5 }
  font_layers(p)
}

# Symmetric colour limits per imaging outcome over all genes in the figure.
MAP_LIM <- sapply(c("subcortical volume", "cortical thickness"), function(m)
  BA[sub("__.*", "", snp) %in% GENES & measure == m &
       analysis == c(`subcortical volume` = "cross-sectional", `cortical thickness` = "longitudinal")[[m]],
     max(abs(std_beta), na.rm = TRUE)], simplify = FALSE)
stopifnot(all(is.finite(unlist(MAP_LIM))))
P <- setNames(lapply(GENES, panels), GENES)
# Map rows are wrapped so their tags sit at the row's left edge, like the line rows'
# (on a fixed-aspect map the tag would follow the panel inward). The wrapped rows
# already carry theme_pgd(), which the final `&` does not reach.
tag_map <- function(row, t) wrap_elements(full = row) + labs(tag = t)
tag_first <- function(row, t) { row[[1]] <- row[[1]] + labs(tag = t); row }
# Per gene: subcortical maps, thickness maps, then cognition and progression. One
# "Protein level" legend, under the last gene: the tertile colours are the same throughout.
rows <- unlist(lapply(seq_along(GENES), function(k) {
  x <- P[[GENES[k]]]
  leg <- if (k == length(GENES)) "bottom" else "none"
  lines <- (x$pm + x$ps + plot_layout(guides = "collect")) & theme(legend.position = leg)
  list(tag_map(x$ra, letters[3 * k - 2]), tag_map(x$rc, letters[3 * k - 1]),
       tag_first(lines, letters[3 * k]))
}), recursive = FALSE)
heights <- rep(c(1.2, 0.72, 1), length(GENES))   # 0.64 clipped the 5-tick thickness bar
heights[length(heights)] <- heights[length(heights)] + 0.12   # room for the legend
fig <- wrap_plots(rows, ncol = 1, heights = heights) & theme_pgd()
# About 80 mm per gene; three genes fill a supplementary page.
save_fig(fig, OUT_DIR, FIG_NAME, 8 + 81 * length(GENES), max_h_mm = FIG_H_MAX_SUPP_MM)
message("12d done: ", file.path(OUT_DIR, paste0(FIG_NAME, ".png")))
