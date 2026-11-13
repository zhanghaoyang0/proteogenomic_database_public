#!/usr/bin/env Rscript
# Integrated clinical main figure from existing summary outputs only.
#   a     proteins by the A/T stage at which their CSF level changes (Pichet Binette et al.
#         2024; 12f_pathology_proteins.R): clinical effect heterogeneity, overall and by outcome type
#   b     NCAM1 brain maps;  c, d  NPTX2 and NCAM1 mPACC trajectories and progression
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ggrepel); library(sf); library(ggseg)
})
source("code/09_statistical_analysis/lib/brain_maps.R")
source("code/09_statistical_analysis/lib/fig_theme.R")   # theme_pgd, save_fig, p_expr
out_dir <- "output"
rd <- function(f) fread(file.path(out_dir, f))
ba <- rd("11a_protein_brain_assoc.csv")
ct <- rd("11c_protein_clinical_assoc.csv")
tj <- rd("11c_protein_clinical_traj.csv")
sv <- rd("11c_protein_clinical_surv.csv")[outcome == "AD-CI"]
rk <- rd("11c_protein_clinical_risk.csv")[outcome == "AD-CI"]
pp <- rd("12f_pathology_proteins.csv"); pt <- rd("12f_tests.csv")
genes <- c("NPTX2", "NCAM1")
for (d in list(ba, ct, tj, sv, rk)) {
  d[, gene := sub("__.*", "", snp)]
  d[, fluid := factor(fifelse(grepl("csf", snp, ignore.case = TRUE), "CSF", "Plasma"),
                      levels = c("CSF", "Plasma"))]
}
# Tertiles on one sequential hue: blue and orange mean CSF and plasma in every other
# figure, and the brain maps use a blue-red diverging scale.
cols <- LEVEL_COL
ter <- c(`0` = "Low", `1` = "Middle", `2` = "High")
for (d in list(tj, sv)) d[, level := factor(ter[as.character(genotype)], levels = names(cols))]
# Drawn at the final size (183 mm wide): text 5-7 pt, panel letters 8 pt bold.
TXT <- 6 / .pt   # geom_text size for 6-pt text
TXT7 <- 7 / .pt  # the A/T schematic beside a
# font_layers() is in lib/fig_theme.R.
# Region-name boxes: the default repel run lets the four subcortical names touch in
# a panel this narrow; push harder and run longer instead.
pad_labels <- function(p) {
  for (l in p$layers) if (inherits(l$geom, "GeomLabelRepel")) {
    l$geom_params$box.padding <- unit(.15, "lines"); l$geom_params$force <- 3
    l$geom_params$max.time <- 3; l$geom_params$max.iter <- 1e5 }
  p
}
base <- theme_bw(base_size = 7) + theme(panel.grid.minor = element_blank(),
  strip.background = element_blank(), strip.text = element_text(size = 7),
  axis.text = element_text(size = 6), legend.text = element_text(size = 6),
  legend.title = element_text(size = 7), legend.key.height = unit(6, "pt"),
  plot.title = element_text(size = 7, face = "plain", hjust = .5), legend.position = "bottom",
  plot.tag = element_text(size = 8, face = "bold"))

# d and e are rows (one outcome each) of two protein blocks, CSF | plasma inside each block;
# the second block shares the first one's y range and drops its y axis text.
mpacc_lim <- tj[gene %in% genes & years <= 5, range(c(fit - 1.96 * se, fit + 1.96 * se))]
curve <- function(g, progression = FALSE, show_y = TRUE) {
  if (progression) {
    d <- sv[gene == g]
    lab <- ct[gene == g & model == "Cox AD-CI, CU" & term == "g",
      .(fluid, lab = sprintf('HR~"%.2f;"~%s*";"~"%d events"', est, p_expr(p), events))]
    p <- ggplot(d, aes(years, surv, colour = level, fill = level)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .13, colour = NA) +
      geom_step(linewidth = .45) + coord_cartesian(xlim = c(0,6), ylim = c(.6,1)) +
      labs(y = "Free of Aβ+\nMCI or dementia")
  } else {
    d <- tj[gene == g]
    lab <- ct[gene == g & model == "mPACC, CU" & term == "years:g",
      .(fluid, lab = sub("italic(P)", "italic(P)[interaction]", p_expr(p), fixed = TRUE))]
    p <- ggplot(d, aes(years, fit, colour = level, fill = level)) +
      geom_ribbon(aes(ymin = fit - 1.96*se, ymax = fit + 1.96*se), alpha = .15, colour = NA) +
      geom_line(linewidth = .45) + coord_cartesian(xlim = c(0,5), ylim = mpacc_lim) + labs(y = "mPACC")
  }
  stopifnot(nrow(lab) == 2, uniqueN(d$fluid) == 2)
  # Bottom left: survival curves are still near 1 there, trajectories lie above it.
  # Inset from the frame so the text clears the bottom tick label and gridline.
  # One line per panel (HR, P and events together), so there is no second line to align.
  # White box with a thin black frame, so the curves and bands behind do not cut through it.
  p <- p + geom_label(data = lab, aes(x = -Inf, y = -Inf, label = lab), inherit.aes = FALSE,
       parse = TRUE, hjust = -.06, vjust = -.6, size = TXT, colour = "grey20", fill = "white",
       border.colour = "black", linewidth = .25, label.padding = unit(1.5, "pt")) +
    # One title row: protein and fluid in each strip ("NPTX2 · CSF"), no protein title above.
    facet_wrap(~fluid, nrow = 1, labeller = as_labeller(function(x) paste(g, "\u00b7", x))) +
    scale_colour_manual(values = cols, name = "Protein level") +
    scale_fill_manual(values = cols, guide = "none") +
    labs(x = "Years from baseline") + base
  # wider x padding, shared with the risk row below, so its edge numbers sit under the ticks
  if (progression) p <- p + scale_x_continuous(expand = expansion(add = RISK_PAD))
  if (!show_y) p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  font_layers(p)
}
# Number at risk under each progression panel: the whole model sample (the curves are
# adjusted predictions for tertile medians, not observed groups), at the axis breaks.
RISK_T <- c(0, 2, 4, 6); RISK_PAD <- .5
n0 <- merge(rk[years == 0, .(snp, n_risk)], ct[model == "Cox AD-CI, CU" & term == "g", .(snp, n)], by = "snp")
stopifnot(nrow(n0) == uniqueN(sv$snp), n0$n_risk == n0$n, all(RISK_T %in% rk$years))
risk_row <- function(g, show_y = TRUE) {
  # One framed row of counts under each curve panel's x title. One row, all tertiles together:
  # the curves are model predictions, not observed groups.
  p <- ggplot(rk[gene == g & years %in% RISK_T], aes(years, 0, label = n_risk)) +
    geom_text(size = TXT, colour = "grey20") +
    # empty strips with no margin: the final theme_pgd() would restore an element_blank()
    facet_wrap(~fluid, nrow = 1, labeller = as_labeller(function(x) rep("", length(x)))) +
    coord_cartesian(xlim = c(0, 6), ylim = c(-1, 1)) +
    scale_x_continuous(expand = expansion(add = RISK_PAD)) +
    labs(x = NULL, y = "N at risk") + base +
    # "N at risk" level in the y-title column
    theme(strip.text = element_text(margin = margin(0, 0, 0, 0)), panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = .5), plot.margin = margin(0, 5.5, 0, 5.5))
  if (!show_y) p <- p + theme(axis.title.y = element_blank())
  font_layers(p)
}
# ---- a, b: proteins by A/T stage -------------------------------------------------------
# Stage names only: the schematic beside a defines them.
CAT <- c(Other = "Other", Early = "Early", Core = "Core", `Late, up` = "Late, up",
         `Late, down` = "Late, down")
stopifnot(setequal(pp$category, names(CAT)), all(names(CAT) %in% pt$group))
pp[, cat := factor(CAT[category], levels = CAT)]
pt <- pt[group %in% names(CAT)][, cat := factor(CAT[group], levels = CAT)]
xlab_n <- setNames(sprintf("%s (%d)", CAT, pt[match(names(CAT), group), n]), CAT)
ymax <- max(pp$P)
# stage colour per box: Late up and down share the Late purple
cat_col <- setNames(c("grey60", STAGE_COL[c("Early", "Core", "Late", "Late")]), CAT)
pa <- ggplot(pp, aes(cat, P)) +
  geom_boxplot(aes(fill = cat), outlier.shape = NA, width = .55, linewidth = .3, colour = "grey30",
               alpha = .35) +
  geom_jitter(data = pp[category != "Other"], aes(colour = cat), width = .15, height = 0,
              size = .5, alpha = .8) +
  scale_fill_manual(values = cat_col, guide = "none") +
  scale_colour_manual(values = cat_col, guide = "none") +
  geom_text(data = pt[group != "Other"], aes(cat, ymax * 2, label = p_expr(p_P)), parse = TRUE,
            size = TXT, colour = "grey20", inherit.aes = FALSE) +
  # Above Other, in the P row: says what each group's P is tested against.
  annotate("text", x = 1, y = ymax * 2, label = 'italic(P)~"versus Other"', parse = TRUE,
           size = TXT, colour = "grey20") +
  scale_y_log10(expand = expansion(mult = c(.05, .12))) + scale_x_discrete(labels = xlab_n) +
  labs(x = NULL, y = "Clinical effect heterogeneity") + base +
  theme(panel.grid.major.x = element_blank())
# b: P within each outcome type (12f_outcome_tests.csv), Late up and down together; outlined
# cells differ from Other at P < 0.05.
OT <- c("Baseline structure", "Structural change", "Cognitive decline", "Progression")
GRP <- c("Other", "Early", "Core", "Late")
ot <- rd("12f_outcome_tests.csv")
stopifnot(setequal(ot$outcome_type, OT), setequal(ot$group, GRP), nrow(ot) == 16)
ot[, `:=`(outcome_type = factor(outcome_type, levels = rev(OT)), group = factor(group, levels = GRP),
          sig = (p_vs_other < 0.05) %in% TRUE)]
pb <- ggplot(ot, aes(group, outcome_type)) +
  geom_tile(aes(fill = median_P), colour = "white", linewidth = .6) +
  geom_tile(data = ot[sig == TRUE], fill = NA, colour = "grey15", linewidth = .4,
            width = .92, height = .88) +
  # white numbers on the darkest cells
  geom_text(aes(label = sprintf("%.1f", median_P), colour = median_P > 5), size = TXT) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_fill_gradient(low = "white", high = STAGE_COL[["Core"]], trans = "log10") +
  scale_x_discrete(position = "bottom") +   # group sizes are on the boxplot axis
  labs(x = NULL, y = NULL) + base +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(),
        legend.position = "none")   # each cell carries its value
pa <- font_layers(pa); pb <- font_layers(pb)
# How the groups were defined, drawn top to bottom beside a: the three A/T stages and the
# change(s) defining each group (vertical, so it needs no row of its own)
stg <- data.table(y = c(3, 2, 1), lab = c("A\u2212T\u2212", "A+T\u2212", "A+T+"),
                  fill = c("grey95", "grey88", "grey80"))
ps <- ggplot() +
  geom_tile(data = stg, aes(1, y, fill = fill), width = 1, height = .42, colour = NA) +
  scale_fill_identity() +
  geom_text(data = stg, aes(1, y, label = lab), size = TXT7, colour = "grey20") +
  annotate("segment", x = 1, xend = 1, y = c(2.76, 1.76), yend = c(2.24, 1.24),
           colour = STAGE_COL[c("Early", "Late")],
           linewidth = .3, arrow = arrow(length = unit(3, "pt"), type = "closed")) +
  annotate("text", x = 1.12, y = c(2.5, 1.5), label = c("Early", "Late"), hjust = 0,
           size = TXT7, colour = STAGE_COL[c("Early", "Late")]) +
  annotate("segment", x = c(1.95, 1.85, 1.85), xend = c(1.95, 1.95, 1.95), y = c(2.5, 2.5, 1.5),
           yend = c(1.5, 2.5, 1.5), colour = STAGE_COL[["Core"]], linewidth = .3) +
  annotate("text", x = 2.05, y = 2, label = "Core\n(both)", hjust = 0, lineheight = .9,
           size = TXT7, colour = STAGE_COL[["Core"]]) +
  coord_cartesian(xlim = c(.45, 2.8), ylim = c(.7, 3.3), expand = FALSE) + theme_void()
ps <- font_layers(ps)

# theme_pgd() first: the maps are wrapped below, where the outer `fig & theme_pgd()`
# does not reach, so without it they fall back to the default font.
theme_maps <- theme_pgd() + theme(plot.title = element_text(size = 7, face = "plain", hjust = .5),
  legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
  legend.key.width = unit(18, "pt"), legend.key.height = unit(5, "pt"))
# One anatomical view, identical colour limits across fluids.
atlas <- aseg; atlas$data <- atlas$data[atlas$data$side == "coronal", ]
d <- ba[gene == "NCAM1" & measure == "subcortical volume" & analysis == "cross-sectional"]
d[, label := to_aseg(region)]
lim <- max(abs(d$std_beta), na.rm = TRUE)
wanted <- unique(as.data.frame(atlas)[grepl("hippocampus|lateral ventricle",
  as.data.frame(atlas)$region), c("label", "region")]); names(wanted) <- c("label", "text")
# One name per structure (left hemisphere), not one per side.
wanted <- wanted[grepl("^Left-", wanted$label), ]
mp <- lapply(c("CSF", "Plasma"), function(fl)
  brain(d[fluid == fl & label %in% atlas$data$label], atlas, lim,
    fl, "Std. beta per SD", label_ns = FALSE, names = wanted, max_labels = 2) |>
    font_layers() |> pad_labels() + theme_maps)
# Each colour bar stands upright to the right of its own maps, its title turned 90 degrees
# (2026-11-12): under the maps, the bar row took height the maps need.
leg_right <- theme(legend.position = "right", legend.direction = "vertical",
  legend.title.position = "left", legend.title = element_text(size = 6, angle = 90, hjust = .5),
  legend.key.width = unit(5, "pt"), legend.key.height = unit(20, "pt"))
leg_sub <- cowplot::get_legend(mp[[1]] + leg_right)
subcortical <- wrap_plots(mp, nrow = 1) & theme(legend.position = "none")
# Surface views show a distinct MRI outcome: longitudinal cortical thinning.
surface_data <- ba[gene == "NCAM1" & measure == "cortical thickness" & analysis == "longitudinal"]
stopifnot(nrow(surface_data) > 0, any(surface_data$fit != "skipped"))
surface_data[, label := to_dk(region)]
surface_lim <- max(abs(surface_data$std_beta), na.rm = TRUE)
surface_plots <- lapply(c("CSF", "Plasma"), function(fl)
  brain(surface_data[fluid == fl], dk, surface_lim,
    fl, "Std. beta per SD per year", label_ns = FALSE, max_labels = 2,
    position = position_brain(hemi ~ side)) |>
    font_layers() + theme_maps)
leg_surf <- cowplot::get_legend(surface_plots[[1]] + leg_right)
# CSF and plasma side by side, each a 2 x 2 block (hemisphere x view), so c is one row.
surface <- wrap_plots(surface_plots, nrow = 1) & theme(legend.position = "none")
# Keep the imaging outcomes and their scales visibly separate.
# Block titles in their own row, so they line up although the two blocks differ in height.
head_txt <- function(x) wrap_elements(full = grid::textGrob(x, gp = grid::gpar(fontsize = 7,
  fontfamily = FIG_FAMILY)))
maps <- wrap_plots(list(head_txt("NCAM1, baseline subcortical volume"), plot_spacer(),
                        head_txt("NCAM1, longitudinal cortical thickness"), plot_spacer(),
                        wrap_elements(full = subcortical), wrap_elements(full = leg_sub),
                        wrap_elements(full = surface), wrap_elements(full = leg_surf)),
  nrow = 2, widths = c(1, 0.1, 1.15, 0.1), heights = c(0.07, 1))
# A single grid aligns cognition and progression rows, including the legend space.
# Row c: mPACC, row d: progression with its number-at-risk row; NPTX2 block left, NCAM1 right.
clinical <- list(curve(genes[1]) + labs(tag = "c"), curve(genes[2], show_y = FALSE),
                 curve(genes[1], progression = TRUE) + labs(tag = "d"),
                 curve(genes[2], progression = TRUE, show_y = FALSE),
                 risk_row(genes[1]), risk_row(genes[2], show_y = FALSE))
# axis_titles = "collect": one "Years from baseline" per row
bottom <- wrap_plots(clinical, ncol = 2, heights = c(1, 1, .16), guides = "collect",
                     axis_titles = "collect") &
  theme(legend.position = "bottom", legend.box.spacing = unit(0, "pt"))   # legend close to the risk row
# top and bottom nested, not wrapped: their tags then sit in the same column as c, and
# the final theme_pgd() reaches their text.
# one panel: b shares the group definitions drawn at its left
top <- (ps + labs(tag = "a")) + pa + pb + plot_layout(widths = c(.55, 1.55, .6))
fig <- wrap_plots(list(top, wrap_elements(full = maps) + labs(tag = "b"), bottom),
  ncol = 1, heights = c(.5, .75, 1)) &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_fig(fig & theme_pgd(), file.path(out_dir, "figures"), "12_clinical_integrated", 170)
