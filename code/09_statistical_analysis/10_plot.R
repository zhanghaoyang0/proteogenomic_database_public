#!/usr/bin/env Rscript
# What the paired individual-level data add. Plotting only, from the summary tables the
# Mac scripts in code/09_statistical_analysis/local/ write (10a, 10b), so it runs anywhere without the
# individual data.
# make_a: CSF against plasma effect, drawn by 03c_plot_annotation.R (sourced).
# make_b / make_ab: the two robustness checks of the paired classes, drawn as a and b of
# 03_supp_single_fluid.png by 03c_plot_annotation.R (sourced):
#   make_b   single-fluid class shares by discovery P and fluid (winner's curse)
#   make_ab  beta in Abeta+ against Abeta-, by fluid (10b)
# APOE +- 1 Mb is left out of everything, as in the 10a/10b summaries.
# Usage:
#   Rscript code/09_statistical_analysis/10_plot.R [out_dir]   (default output/)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

source("code/09_statistical_analysis/lib/fig_theme.R")   # theme_pgd, save_fig, pt_size
IN_DIR <- "output"
args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[1] else IN_DIR   # the tables are always read from output/

read_out <- function(name) {
  f <- file.path(IN_DIR, name)
  if (!file.exists(f)) stop("Missing ", f, ": run the script in code/09_statistical_analysis/local/ that writes it.")
  fread(f)
}

FDR_MAX <- 0.05
# Validated with the dataviz palette checker (CVD and normal-vision separation pass).
# Shared is light grey: in panel a it is the background the single-fluid classes are
# read against, and panel b has no shared signals.
# CLASS_COL (the class colours) is shared, in lib/fig_theme.R.
CLASS_LAB <- c(shared = "Shared (both fluids)", present = "Concordant",
               `effect-specific` = "Divergent", unresolved = "Undetermined")
# Legend order and panel-a draw order (later on top): the classes the figure is about last.
DRAW_ORDER <- c("shared", "unresolved", "present", "effect-specific")

theme_fig <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
        strip.background = element_blank(),
        strip.text = element_text(size = 8.5),
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 7.5),
        plot.tag = element_text(face = "bold", size = 11))
# make_a above keeps theme_fig: its look is set by 03c, which draws it. The supplement
# is drawn at print size (183 mm), so its text is 5-7 pt.
theme_supp <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
        strip.background = element_blank())

# ---- data -------------------------------------------------------------------
fh <- read_out("10a_fluid_heterogeneity.csv")[apoe_region == FALSE]
p_lead <- read_out("S1_pqtl_finemapping.csv.gz")[is_cs_lead == TRUE, .(cs_key, p_gwas = p)]
n0 <- nrow(fh)
fh <- merge(fh, p_lead, by = "cs_key")
# A single-fluid signal either has a same-direction effect in the other fluid or not
# (2026-10-27): 10a's attenuated and detection-limited are merged into "present",
# because their split compared effect sizes across fluids measured with different
# precision, and 111 "attenuated" signals had the larger effect in the other fluid.
fh[class %in% c("attenuated", "detection-limited"), class := "present"]
stopifnot(nrow(fh) == n0, all(fh$class %in% names(CLASS_COL)))
fh[, class := factor(class, levels = names(CLASS_COL))]

cx <- read_out("10b_pqtl_context.csv")[modifier == "abeta" & apoe_region == FALSE &
                                       !is.na(beta_pos) & !is.na(beta_neg)]
cx[, fluid := factor(fluid, levels = c("CSF", "Plasma"))]

# ---- a: CSF against plasma effect -------------------------------------------
MAIN_COL <- setNames(CLASS_COL[c("effect-specific", "present", "shared", "unresolved")],
                     c("Divergent", "Concordant", "Shared (both fluids)", "Undetermined"))
fh[, display_class := factor(fcase(class == "effect-specific", "Divergent",
  class == "present", "Concordant", class == "shared", "Shared (both fluids)",
  default = "Undetermined"), levels = names(MAIN_COL))]
make_a <- function(legend_classes = names(MAIN_COL)) {
  lim <- max(abs(c(fh$beta_csf, fh$beta_plasma)), na.rm = TRUE) * c(-1, 1)
  d <- fh[order(match(display_class, c("Undetermined", "Shared (both fluids)",
                                      "Concordant", "Divergent")))]
  ggplot(d, aes(beta_plasma, beta_csf)) +
    geom_point(aes(colour = display_class), size = 1.15, alpha = .58, stroke = 0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey20", linewidth = .4) +
    scale_colour_manual(values = MAIN_COL, breaks = legend_classes, name = NULL, drop = FALSE) +
    coord_fixed(xlim = lim, ylim = lim, ratio = 1) +
    labs(x = expression(beta[plasma]), y = expression(beta[CSF])) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 3))) +
    theme_fig + theme(axis.title = element_text(size = 13), axis.text = element_text(size = 10))
}
# ---- b: class shares by discovery P -----------------------------------------
make_b <- function() {
  one <- fh[class != "shared"]
  one[, p_bin := factor(fcase(p_gwas > 1e-10, "weak", p_gwas > 1e-20, "mid", default = "strong"),
                        levels = c("weak", "mid", "strong"))]
  d <- one[, .N, by = .(disc_fluid, p_bin, class)]
  d[, share := N / sum(N), by = .(disc_fluid, p_bin)]
  n_lab <- one[, .N, by = .(disc_fluid, p_bin)]
  # Divergent share printed on its segment, which sits on top of unresolved
  # (ggplot stacks the first factor level on top, so unresolved is at the bottom).
  es_lab <- merge(d[class == "effect-specific", .(disc_fluid, p_bin, share)],
                  d[class == "unresolved", .(disc_fluid, p_bin, base = share)],
                  by = c("disc_fluid", "p_bin"), all.x = TRUE)
  stopifnot(nrow(es_lab) == nrow(n_lab))
  es_lab[is.na(base), base := 0]
  es_pct <- 100 * mean(one$class == "effect-specific")
  ggplot(d, aes(p_bin, share, fill = class)) +
    # No legend of its own: panel a's, which also has shared, serves both panels.
    geom_col(width = 0.75, colour = "white", linewidth = 0.3, show.legend = TRUE) +
    geom_text(data = es_lab, aes(p_bin, base + share / 2, label = sprintf("%.0f%%", 100 * share)),
              inherit.aes = FALSE, size = pt_size(6), colour = "white", fontface = "bold", family = FIG_FAMILY) +
    geom_text(data = n_lab, aes(p_bin, 1.04, label = N), inherit.aes = FALSE,
              size = pt_size(6), colour = "grey30", family = FIG_FAMILY) +
    # Strips read "CSF" / "Plasma" as in panel b; the y title says they are discovery fluids.
    facet_wrap(~ disc_fluid) +
    scale_fill_manual(values = CLASS_COL, labels = CLASS_LAB, limits = rev(setdiff(DRAW_ORDER, "shared")),
                      name = NULL, drop = FALSE) +
    # Thresholds in the caption: weak P > 1e-10, moderate 1e-20 to 1e-10, strong P < 1e-20.
    scale_x_discrete(labels = c(weak = "Weak", mid = "Moderate", strong = "Strong")) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%"), breaks = c(0, 0.5, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, 1.08), clip = "off") +
    labs(x = "Discovery signal strength", y = "Single-fluid signals\nby discovery fluid",
         subtitle = sprintf("%.0f%% of %s single-fluid signals are divergent",
                            es_pct, format(nrow(one), big.mark = ","))) +
    theme_supp + theme(panel.grid.major.x = element_blank(),
                       plot.subtitle = element_text(size = 7, face = "bold"))
}

# ---- Supplementary Fig.: Abeta+ against Abeta- --------------------------------
make_ab <- function() {
  # Aligned to the GWAS direction, so a weaker effect in Abeta+ falls below the diagonal.
  # r is on the signed effects, as in 10b and the main text, not on the aligned points
  # (aligning removes the sign agreement: r 0.93 / 0.96 instead of 0.98 / 0.99).
  r_lab <- cx[, .(r = cor(beta_neg, beta_pos), n = .N), by = fluid]
  r_lab[, lab := sprintf("italic(r) == %.2f", r)]
  hit <- cx[fdr_int < FDR_MAX]
  stopifnot(nrow(hit) >= 1)
  lim <- range(c(cx$beta_neg_al, cx$beta_pos_al))
  ggplot(cx, aes(beta_neg_al, beta_pos_al)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_point(size = 0.5, alpha = 0.5, stroke = 0, colour = "grey35") +
    geom_point(data = hit, colour = CLASS_COL[["effect-specific"]], size = 1.8) +
    geom_text_repel(data = hit, aes(label = gene), fontface = "italic", size = pt_size(6), family = FIG_FAMILY,
                    nudge_x = 0.4, nudge_y = -0.3, segment.size = 0.2, min.segment.length = 0) +
    geom_text(data = r_lab, aes(x = -Inf, y = Inf, label = lab), parse = TRUE,
              hjust = -0.15, vjust = 1.4, size = pt_size(6), colour = "grey20", family = FIG_FAMILY) +
    facet_wrap(~ fluid, nrow = 1) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = expression(beta ~ "in A" * beta * "-negative"),
         y = expression(beta ~ "in A" * beta * "-positive")) +
    theme_supp
}

# Sourced by 03c, which draws every figure from these builders (03_variant_annotation a,
# 03_supp_single_fluid a-b); nothing is drawn here.
