#!/usr/bin/env Rscript
# Supplementary figure for section 12 (output/12_fluid_clinical.png), from the tables that
# code/09_statistical_analysis/local/12_fluid_clinical.R and 12b_fluid_clinical_summary.R write. Plotting only.
#   a  cross-fluid genetic effect difference (G) against cross-fluid clinical effect
#      heterogeneity (P), one point per protein
#   b  FDR-significant protein-outcome associations of CSF and of plasma levels, per outcome
#      (a count grid)
# Usage (repo root): Rscript code/09_statistical_analysis/12c_plot.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
source("code/09_statistical_analysis/lib/fig_theme.R")   # theme_pgd, save_fig, pt_size, p_expr
OUT_DIR <- "output"
pr <- fread(file.path(OUT_DIR, "12b_protein_summary.csv"))
tests <- fread(file.path(OUT_DIR, "12b_tests.csv"))
r <- fread(file.path(OUT_DIR, "12_fluid_clinical.csv"))

# Open axes, no panel boxes, as in the other figures (2026-11-12).
theme_fig <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line(linewidth = .3, colour = "grey30"))

prim <- tests[test == "Spearman(G, P)" & P_def == "P"]
stopifnot(nrow(prim) == 1)
plain_num <- function(b) format(b, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)  # 0.001, not 1e-03
# One grey: the result is the overall (absent) correlation, not a contrast between the
# paired-fluid classes (2026-11-13; the classes were coloured before).
pa <- ggplot(pr, aes(G, P)) +
  geom_point(size = 0.9, alpha = 0.6, colour = "grey40") +
  scale_x_log10(labels = plain_num) + scale_y_log10(labels = plain_num) +
  # Top left is empty (few proteins have a small G and a large P); top right held points.
  # One boxed line, framed like the key (lib/fig_theme.R legend_box).
  # A true minus sign: sprintf's "-" prints as a hyphen.
  annotate("label", x = 2 * min(pr$G), y = Inf, hjust = 0, vjust = 1.2, size = pt_size(6), parse = TRUE,
           family = FIG_FAMILY, fill = "white", border.colour = "black", linewidth = 0.25,
           label.padding = unit(1.5, "pt"),
           label = sprintf('"Spearman"~rho=="%s"*";"~~%s',
                           sub("^-", "\u2212", sprintf("%.2f", prim$estimate)), p_expr(prim$p))) +
  # Plain axis titles; the definitions (|z| of the CSF-plasma lead-pQTL difference, mean
  # z² of the association differences over 20 outcomes) are in the caption.
  # Plain-language titles (2026-11-13); the statistics (|z|, mean z²) are in the caption.
  labs(x = "Genetic effect difference between fluids",
       y = "Clinical association difference between fluids") +
  theme_fig

OUT_LAB <- c(HV = "Hippocampal volume", ILV = "Inferior lateral ventricle", mPACC = "mPACC",
             `AD-CI` = "Progression to Aβ+ MCI or dementia", ct_adsign_surfwt_lr = "AD-signature thickness",
             ct_frontal_surfwt_lr = "Frontal thickness", ct_lateralparietal_surfwt_lr = "Lateral parietal thickness",
             ct_lateraltemporal_surfwt_lr = "Lateral temporal thickness", ct_medialparietal_surfwt_lr = "Medial parietal thickness",
             ct_medialtemporal_surfwt_lr = "Medial temporal thickness", ct_occipital_surfwt_lr = "Occipital thickness")
cnt <- fread(file.path(OUT_DIR, "12b_single_fluid_counts.csv"))[, .(model, analysis, outcome, n = n_fdr)]   # 12b
cnt[, fluid := fifelse(model == "CSF only", "CSF", "Plasma")]
if (length(setdiff(cnt$outcome, names(OUT_LAB)))) stop("outcome without a label: ", setdiff(cnt$outcome, names(OUT_LAB))[1])
cnt[, outcome := factor(OUT_LAB[outcome], levels = rev(OUT_LAB))]
cnt[, analysis := factor(analysis, levels = c("baseline", "longitudinal", "progression"),
                         labels = c("Baseline", "Change", "Progression"))]
# A table-like grid (2026-11-13): one row per outcome, a CSF and a plasma column, the
# count in each cell and its shade on a square-root scale, so the 399 of inferior lateral
# ventricle does not flatten the single-digit cells and a zero reads as an empty cell.
cnt[, shade := sqrt(n) / sqrt(max(n))]
pb <- ggplot(cnt, aes(fluid, outcome)) +
  geom_tile(aes(fill = fluid, alpha = shade), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = n, colour = shade > 0.55), size = pt_size(6), family = FIG_FAMILY) +
  scale_fill_manual(values = FLUID_COL, guide = "none") +
  scale_alpha_continuous(range = c(0.04, 1), limits = c(0, 1), guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  # Block names above their rows.
  ggforce::facet_col(~ analysis, scales = "free_y", space = "free", strip.position = "top") +
  labs(x = NULL, y = NULL,
       title = sprintf("Proteins associated at FDR < 0.05 (of %d)", uniqueN(r$ensembl_gene_id))) +
  theme_fig + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    panel.grid = element_blank(), strip.text = element_text(hjust = 0),
                    plot.title = element_text(hjust = 1, size = 7), plot.title.position = "panel")

save_fig(pa + pb + plot_layout(widths = c(1.5, 1)) + plot_annotation(tag_levels = "a") &
           theme_pgd(), file.path(OUT_DIR, "figures"), "12_fluid_clinical", 88)
