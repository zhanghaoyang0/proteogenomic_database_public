#!/usr/bin/env Rscript
# Plotting only, from paired effects, variant annotation and exported tissue/locus summaries.
#
#   03_variant_annotation        integrated cross-fluid effects and regulatory contexts
#     a  paired CSF against plasma effects (data prepared by 10_plot.R)
#     b  conditional GTEx support by protein-level detection class
#     c-d  ROR1 and THY1 regional associations
#     e  complete brain-region/cell-type matrix on the common protein set
# Top row: a+b (one patchwork, so panels and keys align) 40%, c 30%, d 30%; the matrix
# occupies the full second row. Panel titles only where they name a gene (c, d).
#   03_supp_single_fluid          shared and single-fluid signals (supplement)
#     a  single-fluid class shares by discovery strength (10_plot.R make_b)
#     b  VEP consequences of shared and single-fluid signals
#   (beta in Abeta+ against Abeta- moved to Supplementary Fig. 1c, 09f / lib/repl_panel.R)
#     d  GRN regional associations, the shared-regulation example
# All source observations, consequence categories and matrix rows/columns are retained.
# The Roadmap chromatin contrast (03b) is reported in the text only; its supplementary
# figure was dropped on 2026-10-29 (no contrast survives FDR).
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
OUT_DIR <- "output"
source("code/09_statistical_analysis/10_plot.R")           # fh, MAIN_COL, make_a; draws nothing when sourced
source("code/09_statistical_analysis/lib/fig_theme.R")   # theme_pgd, save_fig, pt_size, p_expr

# Only the names are used (class order and checks); no layer is coloured by these.
CLASS_COLS <- c(Shared = "grey45", `CSF-specific` = FLUID_COL[["CSF"]], `Plasma-specific` = FLUID_COL[["Plasma"]])
CLASSES <- names(CLASS_COLS)
CSQ_LV <- c("Protein-altering", "Splice region", "5' UTR", "3' UTR", "Synonymous",
            "Non-coding exon", "Up/downstream", "Intronic", "Intergenic")
CSQ_COLS <- VEP_COL   # lib/fig_theme.R
TYPE_LAB <- c(cis = "cis signals", trans = "trans signals")

bt <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), axis.ticks.length = unit(2, "pt"),
        axis.text = element_text(size = 6))

# ---- supplementary VEP consequences --------------------------------------------
csq <- fread(file.path(OUT_DIR, "03b_consequence.csv"))[set == "proteins assayed in both fluids"]
stopifnot(setequal(unique(csq$class), CLASSES), all(csq$csq %in% CSQ_LV))
# Displayed as CSF only / Plasma only: detection in one fluid, kept apart from the
# effect-specific class of panel a. The data keep their CSF-/Plasma-specific keys.
VEP_CLASS_LAB <- c(Shared = "Shared", `CSF-specific` = "CSF only", `Plasma-specific` = "Plasma only")
stopifnot(setequal(names(VEP_CLASS_LAB), CLASSES))
csq[, `:=`(csq = factor(csq, rev(CSQ_LV)), class = factor(class, rev(CLASSES), labels = VEP_CLASS_LAB[rev(CLASSES)]))]
nlab <- unique(csq[, .(class, type, n = n_signals)])
pa <- ggplot(csq, aes(100 * share, class, fill = csq)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  geom_text(data = nlab, aes(x = 102, y = class, label = n), inherit.aes = FALSE,
            hjust = 0, size = pt_size(5.5)) +
  geom_text(data = csq[csq == "Protein-altering"],
            aes(x = 100 * share / 2, label = sprintf("%.0f%%", 100 * share)),
            colour = "white", size = pt_size(5.5), fontface = "bold") +
  # Column header for the signal counts, above the first row of the cis facet only.
  geom_text(data = data.frame(type = "cis", y = length(CLASSES) + .5), aes(x = 102, y = y),
            label = "n", inherit.aes = FALSE, hjust = 0, size = pt_size(5.5), fontface = "italic") +
  facet_wrap(~ type, ncol = 1, labeller = as_labeller(TYPE_LAB)) +
  scale_fill_manual(values = CSQ_COLS, breaks = CSQ_LV, name = NULL) +
  scale_x_continuous(limits = c(0, 110), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Share of fine-mapped probability (%)", y = NULL) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  bt + theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
             legend.key.size = unit(6, "pt"),
             # plain strip titles, no frames; only the x axis line is kept
             panel.border = element_blank(), axis.line.x = element_line(linewidth = .3),
             strip.background = element_blank(),
             strip.text = element_text(size = 7, margin = margin(1.5, 0, 1.5, 0)),
             panel.spacing.y = unit(8, "pt"))

# The previous standalone regulatory figure exposes its panel definitions when sourced.
an <- list(fh = fh, pa = pa, CSQ_LV = CSQ_LV, csq = csq)
reg <- new.env(parent = globalenv())
sys.source("code/09_statistical_analysis/10d_plot_regulatory_context.R", envir = reg)
FONT <- FIG_FAMILY
theme_set(theme_compact())
invisible(cowplot::set_null_device("agg"))
update_geom_defaults("text", list(family = FONT))
update_geom_defaults("label", list(family = FONT))
fh <- copy(an$fh)
cl_order <- c("Shared (both fluids)", "Undetermined", "Concordant", "Divergent")
setorder(fh, display_class)
fh[, draw_order := match(display_class, cl_order)]
setorder(fh, draw_order)
scatter_col <- c(`Shared (both fluids)` = CLASS_COL[["shared"]], Undetermined = CLASS_COL[["unresolved"]],
                 Concordant = CLASS_COL[["present"]], Divergent = CLASS_COL[["effect-specific"]])
lim <- max(abs(c(fh$beta_csf, fh$beta_plasma))) * 1.04
f4a <- ggplot(fh, aes(beta_plasma, beta_csf)) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = .25) +
  geom_vline(xintercept = 0, colour = "grey88", linewidth = .25) +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linetype = "dashed", linewidth = .3) +
  geom_point(aes(colour = display_class), size = .65, alpha = .6, stroke = 0) +
  scale_colour_manual(values = scatter_col, breaks = c("Divergent", "Concordant", "Shared (both fluids)", "Undetermined"),
                      labels = c("Divergent", "Concordant", "Both fluids", "Undetermined"), name = NULL) +
  # 2 x 2 key: b's three-row key, not a's, then sets the height of the legend row.
  guides(colour = guide_legend(ncol = 2, override.aes = list(alpha = 1, size = 1.4))) +
  # Equal limits and a square panel (CSF and plasma effects on the same scale); patchwork
  # sizes a's column from the row height, so its bottom still lines up with b's.
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
  scale_x_continuous(breaks = c(-2, 0, 2)) + scale_y_continuous(breaks = c(-2, 0, 2)) +
  labs(tag = "a", title = NULL, x = "Plasma pQTL effect (β)", y = "CSF pQTL effect (β)") +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.8), aspect.ratio = 1)

# Retain every consequence category, soften context colours around the PAV highlight.
vep_soft <- VEP_SOFT_COL
vep_supp <- an$pa + theme_compact() +
  scale_fill_manual(values = vep_soft, breaks = an$CSQ_LV, name = NULL) +
  labs(tag = NULL, title = NULL, subtitle = NULL,
       x = "Fine-mapped probability (%)") +
  # Full width (2026-11-12): the key runs in one row under the bars.
  guides(fill = guide_legend(nrow = 1,
                             label.position = "right")) +
  theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text = element_text(size = 6), panel.spacing.y = unit(3, "pt"),
        legend.text = element_text(size = 6), legend.key.size = unit(8, "pt"),
        legend.justification = "center",
        plot.margin = margin(3, 1, 3, 2))

f4b <- reg$bridge + theme_compact() +
  labs(tag = "b", title = NULL, subtitle = NULL,
       y = "Proportion of proteins") +
  scale_x_discrete(labels = c(`Plasma only` = "Plasma\nonly", `Both fluids` = "Both\nfluids", `CSF only` = "CSF\nonly")) +
  # Italic n, as in the text.
  geom_text(data = reg$n_by_class, aes(x = class, y = 1.04, label = paste0("italic(n)~'= ", n_proteins, "'")),
            inherit.aes = FALSE, family = FONT, size = pt_size(6), parse = TRUE) +
  coord_cartesian(ylim = c(0, 1.12), clip = "off") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom", axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 5.8), legend.text = element_text(size = 5.8))

matrix_heat <- reg$heat + theme_compact() +
  labs(tag = "e", title = NULL, subtitle = NULL) +
  guides(fill = "none") +
  theme(axis.line = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x.top = element_text(size = 5.5, angle = 90, hjust = 0, vjust = .5, face = "italic"),
        axis.text.y = element_text(size = 6.2),
        strip.text.y.left = element_text(size = 5.7, angle = 90),
        strip.placement = "outside", panel.spacing.y = unit(1.5, "pt"),
        plot.margin = margin(1, 4, 9, 3))
matrix_bars <- reg$bars + theme_compact() +
  labs(title = NULL, x = NULL, y = NULL) +
  scale_fill_manual(values = reg$STATES, name = NULL,
                    breaks = c("Both fluids", "CSF only", "Plasma only")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.text = element_blank(), panel.spacing.y = unit(1.5, "pt"),
        # Key just above the bars (inside the plot, anchored on the panel's top-right
        # corner), not in the top-axis band of the gene labels where it floated.
        legend.position = "inside", legend.position.inside = c(1, 1),
        legend.justification.inside = c(1, 0), legend.text = element_text(size = 5.5),
        legend.margin = margin(0, 0, 3, 0), legend.background = element_blank(),
        plot.margin = margin(1, 8, 9, 0))
f4e <- (matrix_heat | matrix_bars) + plot_layout(widths = c(4.3, 1))

region <- function(g, tag) {
  d <- copy(reg$lt[gene == g])
  d[, track := factor(track, levels = c(setdiff(unique(track), c("CSF pQTL", "Plasma pQTL")),
                                         "CSF pQTL", "Plasma pQTL"))]
  setorder(d, r2, na.last = FALSE)
  lp <- unique(d[lead == TRUE, pos])
  stopifnot(length(lp) == 1L)
  floor_data <- d[, .(pos = lp, mlp = 10), by = track]
  xp <- range(d$pos)/1e6
  ggplot(d, aes(pos/1e6, mlp)) +
    geom_blank(data = floor_data) +
    geom_hline(yintercept = -log10(5e-8), colour = "grey65", linetype = "dotted", linewidth = .25) +
    geom_vline(xintercept = lp/1e6, colour = LEAD_LINE_COL, linewidth = .25, alpha = .65) +
    # show.legend = TRUE: draw every LD key, also bins without a point in this locus
    geom_point(aes(colour = ld), size = .68, alpha = .8, stroke = 0, show.legend = TRUE) +
    geom_point(data = d[lead == TRUE & p < 5e-8], shape = 23, size = 1.55,
                fill = LEAD_COL, colour = "grey20", stroke = .25) +
    scale_colour_manual(values = setNames(LD_COLS, LD_LABELS), breaks = LD_LABELS,
                        name = "LD r²", na.value = "grey78", drop = FALSE) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 1.5))) +
    # Track titles centred above each track and naming the gene, as in Fig. 3b and 5b,c
    # ("CSF ROR1 pQTL"), so there is no separate gene title above the column; the gene is
    # italic in the eQTL track (expression), plain in the pQTL tracks (protein).
    facet_wrap(~ track, ncol = 1, scales = "free_y", labeller = as_labeller(function(x) ifelse(
      grepl("eQTL$", x), sprintf('%s~italic("%s")~eQTL', sub(" eQTL$", "", x), g),
      sprintf('%s~"%s"~pQTL', sub(" pQTL$", "", x), g)), default = label_parsed)) +
    scale_x_continuous(limits = xp, breaks = scales::breaks_pretty(4), expand = c(0, 0)) +
    scale_y_continuous(breaks = function(l) pretty(c(0, l[2]*.8), 2),
                       expand = expansion(mult = c(.04, .1))) +
    labs(tag = tag, title = NULL,
         x = paste0("Chromosome ", d$chr[1], " (Mb)"), y = "−log₁₀ P") +
    theme(strip.text = element_text(size = 6, hjust = .5, margin = margin(1, 0, 1, 0)),
          strip.background = element_blank(), panel.spacing.y = unit(2.5, "pt"),
          legend.position = "bottom")
}
# The track titles name the gene; the colocalization interpretation is in the caption.
f4c <- region("ROR1", "c")
f4d <- region("THY1", "d")
# a-d in one patchwork row, so the four panels share a top edge, their x axes share a
# bottom edge, and the a, b and LD keys
# share one legend row with their first lines level. Tags at the top left of each plot,
# with no top margin, so all four sit at one height.
row_th <- theme(plot.tag = element_text(size = 8, face = "bold", family = FONT),
                plot.tag.position = c(0.02, 1), plot.tag.location = "plot",
                plot.margin = margin(0, 3, 3, 3),
                legend.justification = c(0.5, 1))
# The LD key belongs to c and d and is drawn once, centred under the pair: it is d's
# legend, anchored at d's left edge (the c/d boundary) and pulled left by half its own
# width, measured here, so it stays in the a/b legend row.
ld_w <- grid::convertWidth(sum(cowplot::get_legend(f4c + theme_pgd())$widths), "pt",
                           valueOnly = TRUE)
row4top <- (f4a + row_th) + (f4b + row_th) +
  (f4c + row_th + guides(colour = "none")) +
  (f4d + row_th + theme(legend.position = "bottom", legend.justification = c(0, 1),
                        legend.margin = margin(0, 0, 0, -ld_w / 2))) +
  plot_layout(nrow = 1, widths = c(2.85, 1.75, 2.7, 2.7))
fig4 <- wrap_plots(wrap_elements(full = row4top), wrap_elements(full = f4e),
                  ncol = 1, heights = c(61, 59))
save_fig(fig4, file.path(OUT_DIR, "figures"), "03_variant_annotation", 120, dpi = 350, raster_device = "ragg")
# Supplementary figure: single-fluid classes by discovery strength (a, 40%) beside the
# VEP consequences (b, 60%), one row (2026-11-13). The amyloid-stratified effects moved to
# Supplementary Fig. 1c (2026-11-12), with the other checks of effect robustness.
supp_tag <- theme(plot.tag = element_text(size = 8, face = "bold", family = FONT))
# Open axes, no panel boxes, as in c, d and the main figures (2026-11-12).
open_axes <- theme(panel.border = element_blank(),
                   axis.line = element_line(linewidth = .3, colour = "grey30"))
s_a <- make_b() + labs(subtitle = NULL, tag = "a") + theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1)) + theme_pgd() + open_axes + supp_tag
# b's nine consequence classes need two rows in 60% of the width.
s_b <- vep_supp + labs(tag = "b") + guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + supp_tag
SUPP_SF_MM <- 85
supp_sf <- (s_a | s_b) + plot_layout(widths = c(2, 3))
save_fig(supp_sf, file.path(OUT_DIR, "figures"), "03_supp_single_fluid", SUPP_SF_MM, dpi = 350,
         raster_device = "ragg", max_h_mm = FIG_H_MAX_SUPP_MM)
message("03c_plot_annotation.R done")
