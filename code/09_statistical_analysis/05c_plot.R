#!/usr/bin/env Rscript
# 05_protein_disease.png and 05_supp_disease.png.
# Plotting only; every table read here is written by 05a_disease_data.R,
# 05b_case_data.R, 05b_drug_support.py, 05d_link_matrix.R and 05f_locus_data.R.
#
# Main figure since 2026-10-27 (assembled at the end of this file):
#   panel a    neurological link matrix (05d_link_matrix.R, sourced)
#   panels b-c CD58 - multiple sclerosis and GRN - Alzheimer disease locus panels
#              (05g_locus_plot.R, sourced)
#   panel d    drug evidence across CSF-only, plasma-only and shared genetic support
# Built here but no longer in the main figure: the category dumbbell (now
# 05_supp_disease.png) and the approved-drug table (named in the text). The forest of
# the six unreported CSF links was dropped on 2026-11-13; Supplementary Table 4 flags
# them (no_prior_report).
#
# A module x disease-category enrichment of the joint proteins was tested as a network
# panel and dropped: nothing survives multiple testing (best CSF cell 4/75 Digestive
# in M4, p = 0.028, BH q = 1.0; best plasma module M1 OR 1.49, p = 0.053 on 368/780
# genes). With 57 CSF and 113 plasma joint genes spread over 9 and 6 modules x 15
# categories, most cells are empty. The network figures stay separate.
#
# Only the joint MR + coloc evidence class is shown. MR alone dominates the counts
# (it tests every pair, while coloc is only attempted where there is a signal to
# colocalize), so a figure showing all three classes side by side is read as "MR
# finds more", which is a property of the two methods rather than of the biology.
# The joint class is also the one the text argues from.
#
# Output:
#   05_protein_disease.png, 05_supp_disease.png
# Usage:
#   Rscript code/09_statistical_analysis/05c_plot.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source("code/09_statistical_analysis/lib/fig_theme.R")   # theme_pgd, save_fig, pt_size
OUT_DIR <- "output"

COLS <- FLUID_COL
FLUIDS <- c("CSF", "Plasma")

out_dis_png <- file.path(OUT_DIR, "figures", "05_protein_disease.png")

s <- fread(file.path(OUT_DIR, "05a_disease_category_summary.csv"))
# Both analyses written by 05a_disease_data.R, drawn as the two sub-panels of the category figure:
# "All testable" is each fluid's discovery yield over its own instrumented proteins,
# and is not a fair comparison because the two denominators hold different proteins;
# "Head-to-head" restricts both fluids to the proteins testable in both within that
# category, so the denominator is one shared protein set (identical on both sides) and
# the disease GWAS is the same. See the comment above summarize_evidence() in
# 05a_disease_data.R.
d <- s[evidence == "MR + coloc",
       .(fluid, disease_category, analysis,
         n = n_unique_proteins, n_tested = n_tested_genes)]
d[, analysis := factor(analysis, levels = c("All testable", "Head-to-head"))]
# Every category must appear in both fluids, or a dumbbell would lose one end.
d <- merge(CJ(fluid = FLUIDS, disease_category = unique(d$disease_category),
              analysis = levels(d$analysis)),
           d, by = c("fluid", "disease_category", "analysis"), all.x = TRUE)
d[, analysis := factor(analysis, levels = c("All testable", "Head-to-head"))]
d[is.na(n), n := 0L]
d[is.na(n_tested), n_tested := 0L]
d[, rate := fifelse(n_tested > 0, 100 * n / n_tested, NA_real_)]

# Ordered by the mean of the two rates, ascending, so the strongest categories sit at
# the top of a rate axis. Mean rather than max: a category carried by one fluid alone
# (Psychiatric, 4.5% CSF against 0% plasma) is not a strong category, and max would
# promote it above categories where both fluids agree.
# Row order is set once, from the inclusive analysis, so the two sub-panels line up.
ord <- d[analysis == "All testable",
         .(m = mean(rate, na.rm = TRUE)), by = disease_category][order(m)]$disease_category
# Categories with fewer than MIN_TESTED proteins tested in either fluid are left out of
# the figure (rates such as 0/2 carry no information); the text does not name them either.
MIN_TESTED <- 10
small_cat <- d[analysis == "All testable", .(n = min(n_tested)), by = disease_category][n < MIN_TESTED, disease_category]
message("Fig. 5a: categories not drawn (< ", MIN_TESTED, " tested proteins in a fluid): ",
        paste(small_cat, collapse = ", "))
d <- d[!disease_category %in% small_cat]
# The order computed above was never applied before 2026-10-24 (rows came out
# reverse-alphabetical); it is now the factor order.
d[, disease_category := factor(disease_category, levels = setdiff(ord, small_cat))]

d[, fluid := factor(fluid, levels = FLUIDS)]

base_theme <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    axis.title.x = element_text(size = 7),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -2),
    # house style: plain strip titles without the grey box, plain centred plot titles
    strip.background = element_blank(),
    plot.title = element_text(face = "plain", hjust = 0.5)
  )

# ---- Supplementary: joint evidence by disease category ------------------------------
# One segment per category joining the two fluids: the gap is the comparison the
# panel is making, so it is drawn rather than left to the eye.
seg <- dcast(d, disease_category + analysis ~ fluid, value.var = c("n", "rate"))
seg_plot <- seg[is.finite(rate_CSF) & is.finite(rate_Plasma)]
seg_plot[, disease_category := factor(disease_category, levels = levels(d$disease_category))]

# x is the rate, not the raw count: plasma has more genes tested in most categories,
# so raw counts tend to put plasma to the right of CSF for a reason that is not
# biological. The count is not lost -- it is printed at each point, so the panel
# carries both without a second copy of itself.
#
# Labels are placed by fluid, not by value: CSF always left, plasma always right.
# Attaching them to the points cannot honour that rule -- in the rows where CSF sits
# to the right of plasma (Digestive, Neoplasms, Respiratory) both labels then point
# inward and land on top of each other. So they go in two fixed columns just outside
# the data range instead, read like a small table beside the plot. Each is written as
# count over the genes with both tests attempted, so the denominator the percentage is
# taken over is visible and the number cannot be misread as a total.
# Label columns are placed per sub-panel, because the two analyses have different
# rate ranges and each gets its own x scale.
d[, XMAX := max(rate, na.rm = TRUE), by = analysis]
d[, `:=`(PAD = 0.05 * XMAX, LAB_W = 0.22 * XMAX)]
d[, `:=`(lab_x = fifelse(fluid == "CSF", -PAD, XMAX + PAD),
         lab_hjust = fifelse(fluid == "CSF", 1, 0),
         lab = fifelse(n_tested > 0, sprintf("%d/%d", n, n_tested), "\u2014"))]

p_cat <- ggplot(d, aes(rate, disease_category)) +
  geom_segment(data = seg_plot, inherit.aes = FALSE,
               aes(x = rate_CSF, xend = rate_Plasma,
                   y = disease_category, yend = disease_category),
               colour = "grey70", linewidth = 0.5) +
  geom_point(data = d[is.finite(rate)], aes(colour = fluid), size = 1.5) +
  geom_text(aes(x = lab_x, label = lab, colour = fluid, hjust = lab_hjust),
            size = pt_size(5.5), show.legend = FALSE, family = FIG_FAMILY) +
  scale_colour_manual(values = COLS, name = NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(4),
                     expand = expansion(mult = 0.26)) +
  facet_wrap(~ analysis, nrow = 1, scales = "free_x") +
  labs(x = "Proteins supported by MR and colocalization (% of tested proteins)",
       y = NULL) +
  base_theme +
  # Key centred under the panel, as in the other supplementary figures (2026-11-12).
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        axis.ticks.length = unit(2, "pt"),
        axis.text.x = element_text(margin = margin(t = 1)))

# ---- panel d: drug evidence by fluid-support class (05b_drug_support.py) --------
# These classes partition the jointly supported gene-disease pairs, unlike the former
# overlapping "CSF"/"Plasma" bars. Counts and percentages are descriptive: the panel
# shows whether a target has a drug for the same disease, not an enrichment test.
db <- fread(file.path(OUT_DIR, "05b_drug_by_fluid_support.csv"))
support_levels <- c("CSF-supported only", "Plasma-supported only", "Supported in both")
db[, `:=`(support_class = factor(support_class, levels = rev(support_levels)),
          drug_evidence = factor(drug_evidence,
                                 levels = c("Clinical-stage drug", "Approved drug")))]
if (!setequal(levels(db$support_class), rev(support_levels)) ||
    !setequal(levels(db$drug_evidence), c("Clinical-stage drug", "Approved drug")))
  stop("unexpected drug-support classes")
# One 100% bar per class, in three parts: pairs with a drug for the same disease (Phase I or
# later), pairs whose gene is the target of a drug approved for another disease only (the
# Results' repurposing count), and the rest. Both drug parts come from the pair table
# (05b_drug_supported.csv, the supported pairs whose gene is a drug target) and must agree
# with the class totals of 05b_drug_by_fluid_support.csv.
dw <- dcast(db, support_class ~ drug_evidence, value.var = c("n_with_drug", "n_pairs"))
setnames(dw, c("support_class", "n_clin", "n_appr", "n", "n_dup"))
stopifnot(dw[, all(n == n_dup & n_appr <= n_clin)])
CLS_LAB <- c(`CSF-supported only` = "CSF only", `Plasma-supported only` = "Plasma only",
             `Supported in both` = "Both")
ds <- fread(file.path(OUT_DIR, "05b_drug_supported.csv"))[supported == TRUE]
ds[, support_class := fifelse(supported_csf & !supported_plasma, "CSF-supported only",
                      fifelse(!supported_csf & supported_plasma, "Plasma-supported only",
                              "Supported in both"))]
ds[, part := fcase(drug_stage >= 3L, "Same disease",
                   approved_drugs_any_indication != "", "Other disease", default = NA_character_)]
dc <- ds[!is.na(part), .N, by = .(support_class, part)]
dc <- dcast(dc, support_class ~ part, value.var = "N", fill = 0L)
dw <- merge(dw, dc, by = "support_class", all.x = TRUE)
for (k in c("Same disease", "Other disease")) set(dw, which(is.na(dw[[k]])), k, 0L)
stopifnot(dw[, all(`Same disease` == n_clin)], dw[, all(`Same disease` + `Other disease` <= n)])
dw[, Neither := n - `Same disease` - `Other disease`]
dw[, support_label := factor(CLS_LAB[as.character(support_class)], levels = rev(CLS_LAB))]
seg <- melt(dw, id.vars = c("support_label", "n"), measure.vars = names(DRUG_ALPHA),
            variable.name = "part", value.name = "k")
seg[, `:=`(part = factor(part, levels = rev(names(DRUG_ALPHA))), pct = 100 * k / n)]
# segment counts inside the bar where they fit; the class total at the bar's end
seg <- seg[order(support_label, -as.integer(part))][, mid := cumsum(pct) - pct / 2, by = support_label]
ROW_COL <- c(`CSF only` = FLUID_COL[["CSF"]], `Plasma only` = FLUID_COL[["Plasma"]], Both = "grey35")
p_drug <- ggplot(seg, aes(pct, support_label, fill = support_label, alpha = part)) +
  geom_col(width = 0.62) +
  geom_text(data = seg[part != "Neither" & pct >= 2], aes(mid, support_label, label = k,
                                                          colour = part == "Same disease"),
            inherit.aes = FALSE, size = pt_size(5.5), family = FIG_FAMILY, show.legend = FALSE) +
  geom_text(data = dw, aes(100, support_label, label = paste0("n = ", n)), inherit.aes = FALSE,
            hjust = -0.1, size = pt_size(5.5), family = FIG_FAMILY) +
  scale_fill_manual(values = ROW_COL, guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none") +
  # The key shows the opacity steps in a neutral grey; each bar carries its class colour.
  scale_alpha_manual(values = DRUG_ALPHA, breaks = names(DRUG_ALPHA),
                     labels = c(`Same disease` = "Drug for the same disease (Phase I+)",
                                `Other disease` = "Approved drug for another disease only",
                                Neither = "Neither"),
                     name = NULL, guide = guide_legend(override.aes = list(fill = "grey30"))) +
  scale_x_continuous(limits = c(0, 112), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0))) +
  # two-line, rotated labels keep the column's left axis narrow, so b's y titles sit by
  # its axis
  scale_y_discrete(labels = function(x) sub(" ", "\n", x)) +
  # "associations", as in the Results; "jointly supported" read as "supported in both fluids"
  labs(x = "Protein\u2013disease associations (%)", y = NULL) +
  guides(alpha = guide_legend(nrow = 1, override.aes = list(fill = "grey30"))) +
  base_theme +
  theme(panel.grid.major = element_blank(), axis.ticks.length = unit(2, "pt"),
        axis.text.y = element_text(size = 6, angle = 90, hjust = 0.5, lineheight = 0.9),
        legend.position = "bottom", legend.justification = "center",
        legend.margin = margin(t = -1))




# ---- assemble ----------------------------------------------------------------
# Redesigned 2026-10-27: a, the neurological protein-disease link matrix (05d); b-c, the
# CD58 (CSF only) and GRN (both fluids) locus panels (05g); d, drug support. The
# figure's argument is the fluid difference, then two mechanisms, then that supported
# pairs are worth something. The novel-link forest above is its own supplementary
# figure; the approved-drug table is named in the text.
load_panels <- function(file) {
  e <- new.env()
  e$PANEL_EXPORT_ONLY <- TRUE
  sys.source(file.path("code/09_statistical_analysis", file), envir = e)
  e
}
p_matrix <- load_panels("05d_link_matrix.R")$p_matrix
loc <- load_panels("05g_locus_plot.R")
p_b <- loc$case_column("CD58")
p_c <- loc$case_column("GRN")
p_b[[1]] <- p_b[[1]] + labs(tag = "b")
p_c[[1]] <- p_c[[1]] + labs(tag = "c")
# the identical r2 scales of the locus panels collect into one legend under them
# wrap_plots() nests each locus column whole; "/" would merge their seven rows into
# the outer layout
# b and c side by side, so each locus keeps the full height for its seven tracks.
p_loci_new <- wrap_plots(p_b, p_c, nrow = 1) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 6),
        legend.title = element_text(size = 6), legend.key.size = unit(6, "pt"))
# d kept out of the locus columns' axis alignment (wrap_elements; free() fails on this nesting): its two-line class labels otherwise set
# the y-label width of every track above and push their titles off the axis
p_right <- wrap_plots(p_loci_new, wrap_elements(full = p_drug + theme_pgd() + theme(plot.margin = margin(10, 2, 1, 6))) +
                        labs(tag = "d"), ncol = 1,
                      heights = c(1, 0.26))
# no "& theme_pgd()" over the whole: it would bring back the matrix's blanked strips
p_dis <- wrap_plots(p_matrix + labs(tag = "a"), p_right, ncol = 2, widths = c(0.4, 0.6)) &
  theme(text = element_text(family = FIG_FAMILY),
        plot.tag = element_text(size = 8, face = "bold"))

if (!exists("PANEL_EXPORT_ONLY")) {
  save_fig(p_dis, file.path(OUT_DIR, "figures"), "05_protein_disease", 168)
  # Slides: drug support by fluid-support class.
  # Supplementary figure: the disease categories (a main-figure panel until 2026-10-25),
  # one panel since the unreported-link forest was dropped (2026-11-13).
  save_fig(p_cat + theme_pgd(), file.path(OUT_DIR, "figures"), "05_supp_disease", 95)
  message("Categories: ", uniqueN(d$disease_category),
          "; joint proteins CSF ", d[fluid == "CSF", sum(n)],
          ", plasma ", d[fluid == "Plasma", sum(n)])
}
