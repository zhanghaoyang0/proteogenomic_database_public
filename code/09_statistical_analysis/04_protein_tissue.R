#!/usr/bin/env Rscript
# Summarize BF2 cis-pQTL x tissue eQTL colocalization across three resources:
#   MetaBrain (5 brain regions), GTEx (49 tissues; 13 brain + 36 peripheral) and
#   SingleBrain (7 single-nucleus brain cell types).
# Reframed as "tissue-of-origin": do CSF cis-pQTLs colocalize preferentially with
# brain eQTLs while plasma maps to peripheral tissues? GTEx (whole-body) carries the
# brain-vs-peripheral contrast; MetaBrain is the brain-specific deep resource.
# GTEx-brain and MetaBrain-brain OVERLAP (same regions) so they are kept as separate
# resource rows and never pooled. The intersection with the protein-disease evidence
# (the triangulation) needs 05a and therefore lives in 06_supplementary_tables.R.
# These coloc tests run on cis windows (gene +/-500kb), so they are cis-pQTL x eQTL
# events. SomaLogic is excluded (plasma-only, trans-inflated; 00_dedup_assays.R drops it too) so the
# CSF-vs-plasma comparison rests on a fair Olink (and CSF NuLISA) basis.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source("code/09_statistical_analysis/lib/common.R")   # assay_platform, assay_fluid, CHR_LEN_HG19
source("code/09_statistical_analysis/lib/eqtl_coloc_scope.R")

args <- commandArgs(trailingOnly = TRUE)
out_class <- if (length(args) >= 1) args[1] else
  "output/04_protein_eqtl_tissue_class.csv"
out_spec <- if (length(args) >= 2) args[2] else
  "output/04_protein_eqtl_specificity.csv"
out_tissue_evidence <-
  "output/04_protein_tissue_pair_evidence.csv.gz"
out_xtab <-
  "output/04_eqtl_support_by_specificity.csv"
out_dumbbell <-
  "output/04_tissue_dumbbell.csv"

# Result paths, scope constants and the cis+cognate filter come from
# lib/eqtl_coloc_scope.R, shared with 04b_protein_celltype.R so both analyses select
# the same assays and the same rows.
mb_file   <- COLOC_RESULT("mb")
gtex_file <- COLOC_RESULT("gtex")
# SingleBrain is a single-nucleus resource, so its "tissue" is a brain CELL TYPE, one
# level of biology below the tissues in the other two resources. It is drawn as its own
# labelled block (04_supp_coloc.png) rather than mixed into the brain rows, so the cell-type
# question -- which brain cell type carries the regulation -- stays separable from this
# script's estimand (tissue-of-origin: does CSF map to brain and plasma to periphery?).
# The cell-type-specificity analysis proper lives in 04b_protein_celltype.R.
sb_file   <- COLOC_RESULT("sb")

annot <- load_assay_annot()

cb <- rbindlist(list(
  read_coloc(mb_file, "MetaBrain"),
  read_coloc(gtex_file, "GTEx"),
  read_coloc(sb_file, "SingleBrain")
), use.names = TRUE)
# GTEx brain regions are prefixed "brain_"; everything else is peripheral.
# MetaBrain and SingleBrain are brain-only by construction.
cb[, tissue_class := fifelse(
  resource == "GTEx" & !grepl("^brain_", tissue), "Peripheral", "Brain"
)]

# gx_attempted_all must be taken BEFORE the cis+cognate filter (see note below), so
# the scope function is applied in two steps here.
cb <- merge(cb, annot, by = "exposure_id", all = FALSE)
# Denominator for the fluid-specificity crosstab further down, taken BEFORE the
# cis+cognate scope filter: every protein for which a GTEx coloc test completed at
# all. Inside the scope nearly every attempted test clears PP4 > 0.8 -- the same
# saturation documented above the figure block -- so a within-scope "attempted"
# denominator would be no denominator at all (203 attempted, 201 colocalizing).
# Taken before the filter it is set by eQTL availability, which is independent of
# whether this protein's cis-pQTL was found.
gx_attempted_all <- cb[resource == "GTEx" & status %chin% c("ok", "no_sig_coloc"),
                       unique(gene)]

cb <- apply_cis_cognate_scope(cb[, setdiff(names(cb), names(annot)[names(annot) != "exposure_id"]), with = FALSE], annot)

# In-scope export: one row per assay x resource x tissue x pQTL credible set, before
# any aggregation. 06_supplementary_tables.R turns this into the protein-tissue
# supplementary table; exporting from here rather than re-deriving it there keeps one
# definition of the cis+cognate scope and of COLOC_PP4.
fwrite(cb, out_tissue_evidence)
message("Saved: ", out_tissue_evidence, " (", nrow(cb), " in-scope coloc rows)")
cb_retained <- cb[coloc_retained == TRUE]
cb_sig <- cb_retained[PP4 > COLOC_PP4]
# Per tissue this is a lower bound: a no_sig_coloc protein carries one row, in the
# tissue of its best comparison (since 2026-10-24; before, no tissue at all), not one
# per tissue tested. n_attempted_genes is not read by the manuscript.
cb_attempted <- cb[
  status %chin% c("ok", "no_sig_coloc") & is_cis == TRUE & cognate == TRUE
]

# ---- Per fluid x resource x tissue: unique proteins colocalizing and tests run ----
summary_dt <- merge(
  cb_retained[, .(n_retained_pairs = .N, n_retained_genes = uniqueN(gene)),
        by = .(fluid, resource, tissue_class, tissue)],
  cb_sig[, .(n_strong_pairs = .N, n_proteins = uniqueN(gene)),
         by = .(fluid, resource, tissue_class, tissue)],
  by = c("fluid", "resource", "tissue_class", "tissue"), all = TRUE
)
summary_dt <- merge(
  cb_attempted[, .(
    n_attempted_pairs = .N, n_attempted_genes = uniqueN(gene)
  ), by = .(fluid, resource, tissue_class, tissue)],
  summary_dt,
  by = c("fluid", "resource", "tissue_class", "tissue"), all = TRUE
)
# Per-resource totals: unique proteins across that resource's tissues (cannot sum the
# per-tissue counts because one protein may colocalize in several tissues).
res_tot <- merge(
  cb_retained[, .(tissue_class = "All", tissue = "All",
            n_retained_pairs = .N, n_retained_genes = uniqueN(gene)),
        by = .(fluid, resource)],
  cb_sig[, .(tissue_class = "All", tissue = "All",
             n_strong_pairs = .N, n_proteins = uniqueN(gene)),
         by = .(fluid, resource)],
  by = c("fluid", "resource", "tissue_class", "tissue"), all = TRUE
)
res_tot <- merge(
  cb_attempted[, .(
    tissue_class = "All", tissue = "All",
    n_attempted_pairs = .N, n_attempted_genes = uniqueN(gene)
  ), by = .(fluid, resource)],
  res_tot,
  by = c("fluid", "resource", "tissue_class", "tissue"), all = TRUE
)
summary_dt <- rbindlist(list(summary_dt, res_tot), use.names = TRUE)
for (col in c(
  "n_attempted_pairs", "n_attempted_genes", "n_retained_pairs",
  "n_retained_genes", "n_strong_pairs", "n_proteins"
)) {
  set(summary_dt, which(is.na(summary_dt[[col]])), col, 0L)
}
setorder(summary_dt, fluid, resource, tissue_class, tissue)
dir.create(dirname(out_class), recursive = TRUE, showWarnings = FALSE)
# Printed, not written: the per-tissue counts are a check on the two rollups below
# (out_class, out_spec), which are what the report actually reads.
print(summary_dt)

# ---- Brain vs Peripheral rollup (unique colocalizing proteins per class) ----
class_dt <- cb_sig[, .(n_proteins = uniqueN(gene)),
                   by = .(fluid, resource, tissue_class)]
class_retained <- cb_retained[, .(n_retained_genes = uniqueN(gene)),
                      by = .(fluid, resource, tissue_class)]
class_dt <- merge(class_dt, class_retained,
                  by = c("fluid", "resource", "tissue_class"), all = TRUE)
class_attempted <- cb_attempted[, .(n_attempted_genes = uniqueN(gene)),
                                by = .(fluid, resource, tissue_class)]
class_dt <- merge(class_dt, class_attempted,
                  by = c("fluid", "resource", "tissue_class"), all = TRUE)
for (col in c("n_proteins", "n_retained_genes", "n_attempted_genes")) {
  set(class_dt, which(is.na(class_dt[[col]])), col, 0L)
}
setorder(class_dt, fluid, resource, tissue_class)
fwrite(class_dt, out_class)

# ---- eQTL support by fluid-specificity class -------------------------------------
# Does a protein's fluid-specificity class predict whether its cis-pQTL colocalizes
# with a BRAIN eQTL? GTEx only, so brain and peripheral come from one resource on one
# signal type. Two denominators, because they answer different questions and the
# manuscript reports both:
#   tested        every protein for which a GTEx coloc was attempted. A rate: how
#                 often does a protein of this class get brain support at all.
#   colocalizing  only proteins that colocalized somewhere in GTEx. A composition:
#                 what share of this class's expression support comes from brain.
# The tested denominator is legitimate here, unlike inside the coloc output used for
# the figure's head-to-head panel, because "attempted" is set by eQTL availability
# rather than by the pQTL already having been found.
spec_file <- "output/03_fluid_specificity_protein.csv"
if (!file.exists(spec_file)) {
  stop("Missing ", spec_file, " -- run 03_fluid_specificity.R first.")
}
spec_cls <- fread(spec_file)[class != "neither", .(gene, class)]
gx_sup <- cb_sig[resource == "GTEx", .(
  brain = any(tissue_class == "Brain"),
  peripheral = any(tissue_class == "Peripheral")
), by = gene]
xt <- merge(spec_cls, gx_sup, by = "gene", all.x = TRUE)
xt[is.na(brain), brain := FALSE][is.na(peripheral), peripheral := FALSE]
xt[, tested := gene %chin% gx_attempted_all]

xtab <- function(d, label) {
  out <- d[, .(
    n_proteins = .N,
    n_brain = sum(brain), pct_brain = round(100 * mean(brain), 1),
    n_peripheral = sum(peripheral), pct_peripheral = round(100 * mean(peripheral), 1)
  ), by = class][order(class)]
  out[, denominator := label]
  if (nrow(out) > 1) {
    out[, fisher_p := signif(fisher.test(as.matrix(d[, table(class, brain)]))$p.value, 3)]
  }
  out[]
}
xtab_dt <- rbindlist(list(
  xtab(xt[tested == TRUE], "attempted in GTEx"),
  xtab(xt[tested == TRUE & (brain | peripheral)], "colocalizing anywhere in GTEx")
), use.names = TRUE, fill = TRUE)
fwrite(xtab_dt, out_xtab)
message("eQTL support by fluid-specificity class:")
print(xtab_dt)

# ---- Tissue-of-origin specificity (GTEx only, brain vs peripheral per protein) ----
# Each colocalizing protein is classed Brain-only / Peripheral-only / Shared by which
# GTEx tissue classes carry a PP4 > 0.8 signal.
gx <- cb_sig[resource == "GTEx",
             .(classes = paste(sort(unique(tissue_class)), collapse = "+")),
             by = .(fluid, gene)]
gx[, spec := fcase(
  classes == "Brain", "Brain-only",
  classes == "Peripheral", "Peripheral-only",
  default = "Shared"
)]
spec_dt <- gx[, .(n_proteins = uniqueN(gene)), by = .(fluid, spec)]
setorder(spec_dt, fluid, spec)
fwrite(spec_dt, out_spec)


# ============================ Figure ============================
fluid_levels <- c("CSF", "Plasma")

# Dumbbell, matching 05_supp_disease.png a (05c_plot.R): one row per tissue, CSF and plasma joined
# by a grey segment.
#
# Two complementary analyses are drawn side by side for every tissue, following the
# design agreed for the manuscript: an inclusive discovery landscape and a restricted
# head-to-head comparison. The inclusive one shows the full scope of what was analysed
# and each fluid's discovery yield; the restricted one is what supports the
# fluid-specific claim, because a pattern that survives there is not driven by protein
# coverage or testability. See the block above `inclusive <-` for the exact
# definitions. The earlier single-panel versions of this figure -- raw counts alone,
# then a rate over each fluid's own testable set -- were both replaced because neither
# separates those two questions.
#
# Counts are printed in two fixed columns (CSF left, plasma right) so exact values are
# readable and the labels cannot collide in the rows where plasma sits left of CSF.

# Two complementary analyses per tissue, drawn as two sub-panels sharing the rows:
#
#   "All testable"     inclusive discovery landscape. Every protein for which a coloc
#                      test ran in that fluid at that tissue. Raw counts. This is
#                      descriptive -- it shows the scope of what was analysed and each
#                      fluid's discovery yield, and is NOT a fair CSF-vs-plasma
#                      comparison, because the two fluids test different protein sets.
#
#   "Shared testable"  restricted head-to-head. Only proteins for which a coloc test
#                      ran in BOTH fluids AT THAT TISSUE, so the denominator is the
#                      same protein set on both sides and eQTL discovery power is
#                      identical (same tissue). Plotted as a percentage of that shared
#                      set, with n/N printed, which is what supports the fluid-specific
#                      claim: a pattern that survives here is not driven by coverage or
#                      testability.
inclusive <- summary_dt[tissue != "All" & n_attempted_genes > 0,
                        .(fluid, resource, tissue_class, tissue, n = n_proteins)]

# Paired (head-to-head) set. The denominator MUST come from outside the coloc output:
# cb_attempted is the set of protein-tissue pairs for which a coloc test ran, and a
# test runs only where a cis-pQTL and an eQTL both exist, so nearly every attempted
# pair also clears PP4 > 0.8. Any denominator taken from inside the coloc results is
# therefore saturated -- at the resource level (187 of 192), and equally at the
# strictest per-tissue protein-pair level (62 of 63 at MetaBrain cortex). Pairing does
# not fix that; the denominator has to be defined before colocalization is attempted.
#
# So the paired denominator is the set of proteins measured in BOTH fluids that carry
# a cis-pQTL in BOTH fluids (cis_both in 03_fluid_specificity.R). It is one
# constant, identical for the two fluids and for every tissue, so each row asks: of
# the same N proteins, how many colocalize in CSF and how many in plasma at this
# tissue? Most of those N have no eQTL in any one tissue, so the rate is nowhere near
# saturated.
#
# Note what this restriction removes: proteins with a cis-pQTL in only one fluid. That
# is deliberate -- they cannot be compared head to head -- but it means this panel
# cannot carry the fluid-specificity claim on its own, because the fluid-specific
# proteins are exactly the ones it excludes. The CSF-only / shared / plasma-only
# gradient in the Results text is the analysis that addresses those.
if (!file.exists(spec_file)) {
  stop("Missing ", spec_file, " -- run 03_fluid_specificity.R first.")
}
shared_genes <- fread(spec_file)[cis_both == TRUE, unique(gene)]
n_shared <- length(shared_genes)
message("Paired denominator: ", n_shared,
        " proteins with a cis-pQTL in both fluids")

shared <- cb_sig[gene %in% shared_genes,
                 .(n = uniqueN(gene)),
                 by = .(fluid, resource, tissue_class, tissue)]
# Complete to both fluids so a tissue where one fluid has none still draws a zero.
shared <- merge(
  unique(inclusive[, .(resource, tissue_class, tissue)])[
    , .(fluid = fluid_levels), by = .(resource, tissue_class, tissue)],
  shared, by = c("fluid", "resource", "tissue_class", "tissue"), all.x = TRUE)
shared[is.na(n), n := 0L]
shared[, rate := 100 * n / n_shared]
# The per-tissue counts the figure draws, so the manuscript can say at how many
# tissues CSF exceeds plasma without rebuilding these two sets.
fwrite(rbind(inclusive[, .(analysis = "inclusive", fluid, resource, tissue_class, tissue, n)],
             shared[, .(analysis = "head_to_head", fluid, resource, tissue_class, tissue, n)]),
       out_dumbbell)

# The supplementary figures use these exported summaries (04c; the main-text panels are
# 10d). Run 04b first when refreshing SingleBrain breadth annotations.
source("code/09_statistical_analysis/04c_plot_tissue_context.R")
plot_tissue_context(dirname(out_dumbbell))

message(sprintf("GTEx tissue-of-origin (PP4 > %g) unique proteins:", COLOC_PP4))
print(dcast(spec_dt, fluid ~ spec, value.var = "n_proteins", fill = 0L))
message("Brain eQTL unique proteins by resource:")
print(dcast(cb_sig[tissue_class == "Brain",
                   .(n = uniqueN(gene)), by = .(fluid, resource)],
            fluid ~ resource, value.var = "n"))
message("Coloc-tested proteins per fluid x resource, and the shared (both-fluids) set:")
print(merge(dcast(cb_attempted[, .(n = uniqueN(gene)), by = .(resource, fluid)],
                  resource ~ fluid, value.var = "n"),
            cb_sig[gene %in% shared_genes, .(shared_pqtl_hits = uniqueN(gene)),
                   by = resource], by = "resource"))
message("Paired denominator (proteins with a cis-pQTL in both fluids): ", n_shared)
message("Saved: ", out_class)
message("Saved: ", out_spec)
message("Saved: ", out_xtab)
