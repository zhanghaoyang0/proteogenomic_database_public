#!/usr/bin/env Rscript
# Section 12: do proteins whose CSF levels track Aβ and tau pathology differ more between
# fluids in their clinical associations? Uses the A/T-stage categories of Pichet Binette et
# al., Nat Neurosci 2024 (BioFINDER-2 CSF Olink Explore; Supplementary Tables 1 and 5):
#   Early  changed from A-T- to A+T- (Aβ stage);  Core  changed at every stage;
#   Late   changed from A+T- to A+T+ (tau stage), up or down.
# Background: the proteins of 12b_protein_summary.csv that the same study measured (its
# Supplementary Table 1), so a protein it could not have selected is not counted as "other".
# Tests (Wilcoxon, two-sided): each category and all categories together against the other
# measured proteins, on clinical effect heterogeneity P and on the CSF-minus-plasma asymmetry.
# Also Early against Late directly, and P recomputed within four outcome types (baseline
# structure, structural change, cognitive decline, progression) from the joint-model CSF-minus-
# plasma differences in 12_fluid_clinical.csv, to show which outcomes carry each group's P.
# Writes output/12f_pathology_proteins.csv (one row per protein), output/12f_tests.csv and
# output/12f_outcome_tests.csv (one row per outcome type x group).
# Usage (repo root): Rscript code/09_statistical_analysis/12f_pathology_proteins.R

suppressPackageStartupMessages(library(data.table))
OUT_DIR <- "output"
XLSX <- "ref/pichetbinette_nn_2024_supp/MOESM3.xlsx"   # 41593_2024_1737_MOESM3_ESM.xlsx
if (!file.exists(XLSX)) stop("missing ", XLSX)
rx <- function(sheet, skip) as.data.table(suppressMessages(
  readxl::read_excel(XLSX, sheet, skip = skip, col_names = FALSE)))

st1 <- rx("ST1_Fig2-5", 3)[!is.na(...1)]      # all analysed proteins (first column: gene)
st5 <- rx("ST5_Fig3-4", 4)[!is.na(...1)]      # the 128 differentially abundant proteins
# Matched on the Olink assay, not the gene symbol (2026-11-13): the study used the same
# BioFINDER-2 CSF Olink data, and its "Olink with Panel" ID is our CSF assay without the
# "_o<k>csf" suffix. Symbols missed renamed genes (our CD97, PVRL4, ATPIF1, SPG20 are its
# ADGRE5, NECTIN4, ATP5IF1, SPART).
dap <- unique(st5[, .(olink_id = ...2, category = ...3)])
stopifnot(uniqueN(dap$olink_id) == 128, !anyDuplicated(dap$olink_id),
          setequal(dap$category, c("Early", "Core", "Late - Upregulated", "Late - Downregulated")),
          all(dap$olink_id %in% st1$...2))
dap[, category := c(Early = "Early", Core = "Core", `Late - Upregulated` = "Late, up",
                    `Late - Downregulated` = "Late, down")[category]]

pr <- fread(file.path(OUT_DIR, "12b_protein_summary.csv"))
n0 <- nrow(pr)
csf_id <- unique(fread(file.path(OUT_DIR, "00_kept_assays.csv"))[fluid == "CSF",
  .(ensembl_gene_id, olink_id = sub("_o[0-9]+csf$", "", protein_assay))])
stopifnot(!anyDuplicated(csf_id$ensembl_gene_id))
pr <- merge(pr, csf_id, by = "ensembl_gene_id")
stopifnot(nrow(pr) == n0)
pr <- pr[olink_id %in% st1$...2]
pr <- merge(pr, dap, by = "olink_id", all.x = TRUE)
stopifnot(nrow(pr) == uniqueN(pr$ensembl_gene_id))
pr[is.na(category), category := "Other"]
message(sprintf("%d of %d clinical proteins measured by Pichet Binette et al.; %d in a category",
                nrow(pr), n0, pr[category != "Other", .N]))

oth <- pr[category == "Other"]
test_one <- function(lab, x) data.table(
  group = lab, n = nrow(x), n_other = nrow(oth),
  median_P = median(x$P), median_P_other = median(oth$P),
  p_P = wilcox.test(x$P, oth$P)$p.value,
  median_csf = median(x$mean_absz_csf), median_plasma = median(x$mean_absz_plasma),
  median_asym = median(x$asym), median_asym_other = median(oth$asym),
  p_asym = wilcox.test(x$asym, oth$asym)$p.value)
tests <- rbind(test_one("Any category", pr[category != "Other"]),
               rbindlist(lapply(c("Early", "Core", "Late, up", "Late, down"),
                                function(k) test_one(k, pr[category == k]))),
               test_one("Late (up or down)", pr[category %in% c("Late, up", "Late, down")]))
tests <- rbind(tests, data.table(group = "Other", n = nrow(oth),
  median_csf = median(oth$mean_absz_csf), median_plasma = median(oth$mean_absz_plasma)), fill = TRUE)
el <- wilcox.test(log10(pr[category %in% c("Late, up", "Late, down"), P]),
                  log10(pr[category == "Early", P]), conf.int = TRUE)
tests <- rbind(tests, data.table(group = "Late vs Early", n = pr[category %like% "^Late", .N],
  n_other = pr[category == "Early", .N], p_P = el$p.value, fold_P = 10^el$estimate,
  fold_lo = 10^el$conf.int[1], fold_hi = 10^el$conf.int[2]), fill = TRUE)
print(tests, digits = 3)

# ---- by outcome type --------------------------------------------------------------------------
r <- fread(file.path(OUT_DIR, "12_fluid_clinical.csv"))[model == "joint" & term == "difference"]
r <- r[ensembl_gene_id %in% pr$ensembl_gene_id][, z := est / se]
r[, outcome_type := fcase(analysis == "baseline", "Baseline structure",
                          analysis == "progression", "Progression",
                          outcome == "mPACC", "Cognitive decline",
                          default = "Structural change")]
stopifnot(r[, uniqueN(outcome_type)] == 4, r[, .N, by = .(ensembl_gene_id, outcome_type)][, all(N >= 1)])
po <- r[, .(P = mean(z^2)), by = .(ensembl_gene_id, outcome_type)]
po <- merge(po, pr[, .(ensembl_gene_id, group = fifelse(category %like% "^Late", "Late", category))],
            by = "ensembl_gene_id")
stopifnot(nrow(po) == 4 * nrow(pr))
ot <- po[, {
  o <- P[group == "Other"]
  .SD[, .(n = .N, median_P = median(P),
          p_vs_other = if (group[1] == "Other") NA_real_ else wilcox.test(P, o)$p.value), by = group]
}, by = outcome_type]
ot <- merge(ot, po[, .(p_late_vs_early = wilcox.test(P[group == "Late"], P[group == "Early"])$p.value),
                   by = outcome_type], by = "outcome_type")
print(ot, digits = 3)
fwrite(ot, file.path(OUT_DIR, "12f_outcome_tests.csv"))

fwrite(pr[, .(gene, ensembl_gene_id, category, P, mean_absz_csf, mean_absz_plasma, asym, class)],
       file.path(OUT_DIR, "12f_pathology_proteins.csv"))
fwrite(tests, file.path(OUT_DIR, "12f_tests.csv"))
stopifnot(file.exists(file.path(OUT_DIR, "12f_tests.csv")))
message("12f done")
