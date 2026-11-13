#!/usr/bin/env Rscript
# Canonical analyzed-assay allowlist (single source of truth for dedup).
#
# Rule: start from BF2 assays, drop SomaLogic, require a fluid label and GWAS
# non-missing N > 1000. Then, within each cohort+fluid, when one protein
# (ensembl_gene_id) is measured by more than one assay, keep only the assay with
# the largest N (ties broken by protein_assay for determinism). This collapses
# the same-protein repeats (mostly CSF Olink+NuLISA) so each protein counts once
# per fluid. All per-protein proteogenomic_database summaries/figures filter to this list.
#
# Usage: Rscript code/09_statistical_analysis/00_dedup_assays.R
# Output: output/00_kept_assays.csv

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/common.R")   # assay_platform, assay_fluid, CHR_LEN_HG19

ENSG_CSV <- "info/bf_protein/ensg/assay_ensg.csv"
OUT_CSV  <- "output/00_kept_assays.csv"

a <- unique(fread(ENSG_CSV)[cohort == "bf2"], by = c("cohort", "protein_assay"))
a[, platform := assay_platform(tag)]
a[, fluid := assay_fluid(tag)]
a[, N := as.double(N)]

analyzed <- a[!is.na(fluid) & platform != "SomaLogic" & N > 1000]
n_pre <- nrow(analyzed)

setorder(analyzed, cohort, fluid, ensembl_gene_id, -N, protein_assay)
kept <- analyzed[, .SD[1L], by = .(cohort, fluid, ensembl_gene_id)]

keep_cols <- c("cohort", "protein_assay", "tag", "platform", "fluid",
               "ensembl_gene_id", "gene", "N")
kept <- kept[, ..keep_cols]

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
fwrite(kept, OUT_CSV)

stopifnot(nrow(kept) > 0)
message("Analyzed assays before dedup: ", n_pre)
message("Kept assays after max-N dedup: ", nrow(kept),
        " (removed ", n_pre - nrow(kept), ")")
message("  CSF: ", kept[fluid == "CSF", .N],
        "  Plasma: ", kept[fluid == "Plasma", .N])
message("Saved: ", OUT_CSV)
