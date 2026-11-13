#!/usr/bin/env Rscript
# Per-assay GWAS sample size and genomic inflation for the analysed assays (server).
# 00_kept_assays.csv carries N = non-missing protein measurements, genotyped or not,
# which exceeds the analysis sample (CSF Olink up to 1,826 against 1,771 participants
# with CSF and genotypes). The GWAS N and lambda_gc are already computed per assay by the
# cleaning step (one <cohort>__<assay>.qc.tsv per GWAS); this script only joins them to the
# kept set. Read from the QC files, not annotation/summary.tsv, which holds only assays with
# a clumped lead. gwas_n = n_max, the participants in that GWAS; n_median is per variant.
# Output: output/00b_gwas_qc.csv (protein_assay, fluid, platform, measured_n, gwas_n,
#   gwas_n_median, n_variants, lambda_gc), read by 06 (gwas_n in S1-S3) and the manuscript.
#   output/00b_panel_summary.csv (cohort, fluid, platform, panel, n_assays, gwas_n_min,
#   gwas_n_max): one row per panel, BF2 over the kept assays above, BF1 over every cleaned
#   BF1 protein GWAS (replication only, no kept-set filter); the supplementary panel table.
# Usage (repo root): Rscript code/09_statistical_analysis/00b_gwas_qc.R

suppressPackageStartupMessages(library(data.table))

OUT_DIR <- "output"
KEPT <- file.path(OUT_DIR, "00_kept_assays.csv")
source("code/09_statistical_analysis/lib/common.R")   # GWAS_BASE
QC_DIR <- file.path(GWAS_BASE, "bf2/clean/protein_qc")
if (!file.exists(KEPT)) stop("missing input: ", KEPT)
if (!dir.exists(QC_DIR)) stop("missing QC directory: ", QC_DIR)

kept <- fread(KEPT)
kept[, qc_file := file.path(QC_DIR, paste0(cohort, "__", protein_assay, ".qc.tsv"))]
miss <- kept[!file.exists(qc_file), .N]
if (miss > 0) stop(miss, " kept assays have no QC file in ", QC_DIR, " -- rerun the clean step with QC_DIR")
qc <- rbindlist(lapply(kept$qc_file, function(f)
  fread(f, select = c("n_variants", "n_median", "n_max", "lambda_gc"))[, qc_file := f]))
if (nrow(qc) != nrow(kept)) stop("expected one QC row per assay, got ", nrow(qc), " for ", nrow(kept))
out <- merge(kept[, .(qc_file, protein_assay, fluid, platform, measured_n = N)], qc, by = "qc_file")
if (anyNA(out$n_max) || anyNA(out$lambda_gc)) stop(out[is.na(n_max) | is.na(lambda_gc), .N], " assays with NA N or lambda_gc")
# A GWAS cannot use more samples than were measured.
if (out[n_max > measured_n, .N] > 0) stop("GWAS N above measured N for ", out[n_max > measured_n, .N], " assays")
out <- out[, .(protein_assay, fluid, platform, measured_n, gwas_n = n_max, gwas_n_median = n_median,
               n_variants, lambda_gc)]
setorder(out, fluid, platform, protein_assay)
fwrite(out, file.path(OUT_DIR, "00b_gwas_qc.csv"))
print(out[, .(n_assays = .N, gwas_n_min = min(gwas_n), gwas_n_max = max(gwas_n),
              lambda_median = round(median(lambda_gc), 3), lambda_sd = round(sd(lambda_gc), 3),
              lambda_max = round(max(lambda_gc), 3)), by = .(fluid, platform)])

# ---- per-panel summary, both cohorts ---------------------------------------------------
# panel = the part after the last "__" (gene names can hold "_" but not "__")
bf2p <- out[, .(cohort = "BioFINDER-2", fluid, platform, panel = sub("^.*__", "", protein_assay), gwas_n)]
BF1_QC <- file.path(GWAS_BASE, "bf1/clean/protein_qc")
f1 <- list.files(BF1_QC, "^bf1__.*\\.qc\\.tsv$", full.names = TRUE)
pq1 <- fread(file.path(BF1_QC, "bf1_pqtl_panel_qc.tsv"))   # the clean step's own panel counts
if (length(f1) != sum(pq1$clean_assays)) stop(length(f1), " BF1 QC files, panel summary says ", sum(pq1$clean_assays))
q1 <- rbindlist(lapply(f1, fread, select = c("panel", "n_max")))
if (anyNA(q1$n_max)) stop("BF1 QC files with NA n_max")
bf1p <- q1[, .(cohort = "BioFINDER-1", fluid = fifelse(grepl("plasma", panel, ignore.case = TRUE), "Plasma", "CSF"),
               platform = fifelse(grepl("nulisa", panel, ignore.case = TRUE), "NuLISA", "Olink"), panel, gwas_n = n_max)]
pan <- rbind(bf2p, bf1p)[, .(n_assays = .N, gwas_n_min = min(gwas_n), gwas_n_max = max(gwas_n)),
                         by = .(cohort, fluid, platform, panel)]
if (pan[cohort == "BioFINDER-2", sum(n_assays)] != nrow(out)) stop("BF2 panel counts do not add up to the kept assays")
setorder(pan, -cohort, fluid, platform, panel)
fwrite(pan, file.path(OUT_DIR, "00b_panel_summary.csv"))
print(pan)
message("00b done: ", nrow(out), " assays -> ", file.path(OUT_DIR, "00b_gwas_qc.csv"),
        "; ", nrow(pan), " panels -> ", file.path(OUT_DIR, "00b_panel_summary.csv"))
