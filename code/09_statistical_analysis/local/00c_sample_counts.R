#!/usr/bin/env Rscript
# LOCAL (Mac). Participants per cohort and fluid in the pQTL analysis, for the Abstract and
# Methods "Study design and participants" (they used to be typed into the qmd).
# A participant counts for a fluid when they entered at least one GWAS of that fluid: in
# {cohort}_covar.phe with every covariate that GWAS used (age, sex, pc1-pc10 and the panel's
# mean_protein_<tag> when the file has it; PLINK drops a sample missing any) and a
# non-missing value for an analysed assay of the panel. Without the covariate condition
# (first run, 0b0493d) BF2 plasma was 1,617 against a largest plasma GWAS N of 1,526.
# Analysed assays: BF2, output/00_kept_assays.csv; BF1, every column of its Olink and
# NuLISA panels (the replication look-up set). SomaLogic panels are left out, as in the
# paper.
#
# Writes aggregate counts only (no participant rows):
#   output/00c_sample_counts.csv         cohort, fluid (CSF / Plasma / any / both), n
#   output/00c_sample_demographics.csv   per cohort, the "any" sample: n, baseline age
#                                        (min, max, mean, SD; years) and men (n, %).
#                                        Nothing reads it yet: kept for the population
#                                        characteristics of the reporting summary.
# sex is gender_baseline_variable as 02_process_covar.py copies it: 0 = man, 1 = woman
# (code/process_wgs/03_bf2_qc/03_bf2_apoe_report.py); any other value stops the script.
#
# Usage (repo root, Mac):
#   Rscript code/09_statistical_analysis/local/00c_sample_counts.R

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/local_data.R")   # PROT_DIR, norm_iid

OUT  <- "output/00c_sample_counts.csv"
OUT_DEMO <- "output/00c_sample_demographics.csv"
kept <- fread("output/00_kept_assays.csv")

tag_fluid <- function(tag) fifelse(grepl("csf", tag, ignore.case = TRUE), "CSF",
                          fifelse(grepl("plasma", tag, ignore.case = TRUE), "Plasma", NA_character_))

count_cohort <- function(cohort) {
  cov <- fread(file.path(PROT_DIR, paste0(cohort, "_covar.phe")), colClasses = list(character = "IID"))
  cov[, IID := norm_iid(IID)]
  pcs <- grep("^pc([1-9]|10)$", names(cov), ignore.case = TRUE, value = TRUE)
  if (!all(c("age", "sex") %in% names(cov)) || length(pcs) != 10)
    stop(cohort, "_covar.phe lacks age, sex or pc1-pc10: ", paste(names(cov), collapse = ", "))
  files <- list.files(PROT_DIR, pattern = paste0("^", cohort, "_.*[.]phe$"))
  files <- setdiff(files, paste0(cohort, "_covar.phe"))
  files <- files[!grepl("somalogic", files, ignore.case = TRUE)]
  if (!length(files)) stop("no panel .phe files for ", cohort, " in ", PROT_DIR)
  seen <- list(CSF = character(), Plasma = character())
  for (f in files) {
    tag <- sub(paste0("^", cohort, "_(.*)[.]phe$"), "\\1", f)
    fl <- tag_fluid(tag)
    if (is.na(fl)) stop("cannot tell the fluid of panel ", tag)
    hdr <- setdiff(names(fread(file.path(PROT_DIR, f), nrows = 0)), c("FID", "#FID", "IID"))
    # base indexing: inside kept[...], `tag` would be kept's own column, not this panel
    want <- if (cohort == "bf2") intersect(hdr, kept$protein_assay[kept$tag == tag]) else hdr
    if (cohort == "bf2" && !length(want)) next   # panel with no analysed assay
    ph <- fread(file.path(PROT_DIR, f), select = c("IID", want), colClasses = list(character = "IID"))
    ph[, IID := norm_iid(IID)]
    has <- ph[rowSums(!is.na(ph[, ..want])) > 0, IID]
    cc <- c("age", "sex", pcs, intersect(paste0("mean_protein_", tag), names(cov)))
    gwas_ids <- cov[complete.cases(cov[, ..cc]), IID]
    seen[[fl]] <- union(seen[[fl]], intersect(has, gwas_ids))
    message(sprintf("%s %s (%s): %d assays, %d with a value, %d of them with complete covariates (%s)",
                    cohort, tag, fl, length(want), length(has), length(intersect(has, gwas_ids)),
                    if (length(cc) > 12) cc[13] else "no mean_protein"))
  }
  if (cohort == "bf2") {   # every analysed BF2 assay must have been found in a .phe file
    miss <- setdiff(kept$protein_assay, unlist(lapply(files, function(f)
      names(fread(file.path(PROT_DIR, f), nrows = 0)))))
    if (length(miss)) stop(length(miss), " kept assays in no bf2 .phe file, e.g. ", miss[1])
  }
  any_ids <- union(seen$CSF, seen$Plasma)
  d <- cov[IID %in% any_ids]
  if (nrow(d) != length(any_ids)) stop(cohort, ": analysis sample not all in the covariate file")
  if (!all(d$sex %in% c(0, 1))) stop(cohort, ": sex values other than 0/1: ", paste(unique(d$sex), collapse = ", "))
  # built outside d[...]: inside it a bare `cohort` would be a column if the file had one
  demo <- data.table(cohort = cohort, n = nrow(d), age_min = floor(min(d$age)), age_max = floor(max(d$age)),
                     age_mean = round(mean(d$age), 1), age_sd = round(sd(d$age), 1),
                     n_men = sum(d$sex == 0), pct_men = round(100 * mean(d$sex == 0), 1))
  list(counts = data.table(cohort = cohort, fluid = c("CSF", "Plasma", "any", "both"),
                           n = c(length(seen$CSF), length(seen$Plasma),
                                 length(any_ids), length(intersect(seen$CSF, seen$Plasma)))),
       demo = demo)
}

out <- lapply(c("bf2", "bf1"), count_cohort)
res  <- rbindlist(lapply(out, `[[`, "counts"))
demo <- rbindlist(lapply(out, `[[`, "demo"))
print(res)
print(demo)
fwrite(res, OUT)
fwrite(demo, OUT_DEMO)
for (f in c(OUT, OUT_DEMO)) if (!file.exists(f)) stop("not written: ", f)
message("00c_sample_counts done -> ", OUT, ", ", OUT_DEMO)
