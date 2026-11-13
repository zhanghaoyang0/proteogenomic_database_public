#!/usr/bin/env Rscript
# LOCAL (Mac). Variant -> regional brain structure in cognitively unimpaired BF2
# participants, following Lu et al.'s rs5848 manuscript (Fig. 1; docs/ref of
# archive/lina_rs5848): rs5848 is run first to check this script against that figure,
# then the same code runs for the variants the atlas prioritizes.
#
#   cross-sectional, first MRI:  y ~ dose + age + sex + APOE e4 carrier + PC1-5      (lm)
#   longitudinal, all MRIs:      y ~ dose*years + age*years + sex*years + e4*years
#                                    + PC1-5 + (years | sid)                       (lmer)
# y: FreeSurfer (longitudinal stream, aseg_vol_* / aparc_*) regional volume / ICV per
# hemisphere, regional cortical thickness, and bilateral hippocampus and inferior
# lateral ventricle (HV, ILV; mean of L and R, / ICV) as in the manuscript.
# CU = baseline cognitive status Normal or SCD.
#
# Writes summary statistics only (no participant rows):
#   output/11a_snp_brain_assoc.csv   one row per variant x analysis x measure x region
#
# Longitudinal models are the slow part (about 1 s each). By default they are fitted only
# for the primary outcomes (HV, ILV and the thickness composites: the FDR family of the
# scan, docs/pipeline.md section 11); the regional ones get fit = "skipped".
# LONG_ALL=1 fits every region, for the variants drawn as brain maps (11b_plot_brain.R, removed 2026-10-28; in git history).
# OUT_NAME sets the output prefix (default 11a_snp_brain), so a scan does not overwrite
# the tables a figure is drawn from.
#
# PROTEINS (comma-separated protein_assay names, e.g. a CSF and a plasma assay of one
# protein) replaces the variants by protein levels: the GWAS-input INT value (mean 0,
# SD 1 already), with that panel's mean_protein as an extra covariate.
#
# Usage (repo root, Mac; plink2 on PATH). SNPS = comma-separated rsid:counted_allele.
#   SNPS=rs5848:T Rscript code/09_statistical_analysis/local/11a_snp_brain.R

suppressPackageStartupMessages({
  library(data.table)
  library(lmerTest)
})
for (pkg in c("readxl")) if (!requireNamespace(pkg, quietly = TRUE)) stop("R package '", pkg, "' is not installed")
source("code/09_statistical_analysis/lib/local_data.R")
source("code/09_statistical_analysis/lib/local_outcomes.R")   # load_mri_cu

OUT_DIR <- "output"
PCS <- paste0("pc", 1:5)
PROTEINS <- Sys.getenv("PROTEINS")
PROT_MODE <- nzchar(PROTEINS)
LONG_ALL <- Sys.getenv("LONG_ALL") == "1"
LONG_PRIMARY <- c("primary", "thickness composite")
OUT_NAME <- Sys.getenv("OUT_NAME", "11a_snp_brain")
MAX_FAIL <- 0.05   # share of longitudinal fits allowed to fail before the run stops

# ---- MRI (lib/local_outcomes.R) ----------------------------------------------------
M <- load_mri_cu()
mri <- M$mri; MEASURES <- M$measures; e4 <- M$e4

# ---- genotype and covariates ------------------------------------------------------
cov <- load_covar()
ex <- load_exposures(cov, PROTEINS, Sys.getenv("SNPS", "rs5848:T"), "11a_dosage")
dose <- ex$dose; snps <- ex$snps; extra_cov <- ex$extra_cov
miss <- setdiff(c("age", "sex", PCS), names(cov))
if (length(miss)) stop(COVAR, " lacks ", paste(miss, collapse = ", "))
MP <- unique(unlist(extra_cov))
ind <- merge(cov[, c("sid", "IID", "age", "sex", PCS, MP), with = FALSE], dose, by = "IID")
# sex as a 0/1 number, so every covariate can be set to its sample mean for predictions
if (uniqueN(na.omit(ind$sex)) != 2) stop("sex does not have two levels in ", COVAR)
ind[, sex := as.numeric(factor(sex)) - 1]
ind <- merge(ind, e4[, .(sid, e4_carrier)], by = "sid")
dat <- merge(mri, ind, by = "sid")
message(sprintf("CU BF2 with MRI and genotype: %d participants, %d scans (%.1f%% with >= 2 scans)",
                uniqueN(dat$sid), nrow(dat), 100 * mean(dat[, .N, by = sid]$N >= 2)))

# ---- models -------------------------------------------------------------------------
fx0 <- c("age", "sex", "e4_carrier", PCS)
fit_cross <- function(d, y) {
  f <- reformulate(c("g", fx), y)
  cf <- summary(lm(f, d))$coefficients
  sdy <- sd(d[[y]], na.rm = TRUE)
  data.table(beta = cf["g", 1], se = cf["g", 2], p = cf["g", 4], std_beta = cf["g", 1] / sdy,
             n = sum(complete.cases(d[, c(y, "g", fx), with = FALSE])))
}
fit_long <- function(d, y) {
  # time interactions for age, sex, e4 and (protein mode) mean_protein; PCs as main effects
  f <- as.formula(paste(y, "~", paste0(c("g", setdiff(fx, PCS)), " * years", collapse = " + "), "+",
                        paste(PCS, collapse = " + "), "+ (years | sid)"))
  m <- tryCatch(suppressMessages(lmer(f, d, control = lmerControl(optimizer = "bobyqa"))),
                error = function(e) NULL)
  if (is.null(m)) return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, std_beta = NA_real_,
                                    n = NA_integer_, fit = "failed"))
  cf <- summary(m)$coefficients
  sdy <- sd(d[years == 0][[y]], na.rm = TRUE)
  data.table(beta = cf["g:years", 1], se = cf["g:years", 2], p = cf["g:years", 5],
             std_beta = cf["g:years", 1] / sdy, n = length(unique(m@frame$sid)),
             fit = if (length(m@optinfo$conv$lme4$messages)) "warning" else "ok")
}

assoc <- list()
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
for (s in snps$snp) {
  t_snp <- Sys.time()
  fx <- c(fx0, extra_cov[[s]])
  d <- copy(dat)[, g := get(s)]
  d <- d[!is.na(g)]
  first <- d[years == 0]
  first[, grp := if (PROT_MODE) cut(g, quantile(g, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = FALSE) - 1L
                 else round(g)]
  message(sprintf("%s: %d CU at first MRI (%s 0/1/2 = %d/%d/%d)", s, nrow(first),
                  if (PROT_MODE) "tertile" else "genotype",
                  first[grp == 0, .N], first[grp == 1, .N], first[grp == 2, .N]))
  # One fork per outcome (the mixed models are the slow part); each returns its two rows.
  out <- parallel::mclapply(seq_len(nrow(MEASURES)), function(i) {
    y <- MEASURES$col[i]
    cr <- fit_cross(first, y)
    lo <- if (LONG_ALL || MEASURES$measure[i] %in% LONG_PRIMARY) fit_long(d, y) else
      data.table(beta = NA_real_, se = NA_real_, p = NA_real_, std_beta = NA_real_, n = NA_integer_, fit = "skipped")
    rbind(
      cbind(data.table(snp = s, analysis = "cross-sectional", measure = MEASURES$measure[i],
                       region = MEASURES$region[i], fit = "ok"), cr),
      cbind(data.table(snp = s, analysis = "longitudinal", measure = MEASURES$measure[i],
                       region = MEASURES$region[i]), lo), fill = TRUE)
  }, mc.cores = N_CORES)
  # A worker killed by the OS (e.g. out of memory) returns NULL, not a try-error.
  bad <- vapply(out, function(o) is.null(o) || inherits(o, "try-error"), logical(1))
  if (any(bad)) stop(sum(bad), " outcome fits crashed in a worker, e.g. ",
                     if (is.null(out[[which(bad)[1]]])) "worker died (NULL)" else out[[which(bad)[1]]])
  assoc <- c(assoc, out)
  message(sprintf("  %s: %d outcomes in %.1f min", s, nrow(MEASURES),
                  as.numeric(difftime(Sys.time(), t_snp, units = "mins"))))
}
res <- rbindlist(assoc, fill = TRUE)

# ---- checks -------------------------------------------------------------------------
n_fail <- res[analysis == "longitudinal" & fit == "failed", .N]
message(sprintf("longitudinal fits: %d ok, %d with a convergence warning, %d failed",
                res[analysis == "longitudinal" & fit == "ok", .N],
                res[analysis == "longitudinal" & fit == "warning", .N], n_fail))
if (n_fail > MAX_FAIL * res[analysis == "longitudinal" & fit != "skipped", .N]) stop("too many failed longitudinal fits")
print(res[measure == "primary", .(snp, analysis, region, beta = signif(beta, 3), se = signif(se, 3),
                                  p = signif(p, 3), n)])
message("nominal P < 0.05 per measure:")
print(res[, .(regions = .N, nominal = sum(p < 0.05, na.rm = TRUE)), keyby = .(snp, analysis, measure)])

f_a <- file.path(OUT_DIR, paste0(OUT_NAME, "_assoc.csv"))
fwrite(res, f_a)
if (!file.exists(f_a)) stop("not written: ", f_a)
message("11a done: ", nrow(snps), " variant(s) -> ", f_a)
