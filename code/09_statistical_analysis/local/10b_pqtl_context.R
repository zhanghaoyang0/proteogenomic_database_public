#!/usr/bin/env Rscript
# LOCAL (Mac). pQTL effects across the AD continuum (docs/pipeline.md, 10b).
# For every distinct signal outside the MHC, in its discovery fluid:
#   protein ~ G + M + G x M + age + sex + pc1-10 + mean_protein
# with M = Abeta status (CSF Abeta42/40 abnormal), and the SNP effect within Abeta- and
# Abeta+ separately. (An APOE e4 dose modifier was dropped on 2026-10-29: not reported,
# author decision 2026-10-27.)
# APOE-region signals are fitted and flagged; the summaries leave them out, because e4
# drives both Abeta and many proteins.
#
# Single interactions are underpowered at this N, so the summary also reports the
# atlas-level agreement of Abeta- and Abeta+ effects and the inflation of the
# interaction P values.
#
# Writes summary statistics only (PILOT=1: 20 signals, to tmp/proteogenomic_database/ instead):
#   output/10b_pqtl_context.csv   one row per signal x modifier
#   output/10b_summary.csv        by fluid x cis/trans x modifier
#
# Usage (repo root, Mac; plink2 on PATH):
#   PILOT=1 Rscript code/09_statistical_analysis/local/10b_pqtl_context.R
#           Rscript code/09_statistical_analysis/local/10b_pqtl_context.R

suppressPackageStartupMessages(library(data.table))
for (pkg in c("readxl"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("R package '", pkg, "' is not installed")
source("code/09_statistical_analysis/lib/local_data.R")

PILOT   <- Sys.getenv("PILOT") == "1"
OUT_DIR <- if (PILOT) "tmp/proteogenomic_database" else "output"
PCS     <- paste0("pc", 1:10)
FDR_MAX <- 0.05
N_MIN   <- 100     # per Abeta stratum; below this the signal is not tested
R_MIN   <- 0.95    # GWAS-beta reproduction the run must reach
MODIFIERS <- c(abeta = "amyloid_status")
dir.create(OUT_DIR, showWarnings = FALSE)

# ---- signals and data -------------------------------------------------------
sig <- atlas_signals()
if (PILOT) sig <- rbind(head(sig[fluid == "CSF"], 10), head(sig[fluid == "Plasma"], 10))
message(sprintf("signals: %d (%d CSF, %d plasma; %d in the APOE region)%s", nrow(sig),
                sig[fluid == "CSF", .N], sig[fluid == "Plasma", .N], sig[apoe_region == TRUE, .N],
                if (PILOT) " [PILOT]" else ""))

sigdose <- signal_dosages(sig, if (PILOT) "10b_dosage_pilot" else "10b_dosage", max_missing = 0.005)
sig <- sigdose$sig; dose <- sigdose$dose

ph  <- load_pheno()
cov <- load_covar()
miss <- setdiff(PCS, names(cov))
if (length(miss)) stop(COVAR, " lacks ", paste(miss, collapse = ", "))
pv  <- load_proteins(unique(sig[, .(protein_assay, panel_tag)]), cov)
dat <- merge(cov, pv, by = "sid")
n0 <- nrow(dat)
dat <- merge(dat, ph$base[, .(sid, amyloid_status, apoe_e4_dose)], by = "sid", all.x = TRUE)
dat <- merge(dat, dose, by = "IID", all.x = TRUE)
if (nrow(dat) != n0) stop("merges changed the row count")
message(sprintf("analysis table: %d participants; Abeta status known for %d (%d positive); e4 dose for %d",
                nrow(dat), sum(!is.na(dat$amyloid_status)), sum(dat$amyloid_status %in% 1),
                sum(!is.na(dat$apoe_e4_dose))))

# ---- models -----------------------------------------------------------------
coef_of <- function(fit, term) {
  cf <- summary(fit)$coefficients
  if (!term %in% rownames(cf)) stop("term ", term, " dropped from the model")
  cf[term, c(1, 2, 4)]
}

fit_one <- function(r) {
  g <- dat[[r$snp]]
  if (r$flip) g <- 2 - g
  y <- dat[[r$protein_assay]]
  X <- cov_matrix(dat, unique(c("age", "sex", PCS, protein_covars(r$panel_tag, cov))))
  ok0 <- complete.cases(g, y, X)
  b_full <- coef_of(lm(y[ok0] ~ g[ok0] + X[ok0, , drop = FALSE]), "g[ok0]")[1]

  rbindlist(lapply(names(MODIFIERS), function(m) {
    M <- dat[[MODIFIERS[[m]]]]
    ok <- ok0 & !is.na(M)
    out <- data.table(modifier = m, n = sum(ok), beta_full = b_full)
    if (m == "abeta" && min(sum(ok & M == 0), sum(ok & M == 1)) < N_MIN) return(out)
    gi <- g[ok]; Mi <- M[ok]; Xi <- X[ok, , drop = FALSE]; yi <- y[ok]
    int <- coef_of(lm(yi ~ gi * Mi + Xi), "gi:Mi")
    out[, `:=`(beta_int = int[1], se_int = int[2], p_int = int[3])]
    if (m == "abeta") {
      neg <- coef_of(lm(yi[Mi == 0] ~ gi[Mi == 0] + Xi[Mi == 0, , drop = FALSE]), "gi[Mi == 0]")
      pos <- coef_of(lm(yi[Mi == 1] ~ gi[Mi == 1] + Xi[Mi == 1, , drop = FALSE]), "gi[Mi == 1]")
      out[, `:=`(n_neg = sum(Mi == 0), n_pos = sum(Mi == 1),
                 beta_neg = neg[1], se_neg = neg[2], beta_pos = pos[1], se_pos = pos[2])]
    }
    out
  }), fill = TRUE)
}

t0 <- Sys.time()
res <- rbindlist(lapply(seq_len(nrow(sig)), function(i) {
  if (i %% 500 == 0) message(sprintf("  %d / %d signals", i, nrow(sig)))
  cbind(sig[i, .(cs_key, gene, protein_assay, platform, fluid, is_cis, apoe_region,
                 snp, chr, pos, ea, beta_gwas = beta)],
        fit_one(sig[i]))
}), fill = TRUE)
message(sprintf("fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- checks -----------------------------------------------------------------
rep_d <- unique(res[, .(cs_key, beta_full, beta_gwas)])
r_rep <- cor(rep_d$beta_full, rep_d$beta_gwas, use = "complete.obs")
message(sprintf("GWAS-beta reproduction: r = %.3f, same sign %.1f%% (%d signals)", r_rep,
                100 * mean(sign(rep_d$beta_full) == sign(rep_d$beta_gwas), na.rm = TRUE), nrow(rep_d)))
if (!(r_rep > R_MIN)) stop("GWAS betas not reproduced (r = ", round(r_rep, 3), " <= ", R_MIN, ")")

# ---- multiplicity and summaries ---------------------------------------------
# BH within modifier x fluid, APOE-region signals excluded from the family.
res[, fdr_int := NA_real_]
res[apoe_region == FALSE & !is.na(p_int), fdr_int := p.adjust(p_int, "BH"), by = .(modifier, fluid)]

# Abeta- vs Abeta+ agreement, on effects aligned to the GWAS direction so that a
# shrinking effect in Abeta+ shows as a ratio below 1 rather than cancelling out.
res[, `:=`(beta_neg_al = sign(beta_gwas) * beta_neg, beta_pos_al = sign(beta_gwas) * beta_pos)]
lambda <- function(p) median(qchisq(p, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
summ <- res[apoe_region == FALSE & !is.na(p_int), .(
  n_signals   = .N,
  n_fdr       = sum(fdr_int < FDR_MAX, na.rm = TRUE),
  n_nominal   = sum(p_int < 0.05),
  lambda_int  = lambda(p_int),
  median_se_int = median(se_int),
  # smallest interaction detectable at 80% power, nominal and Bonferroni alpha
  mde_nominal = median(se_int) * (qnorm(0.975) + qnorm(0.8)),
  mde_bonf    = median(se_int) * (qnorm(1 - 0.025 / .N) + qnorm(0.8)),
  r_neg_pos   = if (modifier[1] == "abeta") cor(beta_neg, beta_pos) else NA_real_,
  median_neg_al = if (modifier[1] == "abeta") median(beta_neg_al) else NA_real_,
  median_pos_al = if (modifier[1] == "abeta") median(beta_pos_al) else NA_real_,
  p_wilcox_pos_vs_neg = if (modifier[1] == "abeta")
    wilcox.test(beta_pos_al, beta_neg_al, paired = TRUE)$p.value else NA_real_
), keyby = .(modifier, fluid, signal = fifelse(is_cis, "cis", "trans"))]
print(summ)

f_res  <- file.path(OUT_DIR, "10b_pqtl_context.csv")
f_summ <- file.path(OUT_DIR, "10b_summary.csv")
fwrite(res, f_res)
fwrite(summ, f_summ)
for (f in c(f_res, f_summ)) if (!file.exists(f)) stop("not written: ", f)
message("10b done: ", uniqueN(res$cs_key), " signals x ", length(MODIFIERS), " modifiers -> ", f_res)
