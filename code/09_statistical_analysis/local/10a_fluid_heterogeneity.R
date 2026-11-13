#!/usr/bin/env Rscript
# LOCAL (Mac). Paired CSF-plasma genetic effects (docs/pipeline.md, 10a).
# Every distinct signal of a protein assayed in both fluids is tested in CSF and in
# plasma in the same participants, and the difference is tested directly. Both fluids
# share one design matrix X = [G, age, sex, pc1-10, mean_protein of both panels], so
# lm(csf - plasma ~ X) gives beta_csf - beta_plasma exactly, and its SE carries the
# within-person correlation of the two measurements without a mixed model.
#
# Effects are on the INT scale of each fluid (the GWAS phenotype). A difference in
# standardized effect can also come from a difference in non-genetic variance, and
# CSF Explore 3072 vs plasma Explore HT can differ by assay; see docs/pipeline.md.
#
# Writes summary statistics only (PILOT=1: 20 signals, to tmp/proteogenomic_database/ instead):
#   output/10a_fluid_heterogeneity.csv   one row per signal
#
# Usage (repo root, Mac; plink2 on PATH):
#   PILOT=1 Rscript code/09_statistical_analysis/local/10a_fluid_heterogeneity.R
#           Rscript code/09_statistical_analysis/local/10a_fluid_heterogeneity.R

suppressPackageStartupMessages(library(data.table))
for (pkg in c("readxl"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("R package '", pkg, "' is not installed")
source("code/09_statistical_analysis/lib/local_data.R")

PILOT   <- Sys.getenv("PILOT") == "1"
OUT_DIR <- if (PILOT) "tmp/proteogenomic_database" else "output"
PCS     <- paste0("pc", 1:10)
FDR_MAX <- 0.05
N_MIN   <- 100     # paired participants below this: not tested
R_MIN   <- 0.95    # GWAS-beta reproduction the run must reach
dir.create(OUT_DIR, showWarnings = FALSE)

# ---- signals ----------------------------------------------------------------
kept <- fread("output/00_kept_assays.csv")
csf  <- kept[fluid == "CSF", .(ensembl_gene_id, csf_assay = protein_assay, csf_tag = tag,
                               csf_platform = platform)]
pl   <- kept[fluid == "Plasma", .(ensembl_gene_id, plasma_assay = protein_assay, plasma_tag = tag)]
if (anyDuplicated(csf$ensembl_gene_id) || anyDuplicated(pl$ensembl_gene_id))
  stop("00_kept_assays.csv is not one assay per protein per fluid")
pairs_assay <- merge(csf, pl, by = "ensembl_gene_id")
message(sprintf("proteins assayed in both fluids: %d", nrow(pairs_assay)))

sig <- atlas_signals()[ensembl_gene_id %in% pairs_assay$ensembl_gene_id]
# A lead SNP called in both fluids for the same protein is one test; it keeps the CSF
# row's effect allele and is labelled discovered in both.
sig[, n_fl := uniqueN(fluid), by = .(ensembl_gene_id, snp)]
setorder(sig, ensembl_gene_id, snp, fluid)
sig <- unique(sig, by = c("ensembl_gene_id", "snp"))
sig[, disc_fluid := fifelse(n_fl == 2, "both", fluid)]
n0 <- nrow(sig)
sig <- merge(sig, pairs_assay, by = "ensembl_gene_id")
if (nrow(sig) != n0) stop("assay pairing changed the signal count")
if (anyNA(sig$shared)) stop("signals of proteins assayed in both fluids without a sharing call")
if (PILOT) sig <- rbind(head(sig[fluid == "CSF"], 10), head(sig[fluid == "Plasma"], 10))
message(sprintf("signals tested: %d (%d CSF, %d plasma, %d in both)%s", nrow(sig),
                sig[disc_fluid == "CSF", .N], sig[disc_fluid == "Plasma", .N],
                sig[disc_fluid == "both", .N], if (PILOT) " [PILOT]" else ""))

# ---- data -------------------------------------------------------------------
# One counted allele per SNP; a SNP that leads two proteins with different effect
# alleles is flipped per signal.
sigdose <- signal_dosages(sig, if (PILOT) "10a_dosage_pilot" else "10a_dosage", max_missing = 0.005)
sig <- sigdose$sig; dose <- sigdose$dose

cov <- load_covar()
miss <- setdiff(PCS, names(cov))
if (length(miss)) stop(COVAR, " lacks ", paste(miss, collapse = ", "))
prot_info <- unique(rbind(sig[, .(protein_assay = csf_assay, panel_tag = csf_tag)],
                          sig[, .(protein_assay = plasma_assay, panel_tag = plasma_tag)]))
pv  <- load_proteins(prot_info, cov)
dat <- merge(cov, pv, by = "sid")
n0 <- nrow(dat)
dat <- merge(dat, dose, by = "IID", all.x = TRUE)
if (nrow(dat) != n0) stop("dosage merge changed the row count")
message(sprintf("analysis table: %d participants; %d without dosages", nrow(dat),
                sum(is.na(dat[[names(dose)[2]]]))))

# ---- models -----------------------------------------------------------------

g_row <- function(y, g, X) {
  cf <- summary(lm(y ~ g + X))$coefficients
  if (rownames(cf)[2] != "g") stop("SNP term dropped from the model")
  cf[2, c(1, 2, 4)]
}

fit_one <- function(r) {
  g <- dat[[r$snp]]
  if (r$flip) g <- 2 - g
  yc <- dat[[r$csf_assay]]
  yp <- dat[[r$plasma_assay]]
  covs <- unique(c("age", "sex", PCS, protein_covars(r$csf_tag, cov),
                   protein_covars(r$plasma_tag, cov)))
  X <- cov_matrix(dat, covs)
  # Reproduction check: the discovery fluid on its full GWAS sample, GWAS covariates.
  y_disc <- if (r$fluid == "CSF") yc else yp
  covs_disc <- unique(c("age", "sex", PCS, protein_covars(if (r$fluid == "CSF") r$csf_tag else r$plasma_tag, cov)))
  Xd <- cov_matrix(dat, covs_disc)
  okd <- complete.cases(g, y_disc, Xd)
  b_full <- g_row(y_disc[okd], g[okd], Xd[okd, , drop = FALSE])[1]

  ok <- complete.cases(g, yc, yp, X)
  out <- data.table(n_full_disc = sum(okd), beta_full_disc = b_full, n_paired = sum(ok))
  if (sum(ok) < N_MIN) return(out)
  g <- g[ok]; X <- X[ok, , drop = FALSE]
  c_ <- g_row(yc[ok], g, X)
  p_ <- g_row(yp[ok], g, X)
  d_ <- g_row(yc[ok] - yp[ok], g, X)
  out[, `:=`(eaf = mean(g) / 2,
             beta_csf = c_[1], se_csf = c_[2], p_csf = c_[3],
             beta_plasma = p_[1], se_plasma = p_[2], p_plasma = p_[3],
             delta = d_[1], se_delta = d_[2], p_het = d_[3],
             resid_r = cor(resid(lm(yc[ok] ~ X)), resid(lm(yp[ok] ~ X))))]
  if (abs(out$delta - (out$beta_csf - out$beta_plasma)) > 1e-8)
    stop("delta != beta_csf - beta_plasma for ", r$cs_key, ": design matrices differ")
  out
}

t0 <- Sys.time()
res <- rbindlist(lapply(seq_len(nrow(sig)), function(i) {
  if (i %% 200 == 0) message(sprintf("  %d / %d signals", i, nrow(sig)))
  cbind(sig[i, .(cs_key, gene, ensembl_gene_id, csf_assay, plasma_assay, csf_platform,
                 fluid, disc_fluid, is_cis, shared, apoe_region, snp, chr, pos, ea,
                 beta_gwas = beta)],
        fit_one(sig[i]))
}), fill = TRUE)
message(sprintf("fitted in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- checks -----------------------------------------------------------------
r_rep <- cor(res$beta_full_disc, res$beta_gwas, use = "complete.obs")
sign_rep <- mean(sign(res$beta_full_disc) == sign(res$beta_gwas), na.rm = TRUE)
message(sprintf("GWAS-beta reproduction: r = %.3f, same sign %.1f%% (%d signals)",
                r_rep, 100 * sign_rep, sum(!is.na(res$beta_full_disc))))
if (!(r_rep > R_MIN)) stop("GWAS betas not reproduced (r = ", round(r_rep, 3), " <= ", R_MIN,
                           "): samples, alleles or transform differ from the GWAS")
message(sprintf("signals with < %d paired participants (not tested): %d",
                N_MIN, res[n_paired < N_MIN, .N]))

# ---- classes ----------------------------------------------------------------
# BH over the signals outside APOE +- 1 Mb only, the set every consumer reads; inside it
# fdr_het stays NA, as fdr_int in 10b (code review 2026-10-30, R11: with the 29 APOE
# signals in the family, VWA1 plasma sat at 0.0501).
res[, fdr_het := NA_real_]
res[apoe_region == FALSE, fdr_het := p.adjust(p_het, "BH")]
# The effect in the fluid where the signal was not discovered, on the discovery
# fluid's direction.
res[, `:=`(
  beta_disc  = fifelse(fluid == "CSF", beta_csf, beta_plasma),
  beta_other = fifelse(fluid == "CSF", beta_plasma, beta_csf),
  se_other   = fifelse(fluid == "CSF", se_plasma, se_csf),
  p_other    = fifelse(fluid == "CSF", p_plasma, p_csf))]
res[, same_dir := sign(beta_other) == sign(beta_disc)]
res[, other_ci_has_0 := abs(beta_other) < qnorm(0.975) * se_other]
# %in% TRUE throughout: an NA comparison must not fall through to a class.
res[, class := fcase(
  apoe_region == TRUE, "APOE region",
  is.na(p_het), "not tested",
  disc_fluid == "both" | shared %in% TRUE, "shared",
  (fdr_het < FDR_MAX) %in% TRUE & (other_ci_has_0 | !same_dir) %in% TRUE, "effect-specific",
  (fdr_het < FDR_MAX) %in% TRUE & same_dir %in% TRUE & !other_ci_has_0 %in% TRUE, "attenuated",
  same_dir %in% TRUE & (p_other < 0.05) %in% TRUE & (fdr_het >= FDR_MAX) %in% TRUE, "detection-limited",
  default = "unresolved")]

# Shared signals can still differ in size between fluids; het_sig splits them.
res[, het_sig := (fdr_het < FDR_MAX) %in% TRUE]
res[apoe_region == TRUE, het_sig := NA]   # not tested, like fdr_het (not "not heterogeneous")
res[, signal := fifelse(is_cis, "cis", "trans")]
print(dcast(res[apoe_region == FALSE], disc_fluid + signal ~ class,
            value.var = "cs_key", fun.aggregate = length))

f_res  <- file.path(OUT_DIR, "10a_fluid_heterogeneity.csv")
fwrite(res, f_res)
if (!file.exists(f_res)) stop("not written: ", f_res)
message("10a done: ", nrow(res), " signals -> ", f_res)
