#!/usr/bin/env Rscript
# LOCAL (Mac). Variant -> cognitive decline and clinical progression in BF2, following
# Lu et al.'s rs5848 manuscript (Fig. 2; archive/lina_rs5848/docs/ref), BF2 only.
#
#   cognition, CU:          mPACC ~ time * (dose + age + sex + e4 + education) + PC1-5 + (time | sid)
#   by Abeta-PET:           mPACC ~ dose * centiloid * time + lower-order terms
#                               + time * (age + sex + e4 + education) + PC1-5 + (time | sid)
#                           in the AD-focused sample (CU, and Abeta+ MCI or dementia), as the
#                           manuscript did; trajectories at the median centiloid of each tertile
#   progression, CU:        Cox, time to MCI or dementia (AC-CI) and to Abeta+ MCI or dementia
#                           (AD-CI), dose + age + sex + e4 + education + PC1-5; AD-CI also with
#                           baseline centiloid added
# mPACC is mPACC_v3 (the manuscript's mPACC-TMT-B); time is years from the baseline visit.
#
# Writes summary statistics only (no participant rows):
#   output/11c_snp_clinical_assoc.csv   one row per variant x model (beta or HR, se, P, n, events)
#   output/11c_snp_clinical_traj.csv    predicted mPACC by genotype x centiloid tertile x year
#   output/11c_snp_clinical_surv.csv    adjusted event-free probability by genotype x year
#   output/11c_snp_clinical_risk.csv    number at risk in each Cox model's sample by whole year
#
# PROTEINS (comma-separated protein_assay names) replaces the variants by protein levels
# (GWAS-input INT values, panel mean_protein added as a covariate); mPACC
# trajectories and AD-CI survival are then drawn at the median of each tertile of the
# level (CU), and the centiloid-interaction models are skipped.
#
# Usage (repo root, Mac; plink2 on PATH). SNPS = comma-separated rsid:counted_allele.
#   SNPS=rs5848:T Rscript code/09_statistical_analysis/local/11c_snp_clinical.R

suppressPackageStartupMessages({
  library(data.table)
  library(lmerTest)
  library(survival)
})
source("code/09_statistical_analysis/lib/local_data.R")
source("code/09_statistical_analysis/lib/local_outcomes.R")   # load_clinical

OUT_DIR <- "output"
PCS <- paste0("pc", 1:5)
PROTEINS <- Sys.getenv("PROTEINS")
PROT_MODE <- nzchar(PROTEINS)

# ---- clinical data (lib/local_outcomes.R) ---------------------------------------------
K <- load_clinical()
base <- K$base; long <- K$long; cu <- K$cu

# ---- genotype and covariates ----------------------------------------------------------
cov <- load_covar()
ex <- load_exposures(cov, PROTEINS, Sys.getenv("SNPS", "rs5848:T"), "11c_dosage")
dose <- ex$dose; snps <- ex$snps; extra_cov <- ex$extra_cov
MP <- unique(unlist(extra_cov))
ind <- merge(cov[, c("sid", "IID", "sex", PCS, MP), with = FALSE], dose, by = "IID")
if (uniqueN(na.omit(ind$sex)) != 2) stop("sex does not have two levels")
ind[, sex := as.numeric(factor(sex)) - 1]
# PERMUTE_SEED: give every participant another participant's genotypes (one shuffle for
# all variants), so each test is under the null -- the calibration check of the models.
if (nzchar(Sys.getenv("PERMUTE_SEED"))) {
  set.seed(as.integer(Sys.getenv("PERMUTE_SEED")))
  idx <- sample(nrow(ind))
  ind[, (snps$snp) := lapply(.SD, function(v) v[idx]), .SDcols = snps$snp]
  message("genotypes permuted across participants (seed ", Sys.getenv("PERMUTE_SEED"), "): null run")
}
base <- merge(base, ind, by = "sid")
cu <- merge(cu, ind, by = "sid")
long <- merge(long, base[, c("sid", "age", "edu", "e4_carrier", "centiloid", "cu", "ad_focused", "sex", PCS, MP, snps$snp),
                         with = FALSE], by = "sid")
message(sprintf("mPACC: %d visits on %d participants (%d CU); AD-focused %d",
                nrow(long), uniqueN(long$sid), uniqueN(long[cu == TRUE, sid]), uniqueN(long[ad_focused == TRUE, sid])))

# ---- models ---------------------------------------------------------------------------
fx0 <- c("age", "sex", "e4_carrier", "edu")
ctrl <- lmerControl(optimizer = "bobyqa")
row_of <- function(cf, term, col_p) data.table(term = term, est = cf[term, 1], se = cf[term, 2], p = cf[term, col_p])

res <- list(); trajs <- list(); survs <- list(); risks <- list()
for (s in snps$snp) {
  fx <- c(fx0, extra_cov[[s]])
  L <- copy(long)[, g := get(s)][!is.na(g)]
  # display levels: genotype 0/1/2, or the median of each tertile of the protein (CU)
  levels <- if (PROT_MODE) {
    g0 <- unique(L[cu == TRUE, .(sid, g)])$g
    as.numeric(tapply(g0, cut(g0, quantile(g0, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = FALSE), median))
  } else 0:2
  # cognition in CU
  f1 <- as.formula(paste("mpacc ~ years * (g +", paste(fx, collapse = " + "), ") +", paste(PCS, collapse = " + "), "+ (years | sid)"))
  m1 <- lmer(f1, L[cu == TRUE], control = ctrl)
  cf <- summary(m1)$coefficients
  res[[length(res) + 1]] <- cbind(data.table(snp = s, model = "mPACC, CU", n = length(unique(m1@frame$sid))),
                                  rbind(row_of(cf, "g", 5), row_of(cf, "years:g", 5)))
  if (PROT_MODE) {
    # mPACC trajectory by protein tertile, from the CU model
    b <- fixef(m1); V <- as.matrix(vcov(m1)); mf <- m1@frame
    grid <- CJ(k = 0:2, years = seq(0, floor(quantile(mf$years, 0.95)), by = 0.5))
    grid[, g := levels[k + 1L]]
    grid[, (c(fx, PCS)) := lapply(mf[, c(fx, PCS)], mean)]
    X <- model.matrix(delete.response(terms(lme4::nobars(f1))), grid)[, names(b), drop = FALSE]
    trajs[[s]] <- grid[, .(snp = s, genotype = k, tertile = "CU", centiloid_median = NA_real_, level = g, years,
                           fit = drop(X %*% b), se = sqrt(rowSums((X %*% V) * X)), n_tertile = NA_integer_)]
  } else {
    # by Abeta-PET (continuous centiloid; tertiles for display)
    A <- L[ad_focused == TRUE & !is.na(centiloid)]
    f2 <- as.formula(paste("mpacc ~ g * centiloid * years + years * (", paste(fx, collapse = " + "), ") +",
                           paste(PCS, collapse = " + "), "+ (years | sid)"))
    m2 <- lmer(f2, A, control = ctrl)
    cf <- summary(m2)$coefficients
    res[[length(res) + 1]] <- cbind(data.table(snp = s, model = "mPACC x centiloid, AD-focused", n = length(unique(m2@frame$sid))),
                                    row_of(cf, "g:centiloid:years", 5))
    # Two checks on the three-way test, whose P values were inflated in the scan
    # (lambda 1.6 over 303 variants while the other outcomes were near 1):
    #  (i) the same model with centiloid rank-normalized, against its long right tail;
    #  (ii) two-stage: each participant's mPACC slope (>= 2 visits, OLS), then
    #       slope ~ g * centiloid + covariates with HC3 standard errors.
    ci <- unique(A[, .(sid, centiloid)], by = "sid")[, centiloid_int := int(centiloid)]   # one value per person
    A <- merge(A, ci[, .(sid, centiloid_int)], by = "sid")
    f2i <- as.formula(sub("g \\* centiloid \\* years", "g * centiloid_int * years", paste(deparse(f2, width.cutoff = 500), collapse = "")))
    m2i <- lmer(f2i, A, control = ctrl)
    cf <- summary(m2i)$coefficients
    res[[length(res) + 1]] <- cbind(data.table(snp = s, model = "mPACC x centiloid (INT), AD-focused", n = length(unique(m2i@frame$sid))),
                                    row_of(cf, "g:centiloid_int:years", 5))
    sl <- A[, if (.N >= 2 && diff(range(years)) > 0) .(slope = coef(lm(mpacc ~ years))[2]), by = sid]
    S <- merge(sl, unique(A[, c("sid", "g", "centiloid", fx, PCS), with = FALSE], by = "sid"), by = "sid")
    ms <- lm(reformulate(c("g * centiloid", fx, PCS), "slope"), S)
    ct <- lmtest::coeftest(ms, vcov = sandwich::vcovHC(ms, type = "HC3"))
    res[[length(res) + 1]] <- data.table(snp = s, model = "mPACC slope x centiloid, two-stage HC3", n = nrow(S),
                                         term = "g:centiloid", est = ct["g:centiloid", 1], se = ct["g:centiloid", 2],
                                         p = ct["g:centiloid", 4])
    b <- fixef(m2); V <- as.matrix(vcov(m2)); mf <- m2@frame
    first <- unique(mf[, c("sid", "centiloid")])
    tert <- cut(first$centiloid, quantile(first$centiloid, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("L", "M", "H"))
    med <- tapply(first$centiloid, tert, median)
    grid <- CJ(g = 0:2, tertile = names(med), years = seq(0, floor(quantile(mf$years, 0.95)), by = 0.5))
    grid[, centiloid := med[tertile]]
    grid[, (c(fx, PCS)) := lapply(mf[, c(fx, PCS)], mean)]
    X <- model.matrix(delete.response(terms(lme4::nobars(f2))), grid)[, names(b), drop = FALSE]
    grid[, `:=`(fit = drop(X %*% b), se = sqrt(rowSums((X %*% V) * X)))]
    trajs[[s]] <- grid[, .(snp = s, genotype = g, tertile, centiloid_median = centiloid, years, fit, se,
                           n_tertile = as.integer(table(tert)[tertile]))]
  }
  # progression in CU
  C <- copy(cu)[, g := get(s)][!is.na(g)]
  rhs <- paste(c("g", fx, PCS), collapse = " + ")
  for (o in c("ac_ci", "ad_ci")) {
    m <- coxph(as.formula(paste("Surv(t_years,", o, ") ~", rhs)), C)
    cf <- summary(m)$coefficients
    res[[length(res) + 1]] <- data.table(snp = s, model = paste0("Cox ", toupper(sub("_", "-", o)), ", CU"),
                                         n = m$n, events = m$nevent, term = "g", est = cf["g", "exp(coef)"],
                                         se = cf["g", "se(coef)"], p = cf["g", "Pr(>|z|)"])
    nd <- as.data.table(lapply(C[, c(fx, PCS), with = FALSE], mean, na.rm = TRUE))[rep(1, 3)][, g := levels]
    sf <- survfit(m, newdata = nd)
    tt <- seq(0, floor(max(C$t_years)), by = 0.25)
    ss <- summary(sf, times = tt, extend = TRUE)
    survs[[length(survs) + 1]] <- rbindlist(lapply(1:3, function(k) data.table(
      snp = s, outcome = toupper(sub("_", "-", o)), genotype = k - 1, years = ss$time,
      surv = ss$surv[, k], lower = ss$lower[, k], upper = ss$upper[, k])))
    # participants still followed at each whole year, in the model's own sample (m$y)
    rk <- summary(survfit(m$y ~ 1), times = 0:floor(max(m$y[, 1])), extend = TRUE)
    risks[[length(risks) + 1]] <- data.table(snp = s, outcome = toupper(sub("_", "-", o)),
                                             years = rk$time, n_risk = rk$n.risk)
  }
  if (!PROT_MODE) {   # mis-calibrated (docs/pipeline.md, 11); kept for the scan tables only
    m <- coxph(as.formula(paste("Surv(t_years, ad_ci) ~", rhs, "+ centiloid")), C)
    cf <- summary(m)$coefficients
    res[[length(res) + 1]] <- data.table(snp = s, model = "Cox AD-CI, CU, + centiloid", n = m$n, events = m$nevent,
                                         term = "g", est = cf["g", "exp(coef)"], se = cf["g", "se(coef)"], p = cf["g", "Pr(>|z|)"])
  }
}
res <- rbindlist(res, fill = TRUE)
print(res[, .(snp, model, term, est = signif(est, 3), se = signif(se, 3), p = signif(p, 3), n, events)])

OUT_NAME <- Sys.getenv("OUT_NAME", "11c_snp_clinical")   # a scan writes its own tables
f_a <- file.path(OUT_DIR, paste0(OUT_NAME, "_assoc.csv"))
f_t <- file.path(OUT_DIR, paste0(OUT_NAME, "_traj.csv"))
f_s <- file.path(OUT_DIR, paste0(OUT_NAME, "_surv.csv"))
f_r <- file.path(OUT_DIR, paste0(OUT_NAME, "_risk.csv"))
fwrite(res, f_a); fwrite(rbindlist(trajs), f_t); fwrite(rbindlist(survs), f_s); fwrite(rbindlist(risks), f_r)
for (f in c(f_a, f_t, f_s, f_r)) if (!file.exists(f)) stop("not written: ", f)
message("11c done: ", nrow(snps), " variant(s) -> ", f_a)
