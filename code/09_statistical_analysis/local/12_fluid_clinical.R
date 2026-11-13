#!/usr/bin/env Rscript
# LOCAL (Mac). Paired fluids and the clinical relevance of each fluid
# (docs/pipeline.md, section 12). For every protein assayed in both fluids
# with a 10a signal outside APOE +- 1 Mb, in CU participants with both measurements:
#   joint model  y ~ CSF + plasma + covariates   -> d = beta_CSF - beta_plasma, z_pheno
#   separate     y ~ CSF + covariates, y ~ plasma + covariates (same participants)
# Outcomes (lib/local_outcomes.R, as 11a/11c): HV, ILV and the thickness composites at the
# first MRI (lm) and x time (lmer, random intercept and slope); mPACC v3 x time (lmer);
# time to AD-CI (Cox). Exposures are the GWAS-input INT levels. Covariates: age, sex,
# APOE e4 carrier, mean_protein of both panels; education for mPACC and Cox.
#
# Writes summary statistics only (no participant rows):
#   output/12_fluid_clinical.csv   one row per protein x outcome
# PILOT=1: the first 10 proteins, to tmp/proteogenomic_database/ instead.
#
# Usage (repo root, Mac):
#   N_CORES=8 Rscript code/09_statistical_analysis/local/12_fluid_clinical.R

suppressPackageStartupMessages({
  library(data.table)
  library(lmerTest)
  library(survival)
})
source("code/09_statistical_analysis/lib/local_data.R")
source("code/09_statistical_analysis/lib/local_outcomes.R")   # load_mri_cu, load_clinical

PILOT   <- Sys.getenv("PILOT") == "1"
OUT_DIR <- if (PILOT) "tmp/proteogenomic_database" else "output"
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
N_MIN   <- 100   # participants with both levels and the outcome; fewer: not fitted

# ---- proteins -------------------------------------------------------------------------
fh <- fread("output/10a_fluid_heterogeneity.csv")[apoe_region == FALSE & !is.na(p_het)]
kept <- fread("output/00_kept_assays.csv")[, .(protein_assay, tag)]
prot <- unique(fh[, .(ensembl_gene_id, gene, csf_assay, plasma_assay)])
if (anyDuplicated(prot$ensembl_gene_id)) stop("10a pairs a protein with more than one assay per fluid")
n0 <- nrow(prot)
prot <- merge(prot, kept[, .(csf_assay = protein_assay, csf_tag = tag)], by = "csf_assay")
prot <- merge(prot, kept[, .(plasma_assay = protein_assay, plasma_tag = tag)], by = "plasma_assay")
if (nrow(prot) != n0) stop("assays missing from 00_kept_assays.csv")
setorder(prot, gene)
if (PILOT) prot <- head(prot, 10)
message(sprintf("proteins: %d%s", nrow(prot), if (PILOT) " [PILOT]" else ""))

# ---- individual data ----------------------------------------------------------------------
cov <- load_covar()
pv <- load_proteins(unique(rbind(prot[, .(protein_assay = csf_assay, panel_tag = csf_tag)],
                                 prot[, .(protein_assay = plasma_assay, panel_tag = plasma_tag)])), cov)
mp_cols <- grep("^mean_protein_", names(cov), value = TRUE)
ind <- merge(cov[, c("sid", "age", "sex", mp_cols), with = FALSE], pv, by = "sid")
if (uniqueN(na.omit(ind$sex)) != 2) stop("sex does not have two levels")
ind[, sex := as.numeric(factor(sex)) - 1]

M <- load_mri_cu()
K <- load_clinical()
e4 <- K$base[, .(sid, e4_carrier, edu)]
ind <- merge(ind, e4, by = "sid", all.x = TRUE)
brain <- merge(M$mri, ind, by = "sid")
cog   <- merge(K$long[, .(sid, years, mpacc)], ind, by = "sid")
cog   <- cog[sid %in% K$base[cu == TRUE, sid]]
surv  <- merge(K$cu[!is.na(ad_ci), .(sid, t_years, ad_ci)], ind, by = "sid")   # NA: no baseline Abeta
OUTCOMES <- M$measures[measure %in% c("primary", "thickness composite"), col]
message(sprintf("CU: %d with MRI and proteins (%d scans); %d with mPACC; %d in the AD-CI analysis (%d events)",
                uniqueN(brain$sid), nrow(brain), uniqueN(cog$sid), nrow(surv), sum(surv$ad_ci)))

# ---- models ----------------------------------------------------------------------------
ctrl <- lmerControl(optimizer = "bobyqa")
# Coefficients for x (the exposure terms named in `terms`) and, with two, their contrast.
grab <- function(b, V, terms) {
  out <- data.table(term = terms, est = b[terms], se = sqrt(diag(V)[terms]))
  if (length(terms) == 2) {
    d <- b[terms[1]] - b[terms[2]]
    sd <- sqrt(V[terms[1], terms[1]] + V[terms[2], terms[2]] - 2 * V[terms[1], terms[2]])
    out <- rbind(out, data.table(term = "difference", est = d, se = sd))
  }
  out[, p := 2 * pnorm(-abs(est / se))]
  out
}
fit <- function(kind, dat, y, x, adj) {
  rhs <- paste(c(x, adj), collapse = " + ")
  m <- tryCatch(switch(kind,
    cross = lm(as.formula(paste(y, "~", rhs)), dat),
    long  = suppressMessages(lmer(as.formula(paste(y, "~ years * (", rhs, ") + (years | sid)")), dat, control = ctrl)),
    cox   = coxph(as.formula(paste("Surv(t_years, ad_ci) ~", rhs)), dat)),
    # A failed fit drops one row; say which, so the coverage check below is not the
    # only trace of it.
    error = function(e) { message("fit failed (", kind, " ", y, " ~ ", rhs, "): ", conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  b <- if (kind == "long") fixef(m) else coef(m)
  V <- as.matrix(vcov(m))
  terms <- if (kind == "long") paste0("years:", x) else x
  r <- grab(b, V, terms)
  r[, n := if (kind == "long") length(unique(m@frame$sid)) else if (kind == "cox") m$n else nobs(m)]
  r
}

fit_protein <- function(k) {
  pr <- prot[k]
  cc <- pr$csf_assay; pp <- pr$plasma_assay
  mp <- unique(intersect(paste0("mean_protein_", c(pr$csf_tag, pr$plasma_tag)), mp_cols))
  base_adj <- c("age", "sex", "e4_carrier", mp)
  rows <- list()
  run <- function(kind, dat, y, adj, label, analysis) {
    dat <- dat[!is.na(get(cc)) & !is.na(get(pp))]
    if (uniqueN(dat$sid) < N_MIN) return(invisible())
    setnames(dat, c(cc, pp), c("csf", "plasma"))
    for (spec in list(list(x = c("csf", "plasma"), model = "joint"),
                      list(x = "csf", model = "CSF only"), list(x = "plasma", model = "plasma only"))) {
      r <- fit(kind, dat, y, spec$x, adj)
      if (!is.null(r)) rows[[length(rows) + 1]] <<- cbind(
        data.table(gene = pr$gene, ensembl_gene_id = pr$ensembl_gene_id, analysis = analysis,
                   outcome = label, model = spec$model), r)
    }
  }
  first <- brain[years == 0]
  for (y in OUTCOMES) {
    run("cross", first[, c("sid", y, cc, pp, base_adj), with = FALSE], y, base_adj, y, "baseline")
    run("long", brain[, c("sid", "years", y, cc, pp, base_adj), with = FALSE], y, base_adj, y, "longitudinal")
  }
  run("long", cog[, c("sid", "years", "mpacc", cc, pp, base_adj, "edu"), with = FALSE], "mpacc",
      c(base_adj, "edu"), "mPACC", "longitudinal")
  run("cox", surv[, c("sid", "t_years", "ad_ci", cc, pp, base_adj, "edu"), with = FALSE], NA,
      c(base_adj, "edu"), "AD-CI", "progression")
  rbindlist(rows)
}

t0 <- Sys.time()
out <- parallel::mclapply(seq_len(nrow(prot)), function(k) tryCatch(fit_protein(k), error = function(e) e),
                          mc.cores = N_CORES)
# A worker killed by the OS (e.g. out of memory) returns NULL, not an error.
bad <- vapply(out, function(o) is.null(o) || inherits(o, "error") || inherits(o, "try-error"), logical(1))
if (any(bad)) stop(sum(bad), " proteins failed in a worker, e.g. ",
                   if (is.null(out[[which(bad)[1]]])) "worker died (NULL)" else conditionMessage(out[[which(bad)[1]]]))
res <- rbindlist(out)
message(sprintf("fitted %d proteins in %.1f min", nrow(prot), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- checks -----------------------------------------------------------------------------------
# Every model of every outcome, not just the joint difference: the single-fluid models
# feed 12b and 12c.
cover <- res[model != "joint" | term == "difference", .(N = uniqueN(gene)), by = .(analysis, outcome, model)]
# Against the full grid of fits fit_protein() asks for: an outcome that failed for every
# protein has no rows at all and would pass a check over `res` alone (code review
# 2026-10-30, R14).
want <- rbind(CJ(analysis = c("baseline", "longitudinal"), outcome = OUTCOMES),
              data.table(analysis = c("longitudinal", "progression"), outcome = c("mPACC", "AD-CI")))
want <- want[, .(model = c("joint", "CSF only", "plasma only")), by = .(analysis, outcome)]
cover <- merge(want, cover, by = c("analysis", "outcome", "model"), all = TRUE)
if (anyNA(cover$model) || nrow(cover) != nrow(want)) stop("fits outside the expected grid")
cover[is.na(N), N := 0L]
print(cover)
if (any(cover$N < 0.9 * nrow(prot))) stop("an outcome was fitted for < 90% of proteins")
print(res[model == "joint" & term == "difference",
          .(proteins = .N, pct_p05 = round(100 * mean(p < 0.05), 1)), by = analysis])

f <- file.path(OUT_DIR, "12_fluid_clinical.csv")
fwrite(res, f)
if (!file.exists(f)) stop("not written: ", f)
message("12 done: ", uniqueN(res$ensembl_gene_id), " proteins -> ", f)
