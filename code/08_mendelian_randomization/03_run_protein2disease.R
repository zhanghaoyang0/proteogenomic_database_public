#!/usr/bin/env Rscript

# Per-chunk MR runner: rows [MR_PAIR_ROW_START, MR_PAIR_ROW_END] of mr_pair_todo.csv.gz.
# IVs: read from info/mr_iv.tsv.gz, written once per run by 02_extract_iv.R.
# Methods (n_iv >=2): IVW + Egger + Weighted median + Weighted mode + Q + Egger intercept.
# Methods (n_iv ==1): Wald ratio.
# Every pair: Steiger directionality (steiger_test(), functions_mr.R) on the same IVs.
# Storage: 1 row per pair -> tmp_chunks/chunk_{start}_{end}.tsv;
#          n_iv == 1: SNP info inlined in TSV (lead_snp/EA/NEA/EAF/beta_e/se_e/p_e/
#                     beta_o/se_o/p_o/N_o); no RDS written.
#          n_iv >= 2: nothing per SNP; the IVs are in mr_iv.tsv.gz and the outcome
#                     side in its clean file.

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

source(path.expand("code/lib/common.R"))   # GWAS_BASE, msg
source(path.expand("code/08_mendelian_randomization/functions_mr.R"))   # harmonize, steiger_test
MR_DIR      <- file.path(MR_ROOT, "protein2disease")
TMP_DIR     <- file.path(MR_DIR, "tmp_chunks")
PAIR_FILE   <- file.path(MR_DIR, "info", "mr_pair_todo.csv.gz")
IV_FILE     <- file.path(MR_DIR, "info", "mr_iv.tsv.gz")

dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

ROW_START <- as.integer(Sys.getenv("MR_PAIR_ROW_START", "1"))
ROW_END   <- as.integer(Sys.getenv("MR_PAIR_ROW_END",   "0"))
if (ROW_END < ROW_START) ROW_END <- ROW_START

t0  <- proc.time()[[3L]]

pair <- fread(PAIR_FILE)[ROW_START:ROW_END]
pair <- pair[!is.na(exposure_id) & !is.na(outcome_id)]
if (!"exposure_N" %in% names(pair) || anyNA(pair$exposure_N))
  stop("exposure_N missing in ", PAIR_FILE, " (Steiger needs it); rerun 02")
msg("rows ", ROW_START, "-", ROW_END, " -> ", nrow(pair), " pairs")

# The IV table must come from this todo list: an exposure it does not list was not
# extracted (an older table, or 02b not run), which is not the same as "no IV".
if (!file.exists(IV_FILE) || file.mtime(IV_FILE) < file.mtime(PAIR_FILE))
  stop(IV_FILE, " is missing or older than ", PAIR_FILE, " -- run 02_extract_iv.R")
iv_all <- fread(IV_FILE, colClasses = c(SNP = "character", EA_e = "character", NEA_e = "character"),
              na.strings = c("", "NA"))   # fwrite writes the no-IV row's NA SNP as ""
missing_expo <- setdiff(pair$exposure_id, iv_all$exposure_id)
if (length(missing_expo)) stop(length(missing_expo), " exposures not in ", IV_FILE,
                               ", e.g. ", missing_expo[1])
iv_all <- split(iv_all[!is.na(SNP)], by = "exposure_id", keep.by = FALSE)

outcome_path <- function(src, id) {
  sub <- if (src == "Other") "other/clean" else "finngen/clean"
  file.path(GWAS_BASE, sub, paste0(id, ".tsv.gz"))
}

oc_cache_id  <- ""
oc_cache_dat <- NULL
load_outcome <- function(id, src) {
  if (id == oc_cache_id) return(oc_cache_dat)
  p <- outcome_path(src, id)
  # Every outcome in the pair list was selected from an existing clean file, so a
  # missing one is a broken input, not a result. It used to become a "no_outcome"
  # row that the ledger then recorded as done. Fail the chunk instead.
  if (!file.exists(p)) stop("outcome file not found: ", p, call. = FALSE)
  d <- fread(p)
  oc_cache_id  <<- id
  oc_cache_dat <<- d
  d
}

run_mr <- function(hm) {
  dat <- data.frame(
    SNP = hm$SNP,
    beta.exposure = hm$BETA_e, se.exposure = hm$SE_e,
    beta.outcome  = hm$BETA_o, se.outcome  = hm$SE_o,
    effect_allele.exposure = hm$EA_e, other_allele.exposure = hm$NEA_e,
    effect_allele.outcome  = hm$EA_e, other_allele.outcome  = hm$NEA_e,
    pval.exposure = hm$P_e, pval.outcome = hm$P_o,
    id.exposure = "x", id.outcome = "y",
    exposure = "x", outcome = "y",
    mr_keep = TRUE
  )
  g <- function(x, f) if (!is.null(x) && nrow(x) >= 1L) x[[f]][1L] else NA_real_
  out <- list()
  if (nrow(dat) == 1L) {
    # Manual Wald ratio — skip TwoSampleMR overhead (data.frame + dispatch).
    b  <- hm$BETA_o[1L] / hm$BETA_e[1L]
    se <- abs(hm$SE_o[1L] / hm$BETA_e[1L])
    out$ivw_b  <- b
    out$ivw_se <- se
    out$ivw_p  <- if (is.finite(b) && is.finite(se) && se > 0) 2 * pnorm(-abs(b/se)) else NA_real_
    # Inline the single SNP, so the result table carries the whole n_iv == 1 analysis.
    out$lead_snp <- hm$SNP[1L]
    out$EA       <- hm$EA_e[1L]
    out$NEA      <- hm$NEA_e[1L]
    out$EAF      <- hm$EA_FREQ_o[1L]
    out$beta_e   <- hm$BETA_e[1L]; out$se_e <- hm$SE_e[1L]; out$p_e <- hm$P_e[1L]
    out$beta_o   <- hm$BETA_o[1L]; out$se_o <- hm$SE_o[1L]; out$p_o <- hm$P_o[1L]
    out$N_o      <- hm$N_o[1L]
  } else {
    # suppressMessages: mr() prints "Analysing ..." per call, four calls per pair (the
    # MR log dir had grown to 111 MB).
    mr <- function(...) suppressMessages(TwoSampleMR::mr(...))
    m_iv <- tryCatch(mr(dat, method_list = "mr_ivw"),               error = function(e) NULL)
    m_eg <- tryCatch(mr(dat, method_list = "mr_egger_regression"),  error = function(e) NULL)
    # Weighted median and mode bootstrap their SEs; a fixed seed per pair makes a rerun
    # reproduce wm/wmode se and p (unseeded up to the 2026-10-26 run, whose values
    # differ from the run before by a median 0.006 in wm_p).
    set.seed(1L)
    m_wm <- tryCatch(mr(dat, method_list = "mr_weighted_median"),   error = function(e) NULL)
    m_wo <- tryCatch(mr(dat, method_list = "mr_weighted_mode"),     error = function(e) NULL)
    hq   <- tryCatch(mr_heterogeneity(dat, method_list = "mr_ivw"), error = function(e) NULL)
    pl   <- tryCatch(mr_pleiotropy_test(dat),                       error = function(e) NULL)
    out$ivw_b      <- g(m_iv,"b"); out$ivw_se      <- g(m_iv,"se"); out$ivw_p      <- g(m_iv,"pval")
    out$egger_b    <- g(m_eg,"b"); out$egger_se    <- g(m_eg,"se"); out$egger_p    <- g(m_eg,"pval")
    out$wm_b       <- g(m_wm,"b"); out$wm_se       <- g(m_wm,"se"); out$wm_p       <- g(m_wm,"pval")
    out$wmode_b    <- g(m_wo,"b"); out$wmode_se    <- g(m_wo,"se"); out$wmode_p    <- g(m_wo,"pval")
    out$Q          <- g(hq, "Q");  out$Q_pval      <- g(hq, "Q_pval")
    out$egger_intercept    <- g(pl, "egger_intercept")
    out$egger_intercept_p  <- g(pl, "pval")
  }
  out
}

# "ok" only when the primary estimate exists. A failed or degenerate IVW fit (e.g. 2
# IVs with NA se) used to be recorded "ok", so the ledger marked it done and never
# retried it. Sensitivity methods may fail without invalidating the row.
validate_mr <- function(x) {
  if (!all(is.finite(c(x$ivw_b, x$ivw_se, x$ivw_p)))) return("invalid_ivw")
  "ok"
}

stub <- function(p, status, n_ivs = NA_integer_) {
  data.table(exposure_id = p$exposure_id, outcome_id = p$outcome_id,
             n_ivs = n_ivs, status = status)
}

results <- vector("list", nrow(pair))
for (i in seq_len(nrow(pair))) {
  p <- pair[i]
  iv <- iv_all[[p$exposure_id]]
  if (is.null(iv) || !nrow(iv)) { results[[i]] <- stub(p, "no_iv"); next }
  oc <- load_outcome(p$outcome_id, p$outcome_source)
  if (!any(iv$SNP %in% oc$SNP))   { results[[i]] <- stub(p, "no_overlap"); next }
  hm <- harmonize(iv, oc)
  if (is.null(hm) || !nrow(hm))   { results[[i]] <- stub(p, "harmo_fail"); next }
  mr_res <- run_mr(hm)
  mr_res$exposure_id <- p$exposure_id
  mr_res$outcome_id  <- p$outcome_id
  mr_res$n_ivs       <- nrow(hm)
  mr_res             <- c(mr_res, steiger_test(hm, p$exposure_N))
  mr_res$status      <- validate_mr(mr_res)
  results[[i]] <- as.data.table(mr_res)
  if (i %% 200L == 0L) msg("done ", i, "/", nrow(pair))
}

dt <- rbindlist(results, fill = TRUE)
setcolorder(dt, c("exposure_id", "outcome_id", "status", "n_ivs"))
chunk_file <- sprintf(file.path(TMP_DIR, "chunk_%07d_%07d.tsv"), ROW_START, ROW_END)
# Written under another name and renamed: 02 and 04 take any chunk_*.tsv as finished.
fwrite(dt, paste0(chunk_file, ".part"), sep = "\t")
if (!file.rename(paste0(chunk_file, ".part"), chunk_file)) stop("could not write ", chunk_file)
msg("wrote ", nrow(dt), " rows -> ", chunk_file)
