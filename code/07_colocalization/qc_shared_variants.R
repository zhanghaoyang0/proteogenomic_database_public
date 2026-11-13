#!/usr/bin/env Rscript
# How well the shared variants cover the two peaks a coloc result rests on.
#
# coloc.susie aligns two saved SuSiE fits by rsid and uses their intersection.
# The two fits are centred on their own clump leads, which may sit up to 500 kb
# apart, so the shared set can thin out exactly around one of the peaks. The
# runner's only guard is MIN_N_OVERLAP = 50 shared variants, which does not
# notice this: across 600 reported protein x disease pairs the outcome peak had
# a median of 0.71 of its neighbouring variants shared (p5 = 0.40), and coverage
# was uncorrelated with PP4 -- a sparse pair looks exactly like a dense one.
#
# For each reported pair this records the fraction of each fit's variants within
# +/-50 kb (and +/-200 kb) of its own reported hit that are in the shared set,
# so a pair resting on thin coverage can be flagged in the supplement or dropped
# in a sensitivity analysis. Reads saved fits only: no SuSiE and no coloc rerun.
#
# Usage: Rscript qc_shared_variants.R [out.tsv.gz] [type ...]
suppressMessages({library(data.table); library(RSQLite)})
SCRIPT_DIR <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
source(file.path(SCRIPT_DIR, "functions_coloc.R"))

a     <- commandArgs(TRUE)
OUT   <- if (length(a) >= 1 && nzchar(a[1])) a[1] else
           file.path(COLOC_ROOT, "info", "shared_variant_qc.tsv.gz")
TYPES <- if (length(a) >= 2) a[-1] else
           c("protein_x_mb","protein_x_gtex","protein_x_sb","protein_x_disease","mb_x_disease")
WIN   <- c(`50kb` = 5e4, `200kb` = 2e5)

# Keyed join, not a named vector: `POS[snps]` on 9.9M names rebuilds a hash on
# every call, which turned a minutes-long pass into an hours-long one.
BIM <- fread(paste0(LD_REF, ".bim"),
             select = c(2, 4), col.names = c("SNP", "BP"))
setkey(BIM, SNP)
pos_of <- function(snps) BIM[.(snps), BP]

# Saved fits are read once and reused: rows are processed exposure-major, so a
# protein locus paired with forty diseases is loaded once, not forty times.
cache <- new.env(parent = emptyenv()); cached <- character()
getfit <- function(p) {
  k <- basename(p)
  if (!is.null(cache[[k]])) return(cache[[k]])
  # Stop, not NULL: a reported row whose fit is missing or unreadable used to vanish
  # from the QC table.
  if (!file.exists(p)) stop("missing fit: ", p)
  f <- rekey_susie(readRDS(p))
  # Positions once per fit, not once per pair: a protein locus is reused across
  # every disease it was tested against.
  if (!is.null(f)) f$pos <- pos_of(f$snp)
  if (length(cached) >= 400L) { rm(list = cached[1L], envir = cache); cached <<- cached[-1L] }
  assign(k, f, cache); cached <<- c(cached, k); f
}
# exposure and outcome sides differ per chain; everything else is shared.
SIDES <- list(
  protein_x_mb      = list(e = function(r) protein_susie_path(protein_src(r$exposure_id), r$exposure_id, r$lead_snp),
                           o = function(r) metabrain_susie_path(r$tissue, r$ENSG, lead_or_na(r$o_lead))),
  protein_x_gtex    = list(e = function(r) protein_susie_path(protein_src(r$exposure_id), r$exposure_id, r$lead_snp),
                           o = function(r) eqtl_susie_path("gtex",        r$tissue, r$ENSG, lead_or_na(r$o_lead))),
  protein_x_sb      = list(e = function(r) protein_susie_path(protein_src(r$exposure_id), r$exposure_id, r$lead_snp),
                           o = function(r) eqtl_susie_path("singlebrain", r$tissue, r$ENSG, lead_or_na(r$o_lead))),
  protein_x_disease = list(e = function(r) protein_susie_path(protein_src(r$exposure_id), r$exposure_id, r$lead_snp),
                           o = function(r) disease_susie_path(r$outcome_source, r$outcome_id, lead_or_na(r$o_lead))),
  mb_x_disease      = list(e = function(r) metabrain_susie_path(r$tissue, r$ENSG, r$lead_snp),
                           o = function(r) disease_susie_path(r$outcome_source, r$outcome_id, lead_or_na(r$o_lead)))
)

# Fraction of a fit's variants within `w` of its own hit that survive into the
# shared set. NA when the hit has no reference position or no neighbours.
covered <- function(f, hit, sh, w) {
  k <- match(hit, f$snp)
  if (is.na(k)) return(NA_real_)
  p0 <- f$pos[k]
  if (is.na(p0)) return(NA_real_)
  near <- !is.na(f$pos) & abs(f$pos - p0) <= w
  if (!any(near)) return(NA_real_)
  mean(f$snp[near] %in% sh)
}

res <- list()
for (ty in TYPES) {
  # The reported rows of the current result (coloc_reported), not every ok ledger row:
  # the ledger also holds BF1, retired panels and non-primary components.
  rf <- file.path(COLOC_ROOT, ty, "result", paste0("coloc_", ty, ".tsv"))
  if (!file.exists(rf)) stop("missing result: ", rf)
  d <- fread(rf)[coloc_reported(status, primary, PP4)]
  if (!nrow(d)) next
  if (ty %in% c("protein_x_disease", "mb_x_disease")) {
    osrc <- unique(fread(file.path(COLOC_ROOT, ty, "info", "coloc_pair.csv"),
                         select = c("outcome_id", "outcome_source")))
    if (anyDuplicated(osrc$outcome_id)) stop("an outcome_id with two outcome_source values")
    d[, outcome_source := osrc$outcome_source[match(outcome_id, osrc$outcome_id)]]
    if (anyNA(d$outcome_source)) stop(ty, ": outcomes missing from coloc_pair.csv")
  }
  keys <- intersect(c("exposure_id","tissue","ENSG","outcome_id","outcome_source","lead_snp","o_lead","e_hit","o_hit","PP4","n_snps"),
                    names(d))
  d <- unique(d[, ..keys])
  setorderv(d, intersect(c("exposure_id","lead_snp"), names(d)))   # exposure-major: cache hits
  side <- SIDES[[ty]]
  t0 <- Sys.time()
  out <- rbindlist(lapply(seq_len(nrow(d)), function(i) {
    r <- d[i]
    fe <- getfit(side$e(r)); fo <- getfit(side$o(r))
    if (is.null(fe) || is.null(fo)) return(NULL)
    sh <- intersect(fe$snp, fo$snp)
    cbind(r, data.table(
      n_shared            = length(sh),
      frac_exposure_kept  = round(length(sh) / length(fe$snp), 3),
      frac_outcome_kept   = round(length(sh) / length(fo$snp), 3),
      shared_frac_e_50kb  = round(covered(fe, r$e_hit, sh, WIN[["50kb"]]),  3),
      shared_frac_o_50kb  = round(covered(fo, r$o_hit, sh, WIN[["50kb"]]),  3),
      shared_frac_e_200kb = round(covered(fe, r$e_hit, sh, WIN[["200kb"]]), 3),
      shared_frac_o_200kb = round(covered(fo, r$o_hit, sh, WIN[["200kb"]]), 3)))
  }), fill = TRUE)
  if (!nrow(out)) next
  out[, coloc_type := ty]
  out[, low_shared_coverage := (shared_frac_e_50kb < 0.7 | shared_frac_o_50kb < 0.7)]
  res[[ty]] <- out
  cat(sprintf("%-18s %6d reported pairs, %6d evaluated  (%.1f min)  flagged %.1f%%\n",
              ty, nrow(d), nrow(out), as.numeric(difftime(Sys.time(), t0, units = "mins")),
              100 * mean(out$low_shared_coverage, na.rm = TRUE)))
}
all <- rbindlist(res, fill = TRUE)
setcolorder(all, c("coloc_type", intersect(names(all), c("exposure_id","tissue","ENSG","outcome_id",
                   "lead_snp","o_lead","e_hit","o_hit","PP4","n_snps"))))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(all, OUT, sep = "\t", compress = "gzip")
cat(sprintf("\n%d pairs | flagged (either peak < 0.7 within 50 kb): %d (%.1f%%)\nwrote -> %s\n",
            nrow(all), sum(all$low_shared_coverage, na.rm = TRUE),
            100 * mean(all$low_shared_coverage, na.rm = TRUE), OUT))
