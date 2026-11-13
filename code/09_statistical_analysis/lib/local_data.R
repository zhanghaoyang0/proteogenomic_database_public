# Individual-level BioFINDER-2 data for the scripts in code/09_statistical_analysis/local/.
# Mac only: everything read here lives under ~/db/secure_data and never enters the
# repo. Only per-protein summary statistics are written back to output/.
#
# Sourced from the repo root, like every proteogenomic_database script:
#   source("code/09_statistical_analysis/lib/local_data.R")
# Requires data.table, loaded by the caller.

source("code/09_statistical_analysis/lib/common.R")      # in_apoe_region
source("code/09_statistical_analysis/lib/cs_signals.R")  # distinct_signals

SECURE <- path.expand("~/db/secure_data")
if (!dir.exists(SECURE)) stop("no ", SECURE, " -- code/09_statistical_analysis/local/ runs on the Mac, not the server")

PET_DIR    <- file.path(SECURE, "bf/pheno/raw/haoyang__20260227_112046_PET")
PET_X      <- file.path(PET_DIR, "haoyang__20260227_112046_x.csv")
TAU_QC     <- file.path(PET_DIR, "v2_tauPET_QC_for_database_June2023.xlsx")
SID_IID    <- file.path(SECURE, "bf/pheno/clean/bf2_sid_iid.csv")
PROT_DIR   <- file.path(SECURE, "bf/array/gwas_input/protein")   # the pQTL GWAS inputs
COVAR      <- file.path(PROT_DIR, "bf2_covar.phe")
PFILE      <- file.path(SECURE, "bf/array/bf2/imputed_qc/bf2_hg19")
WORK_DIR   <- file.path(SECURE, "projects/proteogenomic_database")     # dosages

for (f in c(PET_X, TAU_QC, SID_IID, COVAR, paste0(PFILE, ".pgen")))
  if (!file.exists(f)) stop("missing input: ", f)
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

norm_id <- function(x) trimws(sub("\\.0$", "", as.character(x)))
# IIDs lose their leading zeros in some files (0NNNNNNN vs NNNNNNN).
norm_iid <- function(x) {
  x <- sub("^0+", "", trimws(as.character(x)))
  ifelse(x == "", "0", x)
}
# Dates in the export come in an undocumented format: take the format that parses the
# most non-empty values and say which (only the value shapes are logged, e.g. 99999999).
parse_date <- function(x, what = "date") {
  x <- trimws(x); x[x == ""] <- NA
  fmts <- c("%Y-%m-%d", "%Y%m%d", "%d/%m/%Y", "%m/%d/%Y", "%Y/%m/%d", "%d.%m.%Y", "%d-%m-%Y", "%d%b%Y", "%d-%b-%Y")
  ok <- sapply(fmts, function(f) sum(!is.na(as.IDate(x, format = f))))
  if (grepl("^[0-9]{5}$", na.omit(x)[1])) {   # spreadsheet serial day numbers
    message(what, ": spreadsheet serial days")
    return(as.IDate(as.numeric(x), origin = "1899-12-30"))
  }
  # shape of the values only (digits -> 9, letters -> a), to diagnose a new format
  shapes <- head(sort(table(gsub("[A-Za-z]", "a", gsub("[0-9]", "9", na.omit(x)))), decreasing = TRUE), 3)
  message(what, " shapes: ", paste(names(shapes), shapes, sep = " x", collapse = "; "))
  best <- names(which.max(ok))
  message(sprintf("%s: format %s parses %d of %d non-empty values", what, best, max(ok), sum(!is.na(x))))
  as.IDate(x, format = best)
}

# Rank-based inverse normal, the same transform plink2 --pheno-quantile-normalize
# applied to every protein in the pQTL GWAS.
int <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - 0.5) / sum(!is.na(x)))
}

# Raw export column -> analysis name; the outcome definitions are in docs/pipeline.md (09).
PET_COLS <- c(
  sid = "sid", Study = "Study", Visit = "Visit", excluded = "excluded",
  amyloid_status = "Abnormal_CSF_Ab42_Ab40_Ratio",
  diagnosis      = "diagnosis_baseline_variable",
  cog_status     = "cognitive_status_baseline_variable",
  centiloid      = "CL_bfnc_ber_com_composite",
  tau_pet        = "tnic_cho_com_I_IV",
  thickness      = "ct_adsign_surfwt_lr",
  mmse           = "mmse_score",
  mpacc          = "mPACC_v1",
  apoe           = "apoe_genotype_baseline_variable"
)

# Returns list(base = one row per BF2 sid at Visit 0, long = sid x visit cognition).
load_pheno <- function() {
  hdr <- names(fread(PET_X, nrows = 0))
  miss <- setdiff(PET_COLS, hdr)
  if (length(miss)) stop("PET export is missing: ", paste(miss, collapse = ", "))
  pet <- fread(PET_X, select = unname(PET_COLS), colClasses = list(character = "sid"),
               showProgress = FALSE)
  setnames(pet, unname(PET_COLS), names(PET_COLS))
  pet[, sid := norm_id(sid)]
  if (!"BF2" %in% pet$Study) stop("no Study == 'BF2' rows -- label renamed?")
  pet <- pet[Study == "BF2" & !excluded %in% 1]   # NA is "no flag", not excluded
  pet[, Visit := as.numeric(Visit)]
  if (anyDuplicated(pet[, .(sid, Visit)]))
    stop("several rows per sid x Visit -- pick one per visit before modelling")
  num <- c("amyloid_status", "centiloid", "tau_pet", "thickness", "mmse", "mpacc")
  pet[, (num) := lapply(.SD, function(v) suppressWarnings(as.numeric(v))), .SDcols = num]

  # Tau-PET QC: mask Exclude_tauPET == 1 only, the rule
  # projects/proteogenomic_network/code/02_prepare_clinical.R documents for the Braak I-IV
  # composite.
  qc <- as.data.table(readxl::read_excel(TAU_QC, col_types = "text"))
  bad <- norm_id(qc$sid)[suppressWarnings(as.numeric(qc$Exclude_tauPET)) %in% 1]
  pet[sid %in% bad, tau_pet := NA_real_]

  base <- pet[Visit == 0]
  if (!all(base$amyloid_status %in% c(0, 1, NA))) stop("amyloid_status is not 0/1")
  # AD dementia vs cognitively unimpaired (Normal or SCD); everyone else is NA.
  base[, ad_vs_cu := fcase(trimws(diagnosis) == "Dementia AD", 1,
                           trimws(cog_status) %in% c("Normal", "SCD"), 0,
                           default = NA_real_)]
  if (base[ad_vs_cu == 1, .N] == 0 || base[ad_vs_cu == 0, .N] == 0)
    stop("ad_vs_cu has an empty group -- diagnosis labels changed?")
  # e4 dose from the two-digit genotype code (34 -> 1, 44 -> 2).
  code <- suppressWarnings(as.integer(base$apoe))
  base[, apoe_e4_dose := ifelse(is.na(code), NA_real_,
                                nchar(gsub("[^4]", "", as.character(code))))]

  # Visit is the nominal visit year (0, 2, 4, ...), used as time since baseline.
  long <- pet[!is.na(mmse) | !is.na(mpacc), .(sid, years = Visit, mmse, mpacc)]
  message(sprintf("phenotypes: %d BF2 baseline rows; %d cognition visits on %d sids",
                  nrow(base), nrow(long), uniqueN(long$sid)))
  list(base = base[, .(sid, amyloid_status, ad_vs_cu, centiloid, tau_pet, thickness,
                       mmse, mpacc, apoe_e4_dose)],
       long = long)
}

# GWAS covariates (age, sex, pc1..pc10, mean_protein_<tag>) keyed by sid.
load_covar <- function() {
  map <- fread(SID_IID, colClasses = "character")[, .(sid = norm_id(sid), IID = norm_iid(IID))]
  if (anyDuplicated(map$sid) || anyDuplicated(map$IID)) stop("sid <-> IID map is not 1:1")
  cov <- fread(COVAR, colClasses = list(character = "IID"))
  cov[, IID := norm_iid(IID)]
  cov[, intersect(c("FID", "#FID"), names(cov)) := NULL]
  if (anyDuplicated(cov$IID)) stop("duplicate IIDs in ", COVAR)
  n0 <- nrow(cov)
  cov <- merge(map, cov, by = "IID")
  message(sprintf("covariates: %d of %d GWAS samples map to a sid", nrow(cov), n0))
  cov
}

# INT protein values keyed by sid, one column per protein_assay. The .phe files hold
# exactly the pQTL GWAS samples, so these are the values the pQTLs were found in.
load_proteins <- function(prot_info, cov) {
  out <- cov[, .(sid, IID)]
  for (tag in unique(prot_info$panel_tag)) {
    f <- file.path(PROT_DIR, paste0("bf2_", tag, ".phe"))
    if (!file.exists(f)) stop("missing protein file: ", f)
    want <- prot_info[panel_tag == tag, protein_assay]
    hdr <- names(fread(f, nrows = 0))
    miss <- setdiff(want, hdr)
    if (length(miss)) stop(f, " lacks ", length(miss), " assays, e.g. ", miss[1])
    ph <- fread(f, select = c("IID", want), colClasses = list(character = "IID"))
    ph[, IID := norm_iid(IID)]
    if (anyDuplicated(ph$IID)) stop("duplicate IIDs in ", f)
    ph[, (want) := lapply(.SD, int), .SDcols = want]
    n0 <- nrow(out)
    out <- merge(out, ph, by = "IID", all.x = TRUE)
    if (nrow(out) != n0) stop("protein merge changed the row count for ", tag)
  }
  out[, IID := NULL]
  out
}

# Distinct fine-mapped signals (S1 credible-set leads deduplicated by lib/cs_signals.R),
# with the assay's panel tag, outside the MHC. APOE-region signals are flagged, not
# dropped: each analysis decides what to do with them.
atlas_signals <- function() {
  kept <- fread("output/00_kept_assays.csv")
  s1 <- fread("output/S1_pqtl_finemapping.csv.gz")[
    is_cs_lead == TRUE,
    .(cohort, protein_assay, cs_id, cs_key, fluid, platform, gene, ensembl_gene_id, is_cis,
      snp, chr, pos, ea, beta, se, p, shared = cs_shared_with_other_fluid)]
  sig <- distinct_signals(s1, ".")
  n0 <- nrow(sig)
  sig <- merge(sig, kept[, .(protein_assay, panel_tag = tag)], by = "protein_assay")
  if (nrow(sig) != n0) stop(n0 - nrow(sig), " signals have an assay missing from 00_kept_assays.csv")
  mhc <- fread("code/lib/mhc_hg19.range", header = FALSE,
               col.names = c("chr", "start", "end", "label"))
  in_mhc <- sig$chr == mhc$chr[1] & sig$pos >= mhc$start[1] & sig$pos <= mhc$end[1]
  message(sprintf("signals: %d distinct, %d in the MHC dropped", n0, sum(in_mhc)))
  sig <- sig[!in_mhc]
  sig[, apoe_region := in_apoe_region(chr, pos)]
  sig
}

# Dosage of a counted allele per SNP, keyed by normalized IID, one column per SNP.
# `snps` is data.table(snp, ea) with one row per SNP; a caller that needs another
# allele for the same SNP flips it (2 - dosage). Files go to WORK_DIR, never the repo.
# max_missing: the share of SNPs allowed to be absent from the export (the pfile was
# re-imputed after the GWAS, so a lead can be gone); they are listed and left out of
# the returned table. Above it, or at the default 0, a missing SNP stops the run.
export_dosages <- function(snps, name, max_missing = 0) {
  plink2 <- Sys.which("plink2")
  if (plink2 == "") stop("plink2 not on PATH")
  if (anyDuplicated(snps$snp)) stop("export_dosages: one counted allele per SNP")
  prefix <- file.path(WORK_DIR, name)
  writeLines(snps$snp, paste0(prefix, ".snps"))
  fwrite(snps[, .(snp, ea)], paste0(prefix, ".allele"), sep = "\t", col.names = FALSE)
  rc <- system2(plink2, c("--pfile", PFILE, "--extract", paste0(prefix, ".snps"),
                          "--export", "A", "--export-allele", paste0(prefix, ".allele"),
                          "--out", prefix), stdout = FALSE, stderr = FALSE)
  if (rc != 0 || !file.exists(paste0(prefix, ".raw")))
    stop("plink2 dosage export failed; see ", prefix, ".log")
  raw <- fread(paste0(prefix, ".raw"), colClasses = list(character = "IID"))
  want <- paste0(snps$snp, "_", snps$ea)
  miss <- setdiff(want, names(raw))
  if (length(miss) > max_missing * length(want))
    stop(length(miss), " SNPs missing from the dosage export (absent from the pfile or ",
         "counted allele differs), e.g. ", miss[1])
  if (length(miss)) message(sprintf("dosages: %d of %d SNPs missing from the export, left out: %s",
                                    length(miss), length(want), paste(miss, collapse = ", ")))
  keep <- want %in% names(raw)
  dose <- raw[, c("IID", want[keep]), with = FALSE]
  setnames(dose, want[keep], snps$snp[keep])
  dose[, IID := norm_iid(IID)]
  dose
}

# Lead dosages for a signal table (one row per signal, columns snp and ea): one counted
# allele per SNP, flip = TRUE where a signal's ea is the other allele (use 2 - dosage).
# Signals whose lead is missing from the export are dropped (export_dosages lists them).
# Returns list(sig, dose).
signal_dosages <- function(sig, name, max_missing) {
  snps <- unique(sig[, .(snp, ea)], by = "snp")
  dose <- export_dosages(snps, name, max_missing = max_missing)
  sig <- sig[snp %in% names(dose)]
  sig[snps, flip := ea != i.ea, on = "snp"]
  list(sig = sig, dose = dose)
}

# Covariate matrix of dat[, cols] through model.matrix, so a character or factor column
# (sex) is dummy-coded instead of turning the whole matrix into text; NA rows are kept
# so that rows stay aligned with the caller's outcome and dosage vectors.
cov_matrix <- function(dat, cols) {
  mf <- model.frame(~ ., dat[, ..cols], na.action = na.pass)
  model.matrix(~ ., mf)[, -1, drop = FALSE]
}

# Exposures of 11a/11c, keyed by IID, one column per exposure. Default: dosages of
# `snps_spec` ("rsid:counted_allele,..."). With `proteins` ("assay,...") the exposures
# are those protein levels (the GWAS-input INT), each with its panel's mean_protein as
# an extra covariate. Returns list(dose, snps = data.table(snp, ea), extra_cov).
load_exposures <- function(cov, proteins, snps_spec, name) {
  if (nzchar(proteins)) {
    want <- strsplit(proteins, ",")[[1]]
    kept <- fread("output/00_kept_assays.csv")[, .(protein_assay, panel_tag = tag)]
    pinfo <- kept[protein_assay %in% want]
    if (nrow(pinfo) != length(want)) stop("PROTEINS not all in 00_kept_assays.csv")
    dose <- merge(cov[, .(sid, IID)], load_proteins(pinfo, cov), by = "sid")[, sid := NULL]
    snps <- data.table(snp = pinfo$protein_assay, ea = NA_character_)
    extra_cov <- setNames(lapply(pinfo$panel_tag, function(t) intersect(paste0("mean_protein_", t), names(cov))),
                          pinfo$protein_assay)
  } else {
    spec <- strsplit(snps_spec, ",")[[1]]
    snps <- data.table(snp = sub(":.*", "", spec), ea = sub(".*:", "", spec))
    if (any(snps$snp == snps$ea)) stop("SNPS must be rsid:allele, e.g. rs5848:T")
    dose <- export_dosages(snps, name, max_missing = 0.01)   # a lead can be gone from the re-imputed pfile
    snps <- snps[snp %in% names(dose)]
    if (!nrow(snps)) stop("none of the requested SNPs is in the pfile")
    extra_cov <- setNames(rep(list(character()), nrow(snps)), snps$snp)
  }
  list(dose = dose, snps = snps, extra_cov = extra_cov)
}

# Covariates for one protein: the GWAS used mean_protein_<tag> when the covariate
# file had it, so do the same and say which way it went.
protein_covars <- function(tag, cov) {
  mp <- paste0("mean_protein_", tag)
  c("age", "sex", if (mp %in% names(cov)) mp)
}
