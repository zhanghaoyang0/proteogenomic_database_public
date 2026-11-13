library(data.table)
library(dplyr)
project_root <- "."
setwd(project_root)
source("code/lib/gwas_functions.R")

# Reference afreq (CHR, POS, ID, REF, ALT, ALT_FREQS)
file_freq <- path.expand("~/db/shared_data/ref_geno/1kg_v3/EUR.afreq")
freq <- fread(file_freq)
setnames(freq, c("#CHROM", "ID"), c("CHR", "SNP"))

# Iwaki reference (SNP, RSID, CHR, START, REF, ALT, MAF, ...)
iwaki_ref_file <- file.path(path.expand("~/db/shared_data/gwas/other/raw"), "iwaki_reference.txt.gz")
iwaki_ref <- fread(iwaki_ref_file)
setnames(iwaki_ref, "START", "POS")

raw_dir <- path.expand("~/db/shared_data/gwas/other/raw")
clean_dir <- path.expand("~/db/shared_data/gwas/other/clean")
# PRS-only bases (entry `out_dir = prs_base_dir`): kept out of other/clean, which every
# database step (trait_info, clump, finemap, MR/coloc, app) reads in full.
prs_base_dir <- path.expand("~/db/shared_data/gwas/other/prs_base")
plot_dir <- path.expand("~/db/shared_data/gwas/other/plot")

assert_merge_rows <- function(before, after, label, allow_drop = FALSE) {
  if (!nrow(after)) stop(label, " merge produced zero rows")
  if (nrow(after) > before) stop(label, " merge expanded rows: ", before, " -> ", nrow(after))
  if (!allow_drop && nrow(after) != before) stop(label, " merge changed rows: ", before, " -> ", nrow(after))
}

# List: raw filename -> trait_id; optional rename = c(standard_name = "raw_col_name", ...)
gwas_list <- list(
  list(raw = "2247_1.mtag.txt.gz", trait = "arhi_kalra_2020", n = NULL,
       rename = c(SNP = "snpid", EA = "a1", NEA = "a2", BETA = "mtag_beta", SE = "mtag_se", P = "mtag_pval", N = "n")),
  list(raw = "HD_EA_gwas_sumstats.txt.gz", trait = "hd_wells_2019", n = 250389,
       rename = c(POS = "BP", EA = "ALLELE1", NEA = "ALLELE0", EA_FREQ = "A1FREQ")),
  list(raw = "HAID_EA_gwas_sumstats.txt.gz", trait = "haid_wells_2019", n = 253918,
       rename = c(POS = "BP", EA = "ALLELE1", NEA = "ALLELE0", EA_FREQ = "A1FREQ")),
  list(raw = "ieu-b-7.vcf.gz", trait = "pd_ipdgc_2019", n = 482730,
       rename = c(CHR = "#CHROM", SNP = "ID", EA = "ALT", NEA = "REF")),
  list(raw = "GCST90001390_buildGRCh38.tsv.gz", trait = "lbd_chia_2021", n = 6618, hg38 = TRUE,
       rename = c(
         SNP = "variant_id",
         P = "p_value",
         CHR = "chromosome",
         POS = "base_pair_location",
         EA = "effect_allele",
         NEA = "other_allele",
         EA_FREQ = "effect_allele_frequency",
         BETA = "beta",
         SE = "standard_error"
       )),
  list(raw = "cont_HY.txt.gz", trait = "hy_iwaki_2019", n = NULL,
       rename = NULL),
  list(raw = "cont_MOCA.txt.gz", trait = "moca_iwaki_2019", n = NULL,
       rename = NULL),
  list(raw = "cont_UPDRS_scaled.txt.gz", trait = "updrs_iwaki_2019", n = NULL,
      rename = NULL),
  # Bellenguez 2022 AD (GCST90027158): N varies per variant -> n_cases + n_controls
  list(raw = "GCST90027158_buildGRCh38.tsv.gz", trait = "ad_bellenguez_2022", n = NULL, hg38 = TRUE,
       n_sum = c("n_cases", "n_controls"),
       rename = c(SNP = "variant_id", P = "p_value", CHR = "chromosome", POS = "base_pair_location",
                  EA = "effect_allele", NEA = "other_allele", EA_FREQ = "effect_allele_frequency",
                  BETA = "beta", SE = "standard_error")),
  # Uffelmann 2026 AD (PGC-ALZ3, Nature, cncr.nl): main/case_control EUR only -- no proxy
  # cases (their betas are x2-scaled), ancestry matching the 1KG EUR LD. hg19; 23andMe
  # removed, so N is below the paper's. N varies per variant -> n_case + n_control.
  list(raw = "twent_case_control_eur_neff0.6_nsumstats1.txt.gz", trait = "ad_uffelmann_2026", n = NULL,
       n_sum = c("n_case", "n_control"),
       rename = c(SNP = "rsid", CHR = "chromosome", POS = "base_pair_location",
                  EA = "effect_allele", NEA = "other_allele", EA_FREQ = "effect_allele_frequency",
                  BETA = "beta", SE = "standard_error", P = "p_value")),
  # Kunkle 2019 AD (GCST007511), IGAP stage 1: 21,982 clinical/autopsy cases + 41,944
  # controls, no proxy cases. hg19; no EAF (imputed from 1KG); indels without an rsid
  # are named chrN:pos:I/D and drop out at the EAF step.
  # PRS base only (AMYPAD harmonised PRS): overlaps Bellenguez at 1/8 of its N, not in the database.
  list(raw = "Kunkle_etal_Stage1_results.txt.gz", trait = "ad_kunkle_2019", n = 63926,
       out_dir = prs_base_dir,
       rename = c(CHR = "Chromosome", POS = "Position", SNP = "MarkerName",
                  EA = "Effect_allele", NEA = "Non_Effect_allele", BETA = "Beta", P = "Pvalue")),
  # Jansen 2022 CSF Abeta42 / p-tau (GCST90129599 / GCST90129600): no rsid, Z only.
  # rsid comes from the Catalog's harmonised file; its stats are not used, because it
  # flips EA/NEA and EAF but leaves Zscore unflipped.
  # n_cap = the study's total N: ~500 SNPs per trait, at the top loci (APOE, GMNC, BIN1,
  # CR1, ...), carry replication samples (N up to 13k), which made SuSiE-RSS's
  # "N >= 80% of max" filter drop ~99% of each window. Capping keeps those loci; their
  # Z is ~1.27x what the capped N implies, so treat Jansen coloc as approximate.
  list(raw = "GCST90129599_buildGRCh38.tsv.gz", trait = "csf_ab42_jansen_2022", n = NULL, hg38 = TRUE,
       rsid_map = "GCST90129599.h.tsv.gz", from_z = TRUE, n_cap = 8074,
       rename = c(P = "p_value", CHR = "chromosome", POS = "base_pair_location",
                  EA = "effect_allele", NEA = "other_allele", EA_FREQ = "effect_allele_frequency")),
  list(raw = "GCST90129600_buildGRCh38.tsv.gz", trait = "csf_ptau_jansen_2022", n = NULL, hg38 = TRUE,
       rsid_map = "GCST90129600.h.tsv.gz", from_z = TRUE, n_cap = 7798,
       rename = c(P = "p_value", CHR = "chromosome", POS = "base_pair_location",
                  EA = "effect_allele", NEA = "other_allele", EA_FREQ = "effect_allele_frequency"))
)

# Attach rsid (SNP) by CHR/POS + unordered allele pair from a GWAS Catalog harmonised
# file; keys mapping to more than one rsid are dropped. Stops if under half match.
add_rsid_by_pos <- function(dat, map_file) {
  m <- fread(map_file, select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "rsid"))
  m <- m[!is.na(rsid) & grepl("^rs", rsid)]
  m[, `:=`(CHR = as.character(chromosome), POS = as.integer(base_pair_location),
           a1 = pmin(toupper(effect_allele), toupper(other_allele)),
           a2 = pmax(toupper(effect_allele), toupper(other_allele)))]
  m <- m[, .(SNP = rsid[1], n_rs = uniqueN(rsid)), by = .(CHR, POS, a1, a2)][n_rs == 1, !"n_rs"]
  dat[, `:=`(CHR = as.character(CHR), POS = as.integer(POS), a1 = pmin(EA, NEA), a2 = pmax(EA, NEA))]
  n_in <- nrow(dat)
  dat <- merge(dat, m, by = c("CHR", "POS", "a1", "a2"))
  dat[, c("a1", "a2") := NULL]
  message("  rsid by position: ", nrow(dat), " / ", n_in, " rows matched")
  if (nrow(dat) < 0.5 * n_in) stop("rsid map matched under half the rows: ", map_file)
  dat
}

# Z-only summary statistics -> per-allele effect in SD units:
# SE = 1 / sqrt(2 f (1 - f) (N + Z^2)), BETA = Z * SE
z_to_beta_se <- function(dat) {
  dat[, SE := 1 / sqrt(2 * EA_FREQ * (1 - EA_FREQ) * (N + Zscore^2))]
  dat[, BETA := Zscore * SE]
  dat[, Zscore := NULL]
  dat
}

# Parse BETA, SE, P from a VCF sample column based on FORMAT
parse_vcf_sample_from_format <- function(dat, sample_col) {
  if (!("FORMAT" %in% names(dat)) || !(sample_col %in% names(dat))) {
    return(dat)
  }
  vals <- strsplit(dat[[sample_col]], ":", fixed = TRUE)
  fmt_fields <- strsplit(dat$FORMAT[1], ":", fixed = TRUE)[[1]]
  k <- length(fmt_fields)
  parts <- t(vapply(vals, function(x) { length(x) <- k; x }, character(k)))
  pos_ES <- match("ES", fmt_fields)
  pos_SE <- match("SE", fmt_fields)
  pos_LP <- match("LP", fmt_fields)
  if (!is.na(pos_ES) && !is.na(pos_SE) && !is.na(pos_LP)) {
    dat <- dat %>% mutate(
      BETA = as.numeric(parts[, pos_ES]),
      SE   = as.numeric(parts[, pos_SE]),
      P    = 10^(-as.numeric(parts[, pos_LP]))
    )
  }
  dat
}

for (item in gwas_list) {
  file_in <- file.path(raw_dir, item$raw)
  out_dir <- if (is.null(item$out_dir)) clean_dir else item$out_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_out <- file.path(out_dir, paste0(item$trait, ".tsv.gz"))
  if (!file.exists(file_in)) stop("raw file missing for ", item$trait, ": ", file_in)
  # Skip only a non-empty clean file at least as new as its raw input, as 02_clean_sumstats.R
  # does. A changed gwas_list entry is not seen here: delete the trait's clean file (README).
  if (file.exists(file_out) && file.size(file_out) > 0 && file.mtime(file_out) >= file.mtime(file_in)) {
    message("Skip ", item$trait, ": output is up to date")
    next
  }
  message("Processing ", item$trait, " from ", item$raw)
  is_vcf <- grepl("\\.vcf\\.gz$", file_in, ignore.case = TRUE)
  if (is_vcf) {
    gwas <- fread(cmd = paste("zgrep -v '^##'", shQuote(file_in)), sep = "\t")
  } else {
    gwas <- fread(file_in)
  }
  gwas_clean <- as.data.table(gwas)
  if (!is.null(item$rename)) {
    gwas_clean <- gwas_clean %>% rename(!!!item$rename)
  }
  if (is_vcf) {
    sample_col <- sub("\\.vcf\\.gz$", "", item$raw)
    gwas_clean <- parse_vcf_sample_from_format(gwas_clean, sample_col)
  }
  setDT(gwas_clean)

  # harmonize_to_1kg compares EA/NEA to the 1KG REF/ALT as-is, so lower case never matches
  for (a in intersect(c("EA", "NEA"), names(gwas_clean))) set(gwas_clean, j = a, value = toupper(gwas_clean[[a]]))
  if (!is.null(item$n_sum)) gwas_clean[, N := rowSums(.SD), .SDcols = item$n_sum]
  if (!is.null(item$rsid_map)) gwas_clean <- add_rsid_by_pos(gwas_clean, file.path(raw_dir, item$rsid_map))
  if (isTRUE(item$from_z)) gwas_clean <- z_to_beta_se(gwas_clean)
  if (!is.null(item$n_cap)) {   # after from_z, so SE still matches the published Z
    message("  N capped at ", item$n_cap, ": ", sum(gwas_clean$N > item$n_cap, na.rm = TRUE), " rows")
    gwas_clean[, N := pmin(N, item$n_cap)]
  }

  if (!is.null(item$n)) {
    gwas_clean <- gwas_clean %>% mutate(N = item$n)
  } else if (!("N" %in% names(gwas_clean))) {
    gwas_clean <- gwas_clean %>% mutate(N = NA_real_)
  }

  # Iwaki traits with association stats: fill CHR/POS/alleles/frequency from iwaki_reference, then harmonize to 1KG
  if (item$trait %in% c("hy_iwaki_2019", "moca_iwaki_2019", "updrs_iwaki_2019")) {
    cols_ref <- intersect(c("SNP", "CHR", "POS", "REF", "ALT", "MAF"), names(iwaki_ref))
    iwaki_ref_sub <- iwaki_ref[, ..cols_ref]
    n_pre_merge <- nrow(gwas_clean)
    gwas_clean <- merge(gwas_clean, iwaki_ref_sub, by = "SNP", all.x = TRUE)
    assert_merge_rows(n_pre_merge, gwas_clean, paste(item$trait, "Iwaki reference"))
    if (!("EA" %in% names(gwas_clean))) gwas_clean[, EA := ALT]
    if (!("NEA" %in% names(gwas_clean))) gwas_clean[, NEA := REF]
    if (!("EA_FREQ" %in% names(gwas_clean)) && "MAF" %in% names(gwas_clean)) gwas_clean[, EA_FREQ := MAF]
    # Drop REF/ALT/MAF so harmonize_to_1kg can reuse REF/ALT from 1KG freq without name conflicts
    drop_cols <- intersect(c("REF", "ALT", "MAF"), names(gwas_clean))
    if (length(drop_cols) > 0) gwas_clean[, (drop_cols) := NULL]
  }

  # hg38 traits: take hg19 CHR/POS from 1KG by SNP (no liftover). Rows whose rsid is not
  # in 1KG are dropped: kept, they would carry hg38 coords into an hg19 table and could
  # match an unrelated 1KG variant by position in harmonize_to_1kg.
  if (isTRUE(item$hg38) && "SNP" %in% names(gwas_clean)) {
    kg_pos <- freq[, .(SNP, CHR_1kg = CHR, POS_1kg = POS)]
    n_in <- nrow(gwas_clean)
    gwas_clean <- merge(gwas_clean, kg_pos, by = "SNP")
    assert_merge_rows(n_in, gwas_clean, paste(item$trait, "hg38-to-1KG"), allow_drop = TRUE)
    message("hg38: hg19 CHR/POS from 1KG by SNP; ", nrow(gwas_clean), " / ", n_in, " rows in 1KG")
    gwas_clean[, `:=`(CHR = as.integer(CHR_1kg), POS = as.integer(POS_1kg))]
    gwas_clean[, c("CHR_1kg", "POS_1kg") := NULL]
  }

  if (!"CHR" %in% names(gwas_clean) || !"POS" %in% names(gwas_clean)) {
    message("CHR or POS missing, add from freq by SNP")
    n_pre_merge <- nrow(gwas_clean)
    gwas_clean <- merge(gwas_clean, freq[, .(SNP, CHR, POS)], by = "SNP", all.x = TRUE)
    assert_merge_rows(n_pre_merge, gwas_clean, paste(item$trait, "1KG position"))
  }
  if (!"EA_FREQ" %in% names(gwas_clean)) {
    message("EA_FREQ missing, impute from freq")
    gwas_clean <- harmonize_to_1kg(gwas_clean, freq, impute_frq = TRUE)
  } else {
    message("EA_FREQ found, harmonize to 1kg")
    gwas_clean <- harmonize_to_1kg(gwas_clean, freq)
  }
  gwas_clean <- round_gwas_clean(gwas_clean)
  if (!nrow(gwas_clean)) stop("No rows after harmonization for ", item$trait)
  write_atomic(gwas_clean, file_out)
  message("  Written ", nrow(gwas_clean), " rows -> ", file_out)
}

# No PNGs unless SKIP_PLOT=FALSE: the app draws from manhattan_data/.
skip_plot <- as.logical(Sys.getenv("SKIP_PLOT", "TRUE"))
if (!skip_plot) {
  dir_manhattan <- file.path(plot_dir, "manhattan")
  dir_qq        <- file.path(plot_dir, "qq")
  dir.create(dir_manhattan, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_qq,        recursive = TRUE, showWarnings = FALSE)
  for (item in gwas_list) {
    if (!is.null(item$out_dir)) next   # PRS-only base, not in the app
    trait         <- item$trait
    clean_path    <- file.path(clean_dir, paste0(trait, ".tsv.gz"))
    out_manhattan <- file.path(dir_manhattan, paste0(trait, ".png"))
    out_qq        <- file.path(dir_qq,        paste0(trait, ".png"))
    plot_manhattan_qq(clean_path, out_manhattan, out_qq, skip_if_exists = TRUE)
  }
}
