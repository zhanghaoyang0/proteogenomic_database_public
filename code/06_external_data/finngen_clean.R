#!/usr/bin/env Rscript

# FinnGen R12 summary stats are GRCh38; EUR.afreq is hg19 — lift positions before harmonize.

library(data.table)
library(dplyr)
library(rtracklayer)

args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(repo_root)
source("code/lib/gwas_functions.R")

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(i, default = NULL) {
  if (length(args) >= i && nzchar(args[i])) args[i] else default
}

MODE <- get_arg(1, Sys.getenv("FINNGEN_MODE", "both"))
IDX_START <- as.integer(get_arg(2, Sys.getenv("IDX_START", "1")))
IDX_END <- as.integer(get_arg(3, Sys.getenv("IDX_END", "99999999")))
RAW_DIR <- get_arg(4, Sys.getenv("RAW_DIR", "~/db/shared_data/gwas/finngen/raw/"))
CLEAN_DIR <- get_arg(5, Sys.getenv("CLEAN_DIR", "~/db/shared_data/gwas/finngen/clean/"))
PLOT_DIR   <- get_arg(6, Sys.getenv("PLOT_DIR", "~/db/shared_data/gwas/finngen/plot/"))
# No PNGs unless SKIP_PLOT=FALSE: the app draws from manhattan_data/, and a
# default of drawing made every plain run redraw the whole set.
SKIP_PLOT  <- as.logical(Sys.getenv("SKIP_PLOT", "TRUE"))
INFO_DIR <- get_arg(7, Sys.getenv("INFO_DIR", "~/db/shared_data/gwas/finngen/info/"))
FILE_FREQ <- Sys.getenv("FILE_FREQ", path.expand("~/db/shared_data/ref_geno/1kg_v3/EUR.afreq"))
CHAIN_HG38_TO_HG19 <- Sys.getenv("CHAIN_HG38_TO_HG19", path.expand("~/db/shared_data/ref_geno/liftover/hg38ToHg19.over.chain"))

if (!MODE %in% c("both", "clean", "plot")) {
  stop("MODE must be one of: both, clean, plot")
}

if (is.na(IDX_START) || is.na(IDX_END) || IDX_START < 1 || IDX_END < IDX_START) {
  stop("Invalid IDX_START/IDX_END")
}

# ============================================================
# download info
# ============================================================
# cd ~/db/shared_data/gwas/finngen
# wget -c -P info https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv # fin
# wget -c -P info https://storage.googleapis.com/finngen-public-data-r12/meta_analysis/ukbb/finngen_R12_meta_analysis_mapping_with_definitions.tsv # fin, ukbb
# wget -c -P info https://storage.googleapis.com/finngen-public-data-r12/meta_analysis/mvp_ukbb/FinnGen_R12_MVP_UKBB_mapping.tsv # fin, ukbb, mvp
# wget -c -P info https://storage.googleapis.com/finngen-public-data-r12/meta_analysis/mvp_ukbb/FinnGen_R12_MVP_UKBB_manifest.tsv # fin, ukbb, mvp

# # ============================================================
# # make urls list to download
# # ============================================================
# info_fin = fread(paste0(INFO_DIR, "finngen_R12_manifest.tsv"))
# info_fin_ukbb = fread(paste0(INFO_DIR, "finngen_R12_meta_analysis_mapping_with_definitions.tsv"))
# info_fin_ukbb_mvp = fread(paste0(INFO_DIR, "FinnGen_R12_MVP_UKBB_mapping.tsv"))
# info_fin_ukbb_mvp_manifest = fread(paste0(INFO_DIR, "FinnGen_R12_MVP_UKBB_manifest.tsv"))

# mk_download <- function(dt, pheno_col, name_col, source_name, url_expr) {
#   dt[, .(
#     phenocode = get(pheno_col),
#     name = get(name_col),
#     source = source_name,
#     url = eval(url_expr)
#   )]
# }

# download_info <- rbindlist(list(
#   mk_download(info_fin, "phenocode", "phenotype", "fin", quote(path_https)),
#   mk_download(
#     info_fin_ukbb_mvp_manifest, "fg_phenotype", "name", "fin_ukbb_mvp",
#     quote(paste0(
#       "https://storage.googleapis.com/finngen-public-data-r12/meta_analysis/mvp_ukbb/summary_stats_filtered/",
#       fg_phenotype, "_meta_out_filtered.tsv.gz"
#     ))
#   ),
#   mk_download(
#     info_fin_ukbb, "fg_phenotype", "name", "fin_ukbb",
#     quote(paste0(
#       "https://storage.googleapis.com/finngen-public-data-r12/meta_analysis/ukbb/summary_stats/",
#       "finngen_R12_", fg_phenotype, "_meta_out.tsv.gz"
#     ))
#   )
# ), use.names = TRUE, fill = TRUE)

# source_priority <- c("fin_ukbb_mvp", "fin_ukbb", "fin")
# download_info[, source_rank := match(source, source_priority)]
# setorder(download_info, phenocode, source_rank)
# download_info <- download_info[, .SD[1], by = phenocode][, source_rank := NULL]
# fwrite(download_info, paste0(INFO_DIR, "download_info.csv"))

# # save url lists by source
# fwrite(download_info[source == "fin_ukbb_mvp", .(url)],
#        file = paste0(INFO_DIR, "urls_fin_ukbb_mvp.txt"),
#        col.names = FALSE,
#        quote = FALSE)
# fwrite(download_info[source == "fin_ukbb", .(url)],
#        file = paste0(INFO_DIR, "urls_fin_ukbb.txt"),
#        col.names = FALSE,
#        quote = FALSE)
# fwrite(download_info[source == "fin", .(url)],
#        file = paste0(INFO_DIR, "urls_fin.txt"),
#        col.names = FALSE,
#        quote = FALSE)

# # ============================================================
# # download gwas summary stats with bash
# # ============================================================
# cd ~/db/shared_data/gwas/finngen/
# mkdir -p raw/fin_ukbb_mvp raw/fin_ukbb raw/fin
# cat info/urls_fin_ukbb_mvp.txt | xargs -d '\n' -n 1 -P 5 wget -c -P ./raw/fin_ukbb_mvp/
# cat info/urls_fin_ukbb.txt | grep ALZ | xargs -d '\n' -n 1 -P 5 wget -c -P ./raw/fin_ukbb/
# cat info/urls_fin.txt | xargs -d '\n' -n 1 -P 5 wget -c -P ./raw/fin/

# ============================================================
# helpers
# ============================================================

build_raw_manifest <- function(raw_dir, sources_clean) {
  manifest_list <- list()
  for (i in seq_along(sources_clean)) {
    source_name <- sources_clean[i]
    source_dir <- file.path(path.expand(raw_dir), source_name)
    if (!dir.exists(source_dir)) {
      message("Skip source directory not found: ", source_name)
      next
    }

    files_raw <- sort(list.files(
      source_dir,
      pattern = "\\.tsv\\.gz$",
      full.names = TRUE
    ))
    if (identical(Sys.getenv("FINNGEN_REVERSE", "0"), "1")) {
      files_raw <- rev(files_raw)
      message("FINNGEN_REVERSE=1: reverse order in ", source_name)
    }
    if (length(files_raw) == 0) next

    manifest_list[[length(manifest_list) + 1]] <- data.table(
      source_rank = i,
      source_name = source_name,
      file_in = files_raw
    )
  }

  if (length(manifest_list) == 0) {
    return(data.table(source_rank = integer(), source_name = character(), file_in = character()))
  }

  manifest <- rbindlist(manifest_list, use.names = TRUE)
  manifest[, trait_raw := sub("\\.tsv\\.gz$", "", basename(file_in))]
  manifest[, trait_base := sub("_meta_out_filtered$", "_meta_", trait_raw)]
  manifest[, trait_base := sub("_meta_out$", "_meta_", trait_base)]
  manifest[, trait := paste0(trait_base, source_name)]
  manifest[]
}

build_clean_manifest <- function(clean_dir) {
  files <- sort(list.files(path.expand(clean_dir), pattern = "\\.tsv\\.gz$", full.names = TRUE))
  if (length(files) == 0) {
    return(data.table(file_clean = character(), trait = character()))
  }
  data.table(
    file_clean = files,
    trait = sub("\\.tsv\\.gz$", "", basename(files))
  )
}

slice_manifest <- function(manifest, idx_start, idx_end, label) {
  n <- nrow(manifest)
  if (n == 0) {
    message("No ", label, " files found.")
    return(manifest)
  }
  idx_end <- min(idx_end, n)
  if (idx_start > n) {
    message("Nothing to do for ", label, ": idx_start=", idx_start, " > total=", n)
    return(manifest[0])
  }
  chunk <- manifest[idx_start:idx_end]
  message("Selected ", label, " idx ", idx_start, "-", idx_end, " of ", n, " (n=", nrow(chunk), ")")
  chunk
}

clean_one_trait <- function(file_in, source_name, trait_raw, trait, n_lookup_fin_ukbb, n_lookup_fin_ukbb_mvp,
                            chain_hg38_to_hg19, freq, clean_dir) {
  file_out <- file.path(path.expand(clean_dir), paste0(trait, ".tsv.gz"))
  # Skip only a non-empty clean file at least as new as its raw input, as 02_clean_sumstats.R
  # does. A changed setting (N lookup, filters) is not seen here: delete the clean file (README).
  if (file.exists(file_out) && file.size(file_out) > 0 && file.mtime(file_out) >= file.mtime(file_in)) {
    message("Skipping ", trait, ": output is up to date")
    return(file_out)
  }

  message("Processing ", trait)
  dt <- fread(file_in)

  if (source_name == "fin_ukbb_mvp") {
    phenocode <- sub("_meta_out_filtered$", "", trait_raw)
    n_sample <- as.numeric(n_lookup_fin_ukbb_mvp[phenocode])
    af_col <- "fg_af_alt"
    beta_col <- "all_inv_var_meta_beta"
    se_col <- "all_inv_var_meta_sebeta"
    p_col <- "all_inv_var_meta_p"
    n_col <- "all_meta_N"
  } else {
    phenocode <- sub("^finngen_R12_", "", trait_raw)
    phenocode <- sub("_meta_out$", "", phenocode)
    n_sample <- as.numeric(n_lookup_fin_ukbb[phenocode])
    af_col <- "FINNGEN_af_alt"
    beta_col <- "all_inv_var_meta_beta"
    se_col <- "all_inv_var_meta_sebeta"
    p_col <- "all_inv_var_meta_p"
    n_col <- "all_meta_N"
  }

  if (!all(c(af_col, beta_col, se_col, p_col, n_col) %in% names(dt))) {
    stop("Missing required columns for ", trait)
  }

  gwas <- dt[, .(
    CHR = as.integer(gsub("^chr", "", as.character(`#CHR`), ignore.case = TRUE)),
    POS = as.integer(POS),
    SNP = if ("rsid" %in% names(dt)) as.character(rsid) else as.character(SNP),
    EA = as.character(ALT),
    NEA = as.character(REF),
    EA_FREQ = as.numeric(get(af_col)),
    BETA = as.numeric(get(beta_col)),
    SE = as.numeric(get(se_col)),
    P = as.numeric(get(p_col)),
    N = if (is.na(n_sample)) as.numeric(get(n_col)) else n_sample
  )]

  gwas <- gwas[!is.na(CHR) & CHR >= 1 & CHR <= 22]
  gwas <- gwas[!is.na(POS) & POS > 0]
  gwas <- gwas[!is.na(SNP) & nzchar(SNP)]

  gwas <- lift_to_hg19(gwas, chain_hg38_to_hg19)
  if (nrow(gwas) == 0) {
    stop("No rows after liftover for ", trait)
  }

  gwas_clean <- round_gwas_clean(harmonize_to_1kg(gwas, freq))
  if (!nrow(gwas_clean)) stop("No rows after harmonization for ", trait)
  write_atomic(gwas_clean, file_out)
  file_out
}

plot_one_trait <- function(file_clean, trait, plot_dir) {
  if (!file.exists(file_clean)) {
    return(invisible(NULL))
  }
  out_manhattan <- file.path(path.expand(plot_dir), "manhattan", paste0(trait, ".png"))
  out_qq <- file.path(path.expand(plot_dir), "qq", paste0(trait, ".png"))
  plot_manhattan_qq(file_clean, out_manhattan, out_qq, skip_if_exists = TRUE)
}

# ============================================================
# clean gwas
# ============================================================

dir.create(path.expand(CLEAN_DIR), recursive = TRUE, showWarnings = FALSE)
if (!SKIP_PLOT) dir.create(file.path(path.expand(PLOT_DIR), "manhattan"), recursive = TRUE, showWarnings = FALSE)
if (!SKIP_PLOT) dir.create(file.path(path.expand(PLOT_DIR), "qq"), recursive = TRUE, showWarnings = FALSE)

sources_clean <- c("fin_ukbb_mvp", "fin_ukbb")

# Cases + controls per phenotype. An NA count is a cohort that did not contribute, so it
# adds nothing; but every count NA gives NA, not 0, so clean_one_trait falls back to the
# per-variant N instead of writing N = 0.
sum_na_num <- function(x) {
  m <- as.data.frame(lapply(x, as.numeric))
  s <- rowSums(m, na.rm = TRUE)
  s[rowSums(!is.na(m)) == 0L] <- NA
  s
}

freq <- NULL
chain_hg38_to_hg19 <- NULL
n_lookup_fin_ukbb <- NULL
n_lookup_fin_ukbb_mvp <- NULL

if (MODE %in% c("both", "clean")) {
  freq <- fread(FILE_FREQ)
  chain_hg38_to_hg19 <- rtracklayer::import(CHAIN_HG38_TO_HG19, format = "chain")
  info_fin_ukbb <- fread(file.path(path.expand(INFO_DIR), "finngen_R12_meta_analysis_mapping_with_definitions.tsv"))
  info_fin_ukbb_mvp <- fread(file.path(path.expand(INFO_DIR), "FinnGen_R12_MVP_UKBB_manifest.tsv"))

info_fin_ukbb[, N_sample := sum_na_num(.SD), .SDcols = c(
  "fg_n_cases", "fg_n_controls", "ukbb_n_cases", "ukbb_n_controls"
)]
info_fin_ukbb_mvp[, N_sample := sum_na_num(.SD), .SDcols = c(
  "fg_n_cases", "fg_n_controls", "ukbb_n_cases", "ukbb_n_controls",
  "MVP_AFR_n_cases", "MVP_EUR_n_cases", "MVP_AMR_n_cases",
  "MVP_AFR_n_controls", "MVP_EUR_n_controls", "MVP_AMR_n_controls"
)]

n_lookup_fin_ukbb <- setNames(info_fin_ukbb$N_sample, info_fin_ukbb$fg_phenotype)
n_lookup_fin_ukbb_mvp <- setNames(info_fin_ukbb_mvp$N_sample, info_fin_ukbb_mvp$fg_phenotype)
}

# ============================================================
# Manhattan + QQ
# ============================================================

if (MODE %in% c("both", "clean")) {
  raw_manifest <- build_raw_manifest(RAW_DIR, sources_clean)
  raw_chunk <- slice_manifest(raw_manifest, IDX_START, IDX_END, "raw")

  if (nrow(raw_chunk) > 0) {
    cleaned_files <- character(nrow(raw_chunk))
    clean_errors <- character()
    for (i in seq_len(nrow(raw_chunk))) {
      cleaned_files[i] <- tryCatch(
        clean_one_trait(
          file_in = raw_chunk$file_in[i],
          source_name = raw_chunk$source_name[i],
          trait_raw = raw_chunk$trait_raw[i],
          trait = raw_chunk$trait[i],
          n_lookup_fin_ukbb = n_lookup_fin_ukbb,
          n_lookup_fin_ukbb_mvp = n_lookup_fin_ukbb_mvp,
          chain_hg38_to_hg19 = chain_hg38_to_hg19,
          freq = freq,
          clean_dir = CLEAN_DIR
        ),
        error = function(e) {
          clean_errors <<- c(clean_errors, paste0(raw_chunk$trait[i], ": ", conditionMessage(e)))
          message("ERROR: ", tail(clean_errors, 1L))
          NA_character_
        }
      )
    }

    if (MODE == "both") {
      for (i in seq_len(nrow(raw_chunk))) {
        if (!is.na(cleaned_files[i]) && !SKIP_PLOT) {
          plot_one_trait(cleaned_files[i], raw_chunk$trait[i], PLOT_DIR)
        }
      }
    }
    if (length(clean_errors)) {
      stop(length(clean_errors), " trait(s) failed:\n", paste(clean_errors, collapse = "\n"))
    }
  }
}

if (MODE == "plot") {
  clean_manifest <- build_clean_manifest(CLEAN_DIR)
  clean_chunk <- slice_manifest(clean_manifest, IDX_START, IDX_END, "clean")
  if (nrow(clean_chunk) > 0) {
    for (i in seq_len(nrow(clean_chunk))) {
      if (!SKIP_PLOT) plot_one_trait(clean_chunk$file_clean[i], clean_chunk$trait[i], PLOT_DIR)
    }
  }
}

message("Done.")
