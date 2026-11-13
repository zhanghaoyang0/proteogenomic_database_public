library(data.table)
library(dplyr)
args_all <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- dirname(normalizePath(script_file))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(repo_root)
source("code/lib/gwas_functions.R")

# Args: raw_dir clean_dir idx_start idx_end [plot_dir] [ensg_file] [qc_dir] [manifest_file]
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript 02_clean_sumstats.R raw_dir clean_dir idx_start idx_end [plot_dir] [ensg_file] [qc_dir] [manifest_file]")
raw_dir   = args[1]
clean_dir = args[2]
idx_start = as.integer(args[3])
idx_end   = as.integer(args[4])
plot_dir  = if (length(args) >= 5 && nzchar(args[5])) args[5] else NULL
ensg_file = if (length(args) >= 6 && nzchar(args[6])) args[6] else file.path(repo_root, "info/bf_protein/ensg/assay_ensg.csv")
qc_dir    = if (length(args) >= 7 && nzchar(args[7])) args[7] else NULL
manifest_file = if (length(args) >= 8 && nzchar(args[8])) path.expand(args[8]) else NULL

file_freq = path.expand("~/db/shared_data/ref_geno/1kg_v3/EUR.afreq")
freq = fread(file_freq)

dir.create(path.expand(clean_dir), recursive = TRUE, showWarnings = FALSE)
if (!is.null(qc_dir)) dir.create(path.expand(qc_dir), recursive = TRUE, showWarnings = FALSE)

cohort = sub(".*/(bf[12])/.*", "\\1", raw_dir)
prefix = if (grepl("^bf[12]$", cohort)) paste0(cohort, "__") else ""

ensg_cis = if (file.exists(ensg_file)) fread(ensg_file) else NULL
is_pheno = grepl("/pheno/?$", raw_dir)

# Collect current plink2 zst outputs (linear + logistic), or use an immutable
# manifest so a directory can continue receiving uploads while this job runs.
if (!is.null(manifest_file)) {
  if (!file.exists(manifest_file)) stop("Manifest not found: ", manifest_file)
  manifest_paths = readLines(manifest_file, warn = FALSE)
  manifest_paths = manifest_paths[nzchar(manifest_paths)]
  if (!length(manifest_paths)) stop("Manifest is empty: ", manifest_file)
  if (anyDuplicated(manifest_paths)) stop("Manifest contains duplicated paths: ", manifest_file)
  if (any(!file.exists(manifest_paths))) stop("One or more manifest files no longer exist")
  manifest_names = basename(manifest_paths)
  files_linear = manifest_names[grepl("\\.glm\\.linear\\.zst$", manifest_names)]
  files_logistic = manifest_names[grepl("\\.glm\\.logistic(\\.hybrid)?\\.zst$", manifest_names)]
  if (length(files_linear) + length(files_logistic) != length(manifest_names)) {
    stop("Manifest contains unsupported filenames")
  }
} else {
  files_linear   = list.files(path.expand(raw_dir), pattern = "\\.glm\\.linear\\.zst$")
  files_logistic = list.files(path.expand(raw_dir), pattern = "\\.glm\\.logistic(\\.hybrid)?\\.zst$")
}
files_zst      = c(files_linear, files_logistic)

if (!length(files_zst)) {
  stop("No .glm.linear.zst or .glm.logistic*.zst files found in ", raw_dir)
}

file_tbl = data.table(
  file_raw = files_zst,
  file_type = c(rep("linear", length(files_linear)),
                rep("logistic", length(files_logistic)))
)
file_tbl[, raw_trait := sub("\\.glm\\.(linear|logistic(\\.hybrid)?)\\.zst$", "", file_raw)]
if (is_pheno) {
  file_tbl[, trait := raw_trait]
} else {
  file_tbl[, trait := sub("^[^.]+\\.", "", raw_trait)]
}
files_raw = file_tbl$file_raw
file_types = file_tbl$file_type
traits = file_tbl$trait

idx_end   = min(idx_end, length(traits))
files_raw = files_raw[idx_start:idx_end]
file_types = file_types[idx_start:idx_end]
traits    = traits[idx_start:idx_end]
message("Processing traits ", idx_start, "-", idx_end, " (n=", length(traits), ") prefix='", prefix, "'")

bf1_panels = bf1_gwas_panels()$panel

panel_from_trait = function(trait) {
  hits = bf1_panels[endsWith(trait, paste0("__", bf1_panels))]
  if (length(hits) == 1L) hits else NA_character_
}

clean_qc = function(gwas_clean, trait) {
  expected_cols = c("CHR", "POS", "SNP", "EA", "NEA", "EA_FREQ", "BETA", "SE", "P", "N")
  if (!identical(names(gwas_clean), expected_cols)) {
    stop("Unexpected clean columns for ", trait, ": ", paste(names(gwas_clean), collapse = " "))
  }
  p = as.numeric(gwas_clean$P)
  # P = 0 is plink2 underflow on very strong cis signals (BETA and SE finite): BF2 has
  # it, and BF1 has 140 rows in 7 assays since its N roughly doubled (2026-11-04).
  # lambda_gc_from_p clips it.
  p_ok = is.finite(p) & p >= 0 & p <= 1
  data.table(
    trait = trait,
    panel = panel_from_trait(trait),
    n_variants = nrow(gwas_clean),
    n_min = if (nrow(gwas_clean)) min(gwas_clean$N, na.rm = TRUE) else NA_real_,
    n_median = if (nrow(gwas_clean)) median(gwas_clean$N, na.rm = TRUE) else NA_real_,
    n_max = if (nrow(gwas_clean)) max(gwas_clean$N, na.rm = TRUE) else NA_real_,
    invalid_beta = sum(!is.finite(gwas_clean$BETA)),
    invalid_se = sum(!is.finite(gwas_clean$SE) | gwas_clean$SE <= 0),
    duplicated_snps = sum(duplicated(gwas_clean$SNP)),
    invalid_p = sum(!p_ok),
    lambda_gc = if (any(p_ok)) lambda_gc_from_p(p[p_ok]) else NA_real_
  )
}

for (i in seq_along(traits)) {
  trait    = traits[i]
  ftype    = file_types[i]
  file_in  = path.expand(file.path(raw_dir, files_raw[i]))
  file_out = path.expand(file.path(clean_dir, paste0(prefix, trait, ".tsv.gz")))
  file_qc  = if (!is.null(qc_dir)) path.expand(file.path(qc_dir, paste0(prefix, trait, ".qc.tsv"))) else NULL
  # Skip only a non-empty clean file at least as new as its raw input: a newer
  # raw .zst (rsync -a keeps the Mac mtime) means the GWAS was rerun.
  output_ok = file.exists(file_out) && !is.na(file.info(file_out)$size) && file.info(file_out)$size > 0 &&
    file.mtime(file_out) >= file.mtime(file_in)
  if (output_ok) {
    if (!is.null(file_qc) && !file.exists(file_qc)) {
      message("Clean output exists; creating missing QC: ", trait)
      gwas_clean = fread(file_out, showProgress = FALSE)
      write_atomic(clean_qc(gwas_clean, trait), file_qc)
    } else {
      message("Skipping ", trait)
    }
    next
  }
  message("Processing ", trait, " [", ftype, "]")

  gwas = fread(cmd = paste("zstd -d -c", shQuote(file_in)))
  setnames(gwas, c("#CHROM", "ID",  "A1", "OMITTED", "A1_FREQ", "OBS_CT"),
                 c("CHR",    "SNP", "EA", "NEA",     "EA_FREQ", "N"))
  if ("TEST"    %in% names(gwas)) gwas = gwas[TEST == "ADD"]
  if ("ERRCODE" %in% names(gwas)) gwas = gwas[ERRCODE == "."]
  gwas[, CHR := as.integer(CHR)]

  gwas_clean = round_gwas_clean(harmonize_to_1kg(gwas, freq))
  if (!nrow(gwas_clean)) stop("0 rows after harmonisation: ", trait, " (", file_in, ")")
  write_atomic(gwas_clean, file_out)
  if (!is.null(file_qc)) write_atomic(clean_qc(gwas_clean, trait), file_qc)
}

if (!is.null(plot_dir)) {
  dir_manhattan = file.path(path.expand(plot_dir), "manhattan")
  dir_qq        = file.path(path.expand(plot_dir), "qq")
  dir.create(dir_manhattan, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_qq,        recursive = TRUE, showWarnings = FALSE)

  for (i in seq_along(traits)) {
    trait   = traits[i]
    file_in = path.expand(file.path(clean_dir, paste0(prefix, trait, ".tsv.gz")))
    if (!file.exists(file_in)) next
    out_manhattan = file.path(dir_manhattan, paste0(prefix, trait, ".png"))
    out_qq        = file.path(dir_qq,        paste0(prefix, trait, ".png"))
    cis_info = NULL
    if (!is.null(ensg_cis)) {
      row = ensg_cis[protein_assay == sub("^bf[12]__", "", trait)]
      if (nrow(row) > 0)
        cis_info = list(chr = row$gene_chr[1], start = row$cis_start[1], end = row$cis_end[1])
    }
    plot_manhattan_qq(file_in, out_manhattan, out_qq, skip_if_exists = TRUE, cis_info = cis_info)
  }
}
