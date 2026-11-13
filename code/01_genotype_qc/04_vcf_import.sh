#!/usr/bin/env bash
# Post-imputation QC for TOPMed-imputed genotype data.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/qc_common.sh"

GENO_DIR="${GENO_DIR:-$HOME/db/secure_data/bf/array/bf2}"
IMPUTED_DIR="${IMPUTED_DIR:-${GENO_DIR}/imputed}"
OUT_DIR="${OUT_DIR:-${GENO_DIR}/imputed_qc}"

R2_MIN="${R2_MIN:-0.8}"
MAF_MIN="${MAF_MIN:-0.001}"
THREADS="${THREADS:-5}"
MEM_MB="${MEM_MB:-14000}"
CLEANUP_INTERMEDIATE="${CLEANUP_INTERMEDIATE:-0}"
RSCRIPT_BIN="${RSCRIPT_BIN:-}"

DATASET="${DATASET:-$(basename "${GENO_DIR}" | sed 's/_geno$//')}"
CHR_PREFIX="${CHR_PREFIX:-${DATASET}_topmed_chr}"

RAW_DIR="${OUT_DIR}/raw"
LIST_DIR="${RAW_DIR}/01_lists"
PGEN_DIR="${RAW_DIR}/01_pgen_by_chr"
QC_DIR="${RAW_DIR}/01_topmed_qc"
REPORT_DIR="${OUT_DIR}/report"
LOG_DIR="${OUT_DIR}/log"
TOPMED_VCF_DIR="${IMPUTED_DIR}"

mkdir -p "${LIST_DIR}" "${PGEN_DIR}" "${QC_DIR}" "${REPORT_DIR}" "${LOG_DIR}"

has_valid_topmed_pair() {
  local dir="$1"
  local chr="$2"
  [[ -s "${dir}/chr_${chr}/chr${chr}.dose.vcf.gz" && -s "${dir}/chr_${chr}/chr${chr}.info.gz" ]] || return 1
  gzip -t "${dir}/chr_${chr}/chr${chr}.dose.vcf.gz" "${dir}/chr_${chr}/chr${chr}.info.gz" >/dev/null 2>&1
}


check_inputs() {
  need_cmd gzip
  need_cmd plink2
  setup_rscript

  for chr in {1..22}; do
    if has_valid_topmed_pair "${IMPUTED_DIR}" "${chr}"; then
      continue
    fi

    echo "ERROR: missing or invalid ${IMPUTED_DIR}/chr${chr}.dose.vcf.gz or ${IMPUTED_DIR}/chr${chr}.info.gz" >&2
    echo "Please unzip TOPMed results first, then rerun this script." >&2
    exit 1
  done
}

prepare_topmed_files() {
  echo "==> Step 1: use TOPMed VCF/INFO files from ${TOPMED_VCF_DIR}"
}

make_info_keep_lists() {
  echo "==> Step 2: build per-chromosome INFO/R2 + MAF keep lists..."
  "${RSCRIPT_CMD[@]}" - "${TOPMED_VCF_DIR}" "${LIST_DIR}" "${REPORT_DIR}/01_topmed_info_summary.tsv" "${R2_MIN}" "${MAF_MIN}" <<'RSCRIPT'
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
info_dir <- args[[1]]
list_dir <- args[[2]]
summary_out <- args[[3]]
r2_min <- as.numeric(args[[4]])
maf_min <- as.numeric(args[[5]])

pick_col <- function(dt, candidates, required = TRUE) {
  hit <- intersect(candidates, names(dt))
  if (length(hit) > 0) return(hit[[1]])
  if (required) stop("None of these columns found: ", paste(candidates, collapse = ", "))
  NA_character_
}

info_value <- function(x, key) {
  pattern <- paste0("(^|;)", key, "=([^;]+)")
  out <- sub(pattern, "\\2", regmatches(x, regexpr(pattern, x)))
  hit <- grepl(pattern, x)
  ans <- rep(NA_character_, length(x))
  ans[hit] <- out
  ans
}

summaries <- list()
for (chr in 1:22) {
  f <- file.path(info_dir, paste0("chr_", chr), paste0("chr", chr, ".info.gz"))
  dt <- fread(cmd = paste("gzip -dc", shQuote(f)), showProgress = FALSE)

  id_col <- pick_col(dt, c("SNP", "ID", "MarkerName", "VAR_ID", "Name"))
  if ("INFO" %in% names(dt) && !any(c("Rsq", "R2", "Imputation_R2", "ER2", "LooRsq") %in% names(dt))) {
    r2_col <- "INFO:R2"
    maf_col <- "INFO:MAF/AF"
    r2 <- suppressWarnings(as.numeric(info_value(dt[["INFO"]], "R2")))
    maf <- suppressWarnings(as.numeric(info_value(dt[["INFO"]], "MAF")))
    if (all(is.na(maf))) {
      af <- suppressWarnings(as.numeric(info_value(dt[["INFO"]], "AF")))
      maf <- pmin(af, 1 - af)
    }
  } else {
    r2_col <- pick_col(dt, c("Rsq", "R2", "INFO", "Imputation_R2", "ER2", "LooRsq"))
    maf_col <- pick_col(dt, c("MAF", "MAF(1)", "ALT_Frq", "AF", "ALT_FREQ", "Freq1"), required = FALSE)
    r2 <- suppressWarnings(as.numeric(dt[[r2_col]]))
    if (!is.na(maf_col)) {
      maf_raw <- suppressWarnings(as.numeric(dt[[maf_col]]))
      maf <- if (grepl("frq|freq|af", maf_col, ignore.case = TRUE)) pmin(maf_raw, 1 - maf_raw) else maf_raw
    } else {
      maf <- rep(NA_real_, nrow(dt))
    }
  }

  keep <- !is.na(r2) & r2 >= r2_min
  if (!all(is.na(maf))) keep <- keep & !is.na(maf) & maf >= maf_min

  keep_ids <- unique(as.character(dt[[id_col]][keep]))
  keep_ids <- keep_ids[!is.na(keep_ids) & keep_ids != "" & keep_ids != "."]
  fwrite(data.table(ID = keep_ids), file.path(list_dir, paste0("chr", chr, ".info_r2_maf.keep")), col.names = FALSE)

  summaries[[chr]] <- data.table(
    CHR = chr,
    n_info = nrow(dt),
    n_keep = length(keep_ids),
    id_col = id_col,
    r2_col = r2_col,
    maf_col = ifelse(is.na(maf_col), "", maf_col),
    r2_min_observed = suppressWarnings(min(r2, na.rm = TRUE)),
    r2_median = suppressWarnings(median(r2, na.rm = TRUE)),
    r2_p10 = suppressWarnings(as.numeric(quantile(r2, 0.10, na.rm = TRUE))),
    r2_p90 = suppressWarnings(as.numeric(quantile(r2, 0.90, na.rm = TRUE)))
  )
}

fwrite(rbindlist(summaries), summary_out, sep = "\t")
message("INFO keep lists written for ", length(summaries), " chromosomes")
RSCRIPT

  # Fail loudly if R silently did nothing (conda-run-style stdin swallowing)
  assert_file "${REPORT_DIR}/01_topmed_info_summary.tsv" "make_info_keep_lists R script (Step 2)"
  assert_file "${LIST_DIR}/chr22.info_r2_maf.keep" "make_info_keep_lists R script (Step 2)"
}

import_chr_pfiles() {
  echo "==> Step 3: import dosage VCFs to per-chromosome pfiles (sequential)..."

  # Process in reverse order (22→1) to start with small chromosomes for easier debugging
  for chr in {22..1}; do
    local raw_prefix="${PGEN_DIR}/${CHR_PREFIX}${chr}_raw"
    local final_prefix="${PGEN_DIR}/${CHR_PREFIX}${chr}"

    # Skip if already built from this dose VCF. A pgen older than the VCF comes from
    # an earlier imputation; rebuild it instead of merging old genotypes.
    local dose_vcf="${TOPMED_VCF_DIR}/chr_${chr}/chr${chr}.dose.vcf.gz"
    if [[ -s "${final_prefix}.pgen" && -s "${final_prefix}.pvar" && -s "${final_prefix}.psam" \
          && "${final_prefix}.pgen" -nt "${dose_vcf}" ]]; then
      echo "  chr${chr}: already processed, skipping"
      continue
    fi

    echo "  chr${chr}: importing and standardizing..."

    plink2 \
      --vcf "${TOPMED_VCF_DIR}/chr_${chr}/chr${chr}.dose.vcf.gz" dosage=DS \
      --extract "${LIST_DIR}/chr${chr}.info_r2_maf.keep" \
      --max-alleles 2 \
      --snps-only just-acgt \
      --make-pgen \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${raw_prefix}" \
      > "${LOG_DIR}/01_import_chr${chr}.log" 2>&1

    plink2 \
      --pfile "${raw_prefix}" \
      --set-all-var-ids '@:#:$r:$a' \
      --new-id-max-allele-len 500 \
      --make-pgen \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${final_prefix}" \
      > "${LOG_DIR}/01_standardize_chr${chr}.log" 2>&1

    rm -f "${raw_prefix}".*
  done
}

merge_chr_pfiles() {
  echo "==> Step 4: merge chromosomes and apply final MAF filter..."
  : > "${LIST_DIR}/pmerge_list.txt"
  for chr in {2..22}; do
    echo "${PGEN_DIR}/${CHR_PREFIX}${chr}" >> "${LIST_DIR}/pmerge_list.txt"
  done

  plink2 \
    --pfile "${PGEN_DIR}/${CHR_PREFIX}1" \
    --pmerge-list "${LIST_DIR}/pmerge_list.txt" pfile \
    --rm-dup force-first \
    --make-pgen \
    --threads "${THREADS}" \
    --memory "${MEM_MB}" \
    --out "${QC_DIR}/topmed_merged" \
    > "${LOG_DIR}/01_merge.log" 2>&1

  plink2 \
    --pfile "${QC_DIR}/topmed_merged" \
    --maf "${MAF_MIN}" \
    --make-pgen \
    --threads "${THREADS}" \
    --memory "${MEM_MB}" \
    --out "${QC_DIR}/raw_imputed" \
    > "${LOG_DIR}/01_final_maf.log" 2>&1
}

make_reports() {
  echo "==> Step 5: write TOPMed QC summary reports..."
  "${RSCRIPT_CMD[@]}" - "${QC_DIR}/raw_imputed" "${REPORT_DIR}" <<'RSCRIPT'
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
imputed <- args[[1]]
report_dir <- args[[2]]

read_pvar <- function(prefix) {
  x <- fread(paste0(prefix, ".pvar"), select = c("#CHROM", "POS", "ID", "REF", "ALT"))
  setnames(x, c("CHR", "POS", "ID", "REF", "ALT"))
  x[, CHR := sub("^chr", "", as.character(CHR))]
  x[, key := paste(CHR, POS, toupper(REF), toupper(ALT), sep = ":")]
  x
}

read_psam <- function(prefix) {
  x <- fread(paste0(prefix, ".psam"), check.names = FALSE)
  iid_col <- if ("IID" %in% names(x)) {
    "IID"
  } else if ("#IID" %in% names(x)) {
    "#IID"
  } else {
    names(x)[[2]]
  }
  ids <- as.character(x[[iid_col]])
  ids <- sub("^0+", "", ids)
  ids[ids == ""] <- "0"
  unique(ids)
}

imp_pvar <- read_pvar(imputed)
imp_iid <- read_psam(imputed)

sample_summary <- data.table(
  dataset = "imputed_qc",
  n_samples = length(imp_iid)
)
fwrite(sample_summary, file.path(report_dir, "01_sample_summary.tsv"), sep = "\t")

variant_summary <- data.table(
  dataset = "imputed_qc",
  n_variants = nrow(imp_pvar),
  n_ids = uniqueN(imp_pvar$ID),
  n_keys = uniqueN(imp_pvar$key)
)
fwrite(variant_summary, file.path(report_dir, "01_variant_summary.tsv"), sep = "\t")

chr_counts <- imp_pvar[, .N, by = CHR][, dataset := "imputed_qc"]
chr_counts[, chr_order := suppressWarnings(as.integer(CHR))]
setorder(chr_counts, dataset, chr_order)
chr_counts[, chr_order := NULL]
fwrite(chr_counts, file.path(report_dir, "01_variant_counts_by_chr.tsv"), sep = "\t")
message("Sample summary: ", nrow(sample_summary), " row(s)")
RSCRIPT

  assert_file "${REPORT_DIR}/01_sample_summary.tsv" "make_reports R script (Step 5)"
  assert_file "${REPORT_DIR}/01_variant_summary.tsv" "make_reports R script (Step 5)"
}
cleanup_intermediate() {
  if [[ "${CLEANUP_INTERMEDIATE}" != "1" ]]; then
    return 0
  fi

  echo "==> Cleanup: remove per-chromosome pfiles..."
  rm -rf "${PGEN_DIR}"
  rm -f "${QC_DIR}/topmed_merged.pgen" \
        "${QC_DIR}/topmed_merged.pvar" \
        "${QC_DIR}/topmed_merged.psam" \
        "${QC_DIR}/topmed_merged.log" \
        "${QC_DIR}/topmed_merged-merge.pgen" \
        "${QC_DIR}/topmed_merged-merge.pvar" \
        "${QC_DIR}/topmed_merged-merge.psam" \
        "${QC_DIR}/topmed_merged-merge.log"
}

main() {
  check_inputs
  prepare_topmed_files
  make_info_keep_lists
  import_chr_pfiles
  merge_chr_pfiles
  make_reports
  cleanup_intermediate

  echo "Done."
  echo "Final pfile: ${QC_DIR}/raw_imputed"
  echo "Reports: ${REPORT_DIR}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main "$@"
fi
