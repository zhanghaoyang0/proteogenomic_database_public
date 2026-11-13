library(data.table)

args = commandArgs(trailingOnly = TRUE)
qc_dir = if (length(args) >= 1L) path.expand(args[1]) else path.expand("~/db/shared_data/gwas/bf1/clean/protein_qc")
clean_dir = if (length(args) >= 2L) path.expand(args[2]) else path.expand("~/db/shared_data/gwas/bf1/clean/protein")

args_all = commandArgs(trailingOnly = FALSE)
repo_root = normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", args_all, value = TRUE))), "..", ".."))
source(file.path(repo_root, "code/lib/gwas_functions.R"))
expected = bf1_gwas_panels(file.path(repo_root, "info/bf_protein/qc_panel_summary.csv"))

qc_files = list.files(qc_dir, pattern = "^bf1__.*\\.qc\\.tsv$", full.names = TRUE)
clean_files = list.files(clean_dir, pattern = "^bf1__.*\\.tsv\\.gz$", full.names = TRUE)
if (length(qc_files) != sum(expected$expected_assays)) {
  stop("Expected ", sum(expected$expected_assays), " trait QC files, found ", length(qc_files))
}
if (length(clean_files) != sum(expected$expected_assays)) {
  stop("Expected ", sum(expected$expected_assays), " clean files, found ", length(clean_files))
}

trait_qc = rbindlist(lapply(qc_files, fread), use.names = TRUE, fill = TRUE)
required = c("trait", "panel", "n_variants", "n_min", "n_median", "n_max",
             "invalid_beta", "invalid_se", "duplicated_snps", "invalid_p", "lambda_gc")
missing_cols = setdiff(required, names(trait_qc))
if (length(missing_cols)) stop("QC files lack columns: ", paste(missing_cols, collapse = ", "))
if (anyDuplicated(trait_qc$trait)) stop("Duplicated traits in QC files")
if (anyNA(trait_qc$panel)) stop("One or more traits could not be assigned to a known BF1 panel")
clean_traits = sub("^bf1__(.*)\\.tsv\\.gz$", "\\1", basename(clean_files))
if (!setequal(clean_traits, trait_qc$trait)) {
  stop("Clean-file and trait-QC manifests do not contain the same assays")
}

panel_qc = trait_qc[, .(
  clean_assays = .N,
  empty_assays = sum(n_variants == 0),
  variants = as.numeric(sum(n_variants)),
  n_min = as.numeric(min(n_min, na.rm = TRUE)),
  n_median = as.numeric(median(n_median, na.rm = TRUE)),
  n_max = as.numeric(max(n_max, na.rm = TRUE)),
  invalid_beta = as.numeric(sum(invalid_beta)),
  invalid_se = as.numeric(sum(invalid_se)),
  duplicated_snps = as.numeric(sum(duplicated_snps)),
  invalid_p = as.numeric(sum(invalid_p)),
  lambda_gc_min = min(lambda_gc, na.rm = TRUE),
  lambda_gc_median = median(lambda_gc, na.rm = TRUE),
  lambda_gc_max = max(lambda_gc, na.rm = TRUE)
), by = panel]
panel_qc = merge(expected, panel_qc, by = "panel", all.x = TRUE, sort = FALSE)
panel_qc[, count_ok := clean_assays == expected_assays]

setorder(trait_qc, panel, trait)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(trait_qc, file.path(qc_dir, "bf1_pqtl_trait_qc.tsv.gz"), sep = "\t")
fwrite(panel_qc, file.path(qc_dir, "bf1_pqtl_panel_qc.tsv"), sep = "\t")
print(panel_qc)

if (!all(panel_qc$count_ok)) stop("BF1 panel clean counts do not match expectations")
if (any(trait_qc$n_variants == 0)) stop("One or more BF1 clean outputs contain no variants")
if (trait_qc[, sum(invalid_beta + invalid_se + duplicated_snps + invalid_p)] != 0) {
  stop("Clean-output QC found invalid or duplicated values; inspect bf1_pqtl_trait_qc.tsv.gz")
}
message("BF1 pQTL clean QC passed for all ", nrow(trait_qc), " assays")
