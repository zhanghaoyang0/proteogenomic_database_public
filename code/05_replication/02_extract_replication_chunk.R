#!/usr/bin/env Rscript
library(data.table)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop("Usage: Rscript 02_extract_replication_chunk.R candidates requests start end out_dir")
}
candidate_file = path.expand(args[1])
request_file = path.expand(args[2])
idx_start = as.integer(args[3])
idx_end = as.integer(args[4])
out_dir = path.expand(args[5])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidates = fread(candidate_file)
requests = fread(request_file)
idx_end = min(idx_end, nrow(requests))
if (idx_start > idx_end) quit(status = 0L)
requests = requests[idx_start:idx_end]

source(path.expand("code/lib/common.R"))   # REPO_ROOT
source(file.path(REPO_ROOT, "code/lib/gwas_functions.R"))   # write_atomic

rows = vector("list", nrow(requests))
status = vector("list", nrow(requests))
keep_cols = c("CHR", "POS", "SNP", "EA", "NEA", "EA_FREQ", "BETA", "SE", "P", "N")

for (i in seq_len(nrow(requests))) {
  req = requests[i]
  snps = if (req$cohort == "bf2") {
    unique(candidates[bf2_assay == req$assay, SNP])
  } else {
    unique(candidates[bf1_assay == req$assay & bf1_assay_available == TRUE, SNP])
  }
  message("Request ", req$request_id, ": ", req$cohort, " ", req$assay,
          " (", length(snps), " SNPs)")
  dt = fread(req$clean_path, select = keep_cols, showProgress = FALSE)
  hit = dt[SNP %chin% snps]
  if (nrow(hit)) {
    hit[, `:=`(request_id = req$request_id, cohort = req$cohort, assay = req$assay)]
    rows[[i]] = hit
  }
  status[[i]] = data.table(
    request_id = req$request_id,
    cohort = req$cohort,
    assay = req$assay,
    requested_snps = length(snps),
    found_snps = uniqueN(hit$SNP)
  )
}

stem = sprintf("extract_%05d_%05d", idx_start, idx_end)
extracted = rbindlist(rows, use.names = TRUE, fill = TRUE)
if (!ncol(extracted)) {
  extracted = data.table(
    CHR = integer(), POS = integer(), SNP = character(), EA = character(),
    NEA = character(), EA_FREQ = numeric(), BETA = numeric(), SE = numeric(),
    P = numeric(), N = numeric(), request_id = integer(), cohort = character(),
    assay = character()
  )
}
write_atomic(extracted,
             file.path(out_dir, paste0(stem, ".tsv.gz")))
write_atomic(rbindlist(status),
             file.path(out_dir, paste0(stem, ".status.tsv.gz")))
