#!/usr/bin/env Rscript

# Collect protein2disease MR results. Merge existing collected results with any
# new tmp_chunks/chunk_*.tsv rows, dedupe by (exposure_id, outcome_id), filter to
# the current mr_pair.csv.gz universe, then SPLIT by outcome source:
#   - mr_protein2disease.tsv        : FinnGen only (the manuscript table)
#   - mr_protein2disease_other.tsv  : Other outcomes (data others need; not in article)
# Raw ivw_p only — multiple-testing correction is deliberately NOT applied here
# (correction strategy undecided; do it downstream on the chosen universe).
# Cleanup tmp_chunks/ only if row counts look complete.

suppressPackageStartupMessages(library(data.table))

SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_db.R"))


MR_DIR     <- file.path(MR_ROOT, "protein2disease")
TMP_DIR    <- file.path(MR_DIR, "tmp_chunks")
PAIR_FILE  <- file.path(MR_DIR, "info", "mr_pair.csv.gz")
RESULT     <- file.path(MR_DIR, "result", "mr_protein2disease.tsv.gz")  # gzipped source of truth
RESULT_OTH <- file.path(MR_DIR, "result", "mr_protein2disease_other.tsv")  # plain (small; read as-is downstream)
DB_FILE    <- file.path(MR_DIR, "mr.db")
dir.create(dirname(RESULT), showWarnings = FALSE, recursive = TRUE)

t0  <- proc.time()[[3L]]

chunks <- list.files(TMP_DIR, pattern = "^chunk_.*\\.tsv$", full.names = TRUE)
result_files <- c(RESULT, RESULT_OTH)
result_files <- result_files[file.exists(result_files)]
if (!length(chunks) && !length(result_files)) {
  stop("No chunk files in ", TMP_DIR, " and no existing result files")
}
msg("found ", length(chunks), " chunk files and ", length(result_files), " existing result file(s)")

old <- if (length(result_files)) rbindlist(lapply(result_files, fread), fill = TRUE) else NULL
new <- if (length(chunks)) rbindlist(lapply(chunks, fread), fill = TRUE) else NULL
all <- rbindlist(list(old, new), fill = TRUE)
msg("rbinded: ", nrow(all), " rows")
all <- unique(all, by = c("exposure_id", "outcome_id"), fromLast = TRUE)
msg("after dedup: ", nrow(all), " rows")

cat("\nstatus breakdown:\n"); print(all[, .N, by = status][order(-N)])

# Split by outcome source (chunk rows carry no source -> look it up from mr_pair.csv.gz).
pairs     <- fread(PAIR_FILE, select = c("exposure_id", "outcome_id", "outcome_source"))
if (anyDuplicated(pairs[, .(exposure_id, outcome_id)]))
  stop("mr_pair.csv.gz has duplicate (exposure_id, outcome_id) keys")
pairs[, pair_key := paste(exposure_id, outcome_id, sep = "__")]
all[, pair_key := paste(exposure_id, outcome_id, sep = "__")]
all <- all[pair_key %in% pairs$pair_key]
msg("after filtering to current pair universe: ", nrow(all), " rows")

# A collected row (including an explicit failure status) is required for every
# current pair.  Do this before replacing results, updating the ledger or deleting
# chunks, so a partial array can always be resumed safely.
missing_keys <- setdiff(pairs$pair_key, all$pair_key)
extra_keys <- setdiff(all$pair_key, pairs$pair_key)
if (length(missing_keys) || length(extra_keys) || nrow(all) != nrow(pairs))
  stop("MR collection incomplete: rows=", nrow(all), "/", nrow(pairs),
       " missing=", length(missing_keys), " extra=", length(extra_keys),
       if (length(missing_keys)) paste0("; examples: ", paste(head(missing_keys, 3L), collapse = ", ")) else "")
all[, pair_key := NULL]

other_ids <- unique(pairs[outcome_source == "Other", outcome_id])
is_other  <- all$outcome_id %in% other_ids
all_main  <- all[!is_other]
all_oth   <- all[is_other]
msg("split: FinnGen=", nrow(all_main), " rows, Other=", nrow(all_oth), " rows")

# No multiple-testing correction here — write raw ivw_p; correct downstream.
# Atomic: these are the only copies of the results, and the ledger below says done.
source(file.path(REPO_ROOT, "code/lib/gwas_functions.R"))   # write_atomic
write_atomic(all_main, RESULT);     msg("wrote -> ", RESULT)
write_atomic(all_oth,  RESULT_OTH); msg("wrote -> ", RESULT_OTH)

# ---- ledger: record what has run, then it is safe to drop the chunks ----
con <- db_open(DB_FILE)
n_expo <- done_record(con, all[, .(exposure_id, outcome_id)], replace = TRUE)
msg("ledger updated: ", n_expo, " exposures (", db_done_n(con), " total on file)")
dbDisconnect(con)

pair_n <- nrow(pairs)
if (length(chunks)) {
  unlink(chunks)
  msg("cleaned up ", length(chunks),
      " chunk files after exact pair-key coverage was verified")
}
msg("coverage: 100% (", nrow(all), "/", pair_n, ")")
