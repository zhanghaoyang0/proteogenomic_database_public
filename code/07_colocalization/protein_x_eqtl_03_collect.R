#!/usr/bin/env Rscript

# Collect a protein x molecular-QTL coloc chain into its SQLite ledger, then export the
# flat result tsv that downstream (report scripts, Shiny app) already reads.
#
# Source: staging/job_*.tsv, the current run's output, one file per job (removed
# after import). The one-off migration of legacy data/{id}.tsv.gz files was removed
# on 2026-09-30, so an old data/ dir can no longer leak back into the ledger.
#
# Import is an upsert keyed on exposure_id, so re-running a subset of proteins replaces
# exactly their rows. Multiple testing is not applied here; filter by PP4 downstream.

SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_db.R"))
source(file.path(SCRIPT_DIR, "functions_coloc.R"))   # add_component_flags

# Which molecular-QTL chain to build for. The three differ only in where their
# results live -- the protein selection, the ledger semantics and the todo logic are
# identical -- so they share one script rather than three 89%-identical copies.
#   mb   MetaBrain bulk brain eQTL
#   gtex GTEx multi-tissue eQTL
#   sb   SingleBrain single-nucleus (cell-type) brain eQTL
CHAIN <- commandArgs(trailingOnly = TRUE)[1L]
CHAIN_LABEL <- c(mb = "MetaBrain", gtex = "GTEx eQTL",
                 sb = "SingleBrain single-nucleus eQTL")[CHAIN]
RESOURCE <- c(gtex = "gtex", sb = "singlebrain", mb = "metabrain")[CHAIN]
if (is.na(CHAIN_LABEL)) stop("usage: <script> <mb|gtex|sb>  (got: ",
                             if (is.na(CHAIN)) "nothing" else CHAIN, ")")

COLOC_DIR <- file.path(COLOC_ROOT, sprintf("protein_x_%s", CHAIN))
STAGE_DIR <- file.path(COLOC_DIR, "staging")
DB_FILE   <- file.path(COLOC_DIR, "coloc.db")
DB_TABLE  <- sprintf("coloc_protein_x_%s", CHAIN)
PROT_FULL <- file.path(COLOC_DIR, "info", "protein_list.csv")
RESULT    <- file.path(COLOC_DIR, "result", sprintf("coloc_protein_x_%s.tsv", CHAIN))
KEEP_STAGING <- nzchar(Sys.getenv("KEEP_STAGING", ""))
dir.create(dirname(RESULT), showWarnings = FALSE, recursive = TRUE)

t0  <- proc.time()[[3L]]

# ---- 1. import this run's staging files -------------------------------------------
msg("chain: protein x ", CHAIN_LABEL, " -> ", COLOC_DIR)
stage <- list.files(STAGE_DIR, pattern = "^job_.*\\.tsv$", full.names = TRUE)
if (length(stage)) {
  # A failed chunk blocks the afterok collect, so a resubmit can compute some
  # proteins into a second staging file: keep each protein's newest file only,
  # or db_upsert would append both copies.
  stage <- stage[order(file.mtime(stage))]
  d <- rbindlist(lapply(seq_along(stage), function(i) fread(stage[i])[, .src := i]), fill = TRUE)
  msg("staging: ", length(stage), " file(s), ", nrow(d), " rows, ",
      uniqueN(d$exposure_id), " proteins")
  d <- d[d[, .I[.src == max(.src)], by = exposure_id]$V1][, .src := NULL]
  msg("after keeping the newest file per protein: ", nrow(d), " rows")
  db_upsert(DB_FILE, DB_TABLE, d, "exposure_id")
  msg("imported staging -> ", DB_FILE)
} else {
  msg("no staging files in ", STAGE_DIR)
}

# ---- 3. export the flat table downstream reads -------------------------------------
all <- db_read_all(DB_FILE, DB_TABLE)
if (is.null(all) || !nrow(all)) stop("ledger is empty: ", DB_FILE)
msg("ledger: ", nrow(all), " rows, ", uniqueN(all$exposure_id), " proteins")

# Filter to the current exposure universe before exporting, the way the MR and
# protein_x_disease collects do. The ledger is append-only history: it still holds
# rows for exposures the selection has since dropped (the BF1 cohort, retired
# panels), and exporting it verbatim silently put them in the table downstream
# reads. Keeping them in the DB means a scope change is a re-collect, not a
# recompute -- COLOC_INCLUDE_BF1=1 in 01_pair brings BF1 back for free.
prot_universe <- fread(PROT_FULL, select = "exposure_id")$exposure_id
n_before <- nrow(all)
all <- all[exposure_id %chin% prot_universe]
msg("filtered to the current protein_list: ", n_before, " -> ", nrow(all), " rows, ",
    uniqueN(all$exposure_id), " proteins")
if (!nrow(all)) stop("nothing left after filtering to ", PROT_FULL)

cat("\nstatus breakdown:\n"); print(all[, .N, by = status][order(-N)])

# computed_at is ledger bookkeeping (01_pair staleness), not a result column.
e_src <- protein_src(all$exposure_id)
add_component_flags(all,   # before the drop below: it checks computed_at
  e_rds = protein_susie_path(e_src, all$exposure_id, lead_or_na(all$lead_snp)),
  e_fm  = finemap_path(e_src, all$exposure_id),
  o_rds = eqtl_susie_path(RESOURCE, all$tissue, all$ENSG, lead_or_na(all$o_lead)),
  o_fm  = eqtl_finemap_path(RESOURCE, all$tissue, all$ENSG),
  e_cis_only = TRUE)
if ("computed_at" %in% names(all)) all[, computed_at := NULL]

# ---- 4. coverage check before the table is replaced, as collect_chunked.R does --------
# The manuscript reads this table, so a partial run must not overwrite it.
full <- fread(PROT_FULL)
n_full <- uniqueN(full$exposure_id)   # rows != proteins if the list has duplicates
n_done <- length(intersect(full$exposure_id, unique(all$exposure_id)))
msg("coverage: ", round(100 * n_done / n_full, 2), "% (", n_done, "/", n_full, ")")
if (n_done < n_full && !nzchar(Sys.getenv("ALLOW_PARTIAL", "")))
  stop(n_full - n_done, " protein(s) not computed -- run them, or set ALLOW_PARTIAL=1")
fwrite(all, RESULT, sep = "\t")
msg("wrote -> ", RESULT)

if (length(stage) && !KEEP_STAGING) {
  file.remove(stage)
  msg("removed ", length(stage), " staging file(s) (set KEEP_STAGING=1 to keep)")
}
