#!/usr/bin/env Rscript

# SQLite ledger helpers shared by the coloc sub-project.
#
# Why: results used to be stored as one small file per unit of work (data/{id}.tsv.gz).
# That made "already done?" a directory listing, but cost one inode per protein and made
# collect pay a GPFS open() per file (measured: 149 s to read 5027 files totalling 0.6 MB).
#
# New layout per coloc type:
#   coloc.db            SQLite: one table, all result + status rows. The ledger.
#   staging/job_*.tsv   append-only per-JOB output, imported by collect then removed.
#
# The unit of work (key_col -- one column like exposure_id, or several like
# c("exposure_id","outcome_id")) is what "done" is keyed on -- NOT the row, since one
# unit yields many rows. Re-running a unit replaces all of its rows.
# Chunk boundaries never touch storage, so adding/removing traits reshuffles chunks
# freely without invalidating anything.

suppressPackageStartupMessages({
  library(data.table)
  library(DBI)
  library(RSQLite)
})

# Open the ledger. create = FALSE returns NULL when the db does not exist yet.
db_open <- function(db_file, create = TRUE) {
  if (!create && !file.exists(db_file)) return(NULL)
  dir.create(dirname(db_file), showWarnings = FALSE, recursive = TRUE)
  dbConnect(SQLite(), db_file)
}

# Separator for flattened composite keys; never part of an id.
KEY_SEP <- "\r"

# Distinct values of key_col already in the ledger. Empty vector if the db/table
# is absent -- so a fresh run naturally has an empty done set.
db_done_keys <- function(db_file, table, key_col = "exposure_id") {
  con <- db_open(db_file, create = FALSE)
  if (is.null(con)) return(character(0))
  on.exit(dbDisconnect(con))
  if (!dbExistsTable(con, table)) return(character(0))
  q <- dbGetQuery(con, sprintf("SELECT DISTINCT %s FROM %s",
                               paste(key_col, collapse = ", "), table))
  if (!nrow(q)) return(character(0))
  # A composite key is flattened the same way callers build it (KEY_SEP): "__" occurs
  # inside exposure ids, so it could join two different tuples into one key.
  do.call(paste, c(lapply(key_col, function(k) as.character(q[[k]])), sep = KEY_SEP))
}

# Keys whose ledger rows were computed before one of their inputs changed. `inputs`
# holds the key columns plus input_mtime (epoch seconds of the unit's newest fine-map
# output). Runners stamp every row with computed_at (the job's start time); a key
# with any row lacking it counts as stale, so an unstamped ledger is recomputed rather
# than trusted. Without this, a changed fine-map left its coloc rows "done" until
# someone purged the ledger by hand.
db_stale_keys <- function(db_file, table, key_col, inputs) {
  con <- db_open(db_file, create = FALSE)
  if (is.null(con)) return(character(0))
  on.exit(dbDisconnect(con))
  if (!dbExistsTable(con, table)) return(character(0))
  kc <- paste(key_col, collapse = ", ")
  q <- if ("computed_at" %in% dbListFields(con, table)) {
    dbGetQuery(con, sprintf(
      "SELECT %s, MIN(computed_at) AS computed_at, SUM(computed_at IS NULL) AS n_null FROM %s GROUP BY %s",
      kc, table, kc))
  } else {
    dbGetQuery(con, sprintf("SELECT DISTINCT %s, NULL AS computed_at, 1 AS n_null FROM %s", kc, table))
  }
  q <- as.data.table(q)
  for (k in key_col) q[, (k) := as.character(get(k))]
  inp <- copy(inputs)
  for (k in key_col) inp[, (k) := as.character(get(k))]
  m <- merge(q, inp, by = key_col)
  m <- m[n_null > 0 | (!is.na(input_mtime) & input_mtime > computed_at)]
  if (!nrow(m)) return(character(0))
  do.call(paste, c(lapply(key_col, function(k) m[[k]]), sep = KEY_SEP))
}

# Upsert keyed on key_col: every existing row for the incoming keys is deleted, then
# the new rows are inserted. Idempotent -- re-importing the same staging file is a no-op
# in effect, which is what makes a failed-then-resubmitted job safe.
db_upsert <- function(db_file, table, dt, key_col = "exposure_id") {
  stopifnot(all(key_col %in% names(dt)))
  if (!nrow(dt)) return(invisible(0L))
  con <- db_open(db_file)
  on.exit(dbDisconnect(con))

  if (!dbExistsTable(con, table)) {
    dbWriteTable(con, table, as.data.frame(dt))
    dbExecute(con, sprintf("CREATE INDEX IF NOT EXISTS idx_%s_key ON %s(%s)",
                           table, table, paste(key_col, collapse = ", ")))
    return(invisible(nrow(dt)))
  }

  # Reconcile schema drift in both directions: new columns are added to the table,
  # columns the table already has but this batch lacks are filled with NA.
  have <- dbListFields(con, table)
  for (nm in setdiff(names(dt), have)) {
    dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN %s", table, nm))
    have <- c(have, nm)
  }
  dt <- copy(dt)
  for (nm in setdiff(have, names(dt))) dt[, (nm) := NA]
  setcolorder(dt, have)

  keys <- unique(as.data.frame(dt)[, key_col, drop = FALSE])
  # Composite keys are matched column-by-column against a temp table rather than on a
  # pasted string, so a value containing the separator can never collide.
  cond <- paste(sprintf("t.%s = k.%s", key_col, key_col), collapse = " AND ")
  dbBegin(con)
  ok <- tryCatch({
    dbWriteTable(con, "stg_keys", keys, overwrite = TRUE, temporary = TRUE)
    dbExecute(con, sprintf(
      "DELETE FROM %s AS t WHERE EXISTS (SELECT 1 FROM stg_keys k WHERE %s)", table, cond))
    dbAppendTable(con, table, as.data.frame(dt))
    TRUE
  }, error = function(e) { message("db_upsert failed: ", conditionMessage(e)); FALSE })
  if (ok) dbCommit(con) else { dbRollback(con); stop("db_upsert rolled back for ", table) }
  invisible(nrow(dt))
}

# Whole ledger as a data.table (used to export the flat result tsv downstream reads).
db_read_all <- function(db_file, table) {
  con <- db_open(db_file, create = FALSE)
  if (is.null(con)) return(NULL)
  on.exit(dbDisconnect(con))
  if (!dbExistsTable(con, table)) return(NULL)
  as.data.table(dbGetQuery(con, sprintf("SELECT * FROM %s", table)))
}

# Append one unit of work to this job's staging file, normalised to a fixed column set
# so the file stays a valid tsv. Appending per unit (not per job) keeps the old
# crash-resume granularity: a job that dies mid-chunk still has its finished units.
stage_append <- function(stage_file, dt, cols) {
  dt <- copy(dt)
  for (nm in setdiff(cols, names(dt))) dt[, (nm) := NA]
  dt <- dt[, ..cols]
  fwrite(dt, stage_file, sep = "\t", append = file.exists(stage_file),
         col.names = !file.exists(stage_file))
}

# Staging file for the current job. SLURM_JOB_ID is unique per submission, so chunk
# boundaries changing between submissions can never collide.
stage_file_path <- function(stage_dir, fallback = "local") {
  jid <- Sys.getenv("SLURM_JOB_ID", "")
  if (!nzchar(jid)) jid <- fallback
  dir.create(stage_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(stage_dir, paste0("job_", jid, ".tsv"))
  # Afresh at job start: a requeued job keeps its SLURM_JOB_ID and reruns its whole
  # slice, so appending would stage its finished units twice.
  if (file.exists(f) && !file.remove(f)) stop("could not reset ", f)
  f
}

source(path.expand("code/lib/common.R"))   # msg, assert_clump_dirs, GWAS_BASE, CLUMP_SUBDIR
