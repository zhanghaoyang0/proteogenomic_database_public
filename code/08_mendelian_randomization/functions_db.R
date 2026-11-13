# SQLite "done" ledger for the MR arms.
#
# Why not mirror the whole result table (as coloc/functions_db.R does): the MR result
# is 9 slim columns over ~10.8M pairs -- 390 MB gzipped, ~1.5 GB raw. Putting that in
# SQLite would cost ~1 GB plus index, i.e. more than the gzip it replaces. The result
# table stays a gzipped tsv; only "which pairs have run" moves into the database.
#
# Storing that as one row per pair would itself be ~1 GB of TEXT keys. But the pair set
# is (near-)rectangular -- 19,479 exposures x 555 outcomes -- so the ledger keeps one
# row per EXPOSURE holding a bitmap over outcomes:
#
#   outcome_dict(outcome_id, idx)   append-only; an idx is never reused or renumbered,
#                                   so bitmaps written by older runs stay valid
#   done_expo(exposure_id, bits)    bits[idx] = 1 means that pair has a result row
#
# ~19.5k rows / a few MB, and it is exact -- no "count matched so assume done"
# shortcut that could silently skip a swapped-in outcome.

suppressPackageStartupMessages({
  library(data.table)
  library(DBI)
  library(RSQLite)
})

db_open <- function(db_file, create = TRUE) {
  if (!create && !file.exists(db_file)) return(NULL)
  dir.create(dirname(db_file), showWarnings = FALSE, recursive = TRUE)
  con <- dbConnect(SQLite(), db_file)
  # DELETE, not MEMORY: a collect killed mid-transaction then rolls back instead of
  # corrupting mr.db; like MEMORY it writes no WAL side files on GPFS.
  dbExecute(con, "PRAGMA journal_mode = DELETE")
  dbExecute(con, "CREATE TABLE IF NOT EXISTS outcome_dict (
                    outcome_id TEXT PRIMARY KEY, idx INTEGER NOT NULL)")
  dbExecute(con, "CREATE TABLE IF NOT EXISTS done_expo (
                    exposure_id TEXT PRIMARY KEY, bits BLOB NOT NULL)")
  con
}

# Assign an idx to any outcome_id not yet in the dictionary. Existing ids keep theirs.
dict_sync <- function(con, outcome_ids) {
  d <- as.data.table(dbGetQuery(con, "SELECT outcome_id, idx FROM outcome_dict"))
  new <- setdiff(unique(as.character(outcome_ids)), d$outcome_id)
  if (length(new)) {
    start <- if (nrow(d)) max(d$idx) else 0L
    add <- data.table(outcome_id = new, idx = start + seq_along(new))
    dbAppendTable(con, "outcome_dict", as.data.frame(add))
    d <- rbind(d, add)
  }
  setkey(d, outcome_id)
  d
}

dict_read <- function(con) {
  d <- as.data.table(dbGetQuery(con, "SELECT outcome_id, idx FROM outcome_dict"))
  if (nrow(d)) setkey(d, outcome_id)
  d
}

# Trailing bits an older (shorter) bitmap never had are read as FALSE, not NA.
bits_unpack <- function(raw, n) {
  v <- as.logical(rawToBits(raw))
  length(v) <- n
  !is.na(v) & v
}

# done_lookup: TRUE/FALSE per row of `pair` (needs exposure_id, outcome_id).
# Builds an n_exposure x n_outcome logical matrix -- 10.8M cells is ~43 MB, versus
# ~1.5 GB for the pasted-key vector the old scan-the-result-file approach allocated.
done_lookup <- function(con, pair) {
  dict <- dict_read(con)
  if (!nrow(dict)) return(rep(FALSE, nrow(pair)))
  n_out <- max(dict$idx)

  expo <- unique(as.character(pair$exposure_id))
  rows <- as.data.table(dbGetQuery(con, "SELECT exposure_id, bits FROM done_expo"))
  rows <- rows[exposure_id %in% expo]
  if (!nrow(rows)) return(rep(FALSE, nrow(pair)))

  m <- matrix(FALSE, nrow = length(expo), ncol = n_out,
              dimnames = list(expo, NULL))
  ei <- match(rows$exposure_id, expo)
  for (i in seq_len(nrow(rows))) m[ei[i], ] <- bits_unpack(rows$bits[[i]], n_out)

  oi <- dict$idx[match(as.character(pair$outcome_id), dict$outcome_id)]
  ri <- match(as.character(pair$exposure_id), expo)
  ok <- !is.na(oi) & !is.na(ri)
  out <- rep(FALSE, nrow(pair))
  out[ok] <- m[cbind(ri[ok], oi[ok])]
  out
}

# done_record: OR the (exposure_id, outcome_id) pairs in `dt` into the ledger.
# replace = TRUE makes the ledger exactly `dt`: 04 passes the whole current result, so
# pairs dropped from mr_pair.csv.gz stop counting as done (they would never be rerun).
done_record <- function(con, dt, replace = FALSE) {
  dt <- unique(as.data.table(dt)[, .(exposure_id = as.character(exposure_id),
                                     outcome_id  = as.character(outcome_id))])
  if (!nrow(dt)) return(invisible(0L))
  dict <- dict_sync(con, dt$outcome_id)
  n_out <- max(dict$idx)
  dt[, idx := dict$idx[match(outcome_id, dict$outcome_id)]]

  old <- if (replace) data.table(exposure_id = character(), bits = list()) else
    as.data.table(dbGetQuery(con, "SELECT exposure_id, bits FROM done_expo"))
  old <- old[exposure_id %in% unique(dt$exposure_id)]

  agg <- dt[, .(idx = list(idx)), by = exposure_id]
  # Row in `old` carrying this exposure's previous bitmap, or NA when it is new.
  pi <- if (nrow(old)) match(agg$exposure_id, old$exposure_id) else rep(NA_integer_, nrow(agg))
  # Kept out of the data.table: a one-row list column gets simplified back to a bare
  # raw vector, which RSQLite then refuses to bind.
  bl <- lapply(seq_len(nrow(agg)), function(i) {
    v <- logical(ceiling(n_out / 8) * 8)
    if (!is.na(pi[i])) v[seq_len(n_out)] <- bits_unpack(old$bits[[pi[i]]], n_out)
    v[agg$idx[[i]]] <- TRUE
    packBits(v, type = "raw")
  })

  dbBegin(con)
  ok <- tryCatch({
    dbWriteTable(con, "stg_expo", data.frame(exposure_id = agg$exposure_id),
                 overwrite = TRUE, temporary = TRUE)
    if (replace) dbExecute(con, "DELETE FROM done_expo")
    else dbExecute(con, "DELETE FROM done_expo WHERE exposure_id IN (SELECT exposure_id FROM stg_expo)")
    df <- data.frame(exposure_id = agg$exposure_id, stringsAsFactors = FALSE)
    df$bits <- bl
    dbAppendTable(con, "done_expo", df)
    TRUE
  }, error = function(e) { message("done_record failed: ", conditionMessage(e)); FALSE })
  if (ok) dbCommit(con) else { dbRollback(con); stop("done_record rolled back") }
  invisible(nrow(agg))
}

# fwrite(compress = "gzip") leaves the gzip stream unterminated when the table has no
# data rows: fread still reads it, but zcat/gunzip -t report "unexpected end of file",
# which aborts the submit scripts under `set -o pipefail`. Write the header through a
# gz connection in that case so the file is always a valid archive.
fwrite_gz <- function(dt, path) {
  if (nrow(dt)) {
    fwrite(dt, path, compress = "gzip")
  } else {
    con <- gzfile(path, "wb")
    on.exit(close(con))
    writeLines(paste(names(dt), collapse = ","), con)
  }
  invisible(path)
}

db_done_n <- function(con) {
  as.integer(dbGetQuery(con, "SELECT COUNT(*) n FROM done_expo")$n)
}

# A clumped file with >= 2 lines has at least one lead SNP (line 1 is the header).
# Used by 01_pair_protein2disease.R.
clump_has_lead <- function(dir) {
  files <- list.files(dir, pattern = "\\.clumped$", full.names = TRUE)
  ids   <- sub("\\.clumped$", "", basename(files))
  keep <- vapply(files, function(f) length(readLines(f, n = 2)) >= 2, logical(1L))
  ids[keep]
}

# Write the not-yet-done pair list, and return its path. Shared by
# 01_pair_protein2disease.R.
#
# Two reasons done comes from the ledger rather than from files:
#  - tmp_chunks/ alone is not enough: 04_collect deletes it on success, so a finished
#    run would report its entire pair list as todo and a resubmit would recompute all
#    ~10.8M pairs. The ledger survives that cleanup.
#  - the result tables are not enough either: deriving done from them meant pasting
#    ~1.5M string keys. The ledger answers the same question from an indexed table;
#    the result files are re-read only once, by 04_collect, to seed an empty ledger.
# Chunks are still unioned in, to cover a run that finished but was not collected yet.
#
# `pair` gets a temporary `done` column, dropped again before returning (data.table
# updates by reference, so the caller's `pair` is left with its original columns).
write_pair_todo <- function(pair, db_file, tmp_dir, info_dir) {
  con <- db_open(db_file)
  pair[, done := done_lookup(con, pair)]
  msg("ledger done: ", sum(pair$done), " (", db_done_n(con), " exposures on file)")

  chunk_files <- list.files(tmp_dir, pattern = "^chunk_.*\\.tsv$", full.names = TRUE)
  if (length(chunk_files)) {
    msg("scanning ", length(chunk_files), " uncollected chunk files ...")
    ck <- rbindlist(lapply(chunk_files, fread, select = c("exposure_id", "outcome_id")), fill = TRUE)
    pair[ck, on = .(exposure_id, outcome_id), done := TRUE]
  }
  dbDisconnect(con)

  todo <- pair[done == FALSE]
  todo[, done := NULL]
  pair[, done := NULL]
  out_todo <- file.path(info_dir, "mr_pair_todo.csv.gz")
  fwrite_gz(todo, out_todo)
  msg("todo: ", nrow(todo), "/", nrow(pair), " -> ", out_todo)
  invisible(out_todo)
}

source(path.expand("code/lib/common.R"))   # msg, assert_clump_dirs, GWAS_BASE, CLUMP_SUBDIR
