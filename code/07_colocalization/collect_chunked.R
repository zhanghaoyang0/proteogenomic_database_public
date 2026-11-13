#!/usr/bin/env Rscript
# Shared body of the chunk-based collectors (protein x disease, mb x disease).
# Those two differ only in a directory, a table name and the key columns; the
# other ~75 lines were duplicated, so a fix to one had to be remembered for the
# other. protein_x_eqtl_03_collect.R is NOT one of these -- the eQTL chains
# stage into staging/ + data/ rather than tmp_chunks/, so they keep their own
# collector.
#
# Source: tmp_chunks/chunk_*.tsv, this run's output, one file per job (removed
# after import). The one-time seeding from a pre-ledger result/*.tsv was removed
# on 2026-09-30: on an empty ledger it re-imported whatever old table was there.
#
# Import is an upsert on `key_cols`, so re-running a subset of pairs replaces
# exactly their rows. Chunks are deleted only AFTER they are in the ledger: the
# code this replaced deleted them while they were still the sole record.
#
# annotate: optional function(all, pair) run on the exported rows just before the
# write -- the chain's add_component_flags() call, which needs the pair list for
# outcome_source.
collect_chunked <- function(coloc_dir, db_table, key_cols, annotate = NULL) {
  COLOC_DIR <- path.expand(coloc_dir)
  TMP_DIR   <- file.path(COLOC_DIR, "tmp_chunks")
  DB_FILE   <- file.path(COLOC_DIR, "coloc.db")
  PAIR_FILE <- file.path(COLOC_DIR, "info", "coloc_pair.csv")
  RESULT    <- file.path(COLOC_DIR, "result", paste0(db_table, ".tsv"))
  KEEP_CHUNKS <- nzchar(Sys.getenv("KEEP_CHUNKS", ""))
  dir.create(dirname(RESULT), showWarnings = FALSE, recursive = TRUE)

  pair_key <- function(d)
    do.call(paste, c(lapply(key_cols, function(k) as.character(d[[k]])), sep = KEY_SEP))   # functions_db.R

  # ---- 2. import this run's chunks --------------------------------------------
  chunks <- list.files(TMP_DIR, pattern = "^chunk_.*\\.tsv$", full.names = TRUE)
  if (length(chunks)) {
    # A pair has many legitimate rows (one per disease locus x CS pair), so dedup
    # by pair: if a pair was written by more than one chunk (a resubmit), keep
    # only the newest chunk's rows. Deduping on key + cs_id dropped real rows.
    chunks <- chunks[order(file.mtime(chunks))]
    all <- rbindlist(lapply(seq_along(chunks), function(i) fread(chunks[i])[, .src := i]), fill = TRUE)
    msg("chunks: ", length(chunks), " file(s), ", nrow(all), " rows")
    all[, .pk := pair_key(all)]
    all <- all[all[, .I[.src == max(.src)], by = .pk]$V1]
    all[, c(".src", ".pk") := NULL]
    msg("after keeping the newest chunk per pair: ", nrow(all))
    db_upsert(DB_FILE, db_table, all, key_cols)
    msg("imported chunks -> ", DB_FILE)
  } else {
    msg("no chunk files in ", TMP_DIR)
  }

  # ---- 3. export the current pair universe, only when it is complete -----------
  # The ledger is add-only, so it can hold pairs that are no longer in the pair list
  # (a changed outcome set); exporting them would inflate coverage and the results.
  # Coverage is checked before the table is written, so a partial run never replaces
  # the result downstream reads. The chunks are already in the ledger either way.
  if (!file.exists(PAIR_FILE)) stop("missing pair list: ", PAIR_FILE)
  pair  <- fread(PAIR_FILE)
  pairs <- unique(pair_key(pair))
  all <- db_read_all(DB_FILE, db_table)
  if (is.null(all) || !nrow(all)) stop("ledger is empty: ", DB_FILE)
  all[, .pk := pair_key(all)]
  n_stale <- uniqueN(all[!.pk %in% pairs, .pk])
  if (n_stale) msg("ledger pairs not in the current pair list (not exported): ", n_stale)
  all <- all[.pk %in% pairs]
  got_n <- uniqueN(all$.pk)
  all[, .pk := NULL]
  msg("ledger: ", nrow(all), " rows | coverage: ", round(100 * got_n / length(pairs), 2),
      "% (", got_n, "/", length(pairs), ")")
  if (got_n < length(pairs) && !nzchar(Sys.getenv("ALLOW_PARTIAL", "")))
    stop(length(pairs) - got_n, " pair(s) not computed -- run them, or set ALLOW_PARTIAL=1")

  cat("\nstatus breakdown:\n"); print(all[, .N, by = status][order(-N)])

  # computed_at is ledger bookkeeping (01_pair staleness), not a result column.
  if (!is.null(annotate)) annotate(all, pair)   # before the drop: it checks computed_at
  if ("computed_at" %in% names(all)) all[, computed_at := NULL]
  fwrite(all, RESULT, sep = "\t")
  msg("wrote -> ", RESULT)

  # ---- 4. drop the chunks (safe now: they are in the ledger) --------------------
  if (length(chunks) && !KEEP_CHUNKS) {
    unlink(chunks)
    msg("removed ", length(chunks), " chunk file(s) (set KEEP_CHUNKS=1 to keep)")
  }
  invisible(RESULT)
}
