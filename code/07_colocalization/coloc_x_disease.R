#!/usr/bin/env Rscript
# Shared body of the two "x disease" coloc runners (protein x disease,
# MetaBrain eQTL x disease). They were 77% identical: same loop, same
# coloc.susie call, same status vocabulary, differing only in where the exposure
# side comes from and which id columns travel with a row.
#
# The eQTL chains (coloc_protein_x_qtl.R) keep their own loop: they iterate proteins
# rather than pairs, prefilter candidates from a (tissue, ENSG) manifest, and append
# to staging/ rather than writing one chunk per row range. What is kept from each
# comparison is shared (coloc_acc_* in functions_coloc.R).
#
# coloc.susie only: both sides must be fine-mapped (>= 1 credible set). Per
# exposure locus x overlapping disease locus, run coloc.susie on the two saved
# SuSiE fits, aligned by rsid (see functions_coloc.R). There is no coloc.abf
# fallback -- where the disease side has no credible set in the region there is
# nothing to colocalize (PP4 would be ~0), so a status row is recorded rather
# than mixing in a single-causal-variant method. Rows with PP4 >= PP4_KEEP are kept
# as status "ok"; a pair tested only below it gets one "no_sig_coloc" row carrying
# its highest-PP4 comparison, so every tested pair has a PP4. "ok" still means
# retained at PP4 >= PP4_KEEP.
#
# cfg:
#   coloc_dir       output/input root under ~/db/shared_data/coloc
#   pair_ok         function(pair) -> filtered pair table, exposure_id derived
#   exposure_loci   function(p) -> finemap_loci for this pair's exposure
#   exposure_rds    function(p, lead) -> path to the exposure's saved SuSiE fit
#   id_cols         function(p) -> named list of id columns carried on every row
#   progress_every  rows between progress lines
run_coloc_x_disease <- function(cfg) {
  # Every row is stamped with the job's start; 01_pair redoes a pair whose fine-maps
  # are newer (db_stale_keys in functions_db.R).
  RUN_STARTED <- as.integer(Sys.time())
  COLOC_DIR <- path.expand(cfg$coloc_dir)
  # Overridable so a run can be compared against another build without touching
  # the real chunk directory.
  TMP_DIR   <- Sys.getenv("COLOC_TMP_DIR", file.path(COLOC_DIR, "tmp_chunks"))
  # COLOC_TODO: this submission's own copy (submit_coloc_chain.sh).
  PAIR_FILE <- Sys.getenv("COLOC_TODO", file.path(COLOC_DIR, "info", "coloc_pair_todo.csv"))
  dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

  ROW_START <- as.integer(Sys.getenv("COLOC_PAIR_ROW_START", "1"))
  ROW_END   <- as.integer(Sys.getenv("COLOC_PAIR_ROW_END",   "0"))
  if (ROW_END < ROW_START) ROW_END <- ROW_START

  t0  <- proc.time()[[3L]]
  msg <- function(...) message(sprintf("[%.1fs] ", proc.time()[[3L]] - t0), ...)

  todo <- fread(PAIR_FILE)
  # A slice past the end used to read NA rows and write a chunk of NA keys.
  if (ROW_START > nrow(todo)) stop("rows ", ROW_START, "-", ROW_END, " are past the ", nrow(todo), " rows of ", PAIR_FILE)
  pair <- cfg$pair_ok(todo[ROW_START:min(ROW_END, nrow(todo))])
  if (!"o_leads" %in% names(pair)) stop(PAIR_FILE, " has no o_leads column; rerun the 01 pair script")
  msg("rows ", ROW_START, "-", ROW_END, " -> ", nrow(pair), " pairs")

  # Disease loci for the most recent outcome; pairs are sorted by outcome_id, so
  # a one-deep cache is enough.
  oloc_id <- ""; oloc_dat <- NULL
  get_oloc <- function(id, src) {
    if (id == oloc_id) return(oloc_dat)
    oloc_dat <<- finemap_loci(disease_finemap_path(src, id), cis_only = FALSE)
    oloc_id <<- id
    oloc_dat
  }

  stub <- function(p, status, cs_id = NA_integer_)
    do.call(data.table, c(cfg$id_cols(p),
                          list(outcome_id = p$outcome_id, cs_id = cs_id, status = status)))

  results <- list()
  n_keep  <- 0L
  for (i in seq_len(nrow(pair))) {
    p <- pair[i]
    regs <- cfg$exposure_loci(p)
    if (is.null(regs) || !nrow(regs)) {
      results[[length(results) + 1L]] <- stub(p, "no_region"); next
    }
    oloc <- get_oloc(p$outcome_id, p$outcome_source)
    acc  <- coloc_acc()

    for (k in seq_len(nrow(regs))) {
      r     <- regs[k]
      rds_e <- cfg$exposure_rds(p, r$clump_lead_snp)
      # The locus has credible-set rows, so its SuSiE fit must exist; skipping it
      # would record the pair as done with no coloc.
      if (!file.exists(rds_e)) stop("finemapped locus without its SuSiE fit: ", rds_e)
      cand <- loci_near(oloc, r)
      if (is.null(cand) || !nrow(cand)) next

      # Exposure fit is constant across the candidate loop -- load once per region.
      fe_k <- rekey_susie(readRDS(rds_e))
      if (is.null(fe_k)) next

      for (j in seq_len(nrow(cand))) {
        rds_o <- disease_susie_path(p$outcome_source, p$outcome_id, cand$clump_lead_snp[j])
        sres  <- run_coloc_susie(rds_e, rds_o, fe = fe_k)
        annot <- function(d, status) {
          d[, names(cfg$id_cols(p)) := cfg$id_cols(p)]
          d[, `:=`(outcome_id = p$outcome_id, cs_id = k, status = status,
                   coloc_method = "susie",
                   lead_snp = r$clump_lead_snp, o_lead = cand$clump_lead_snp[j],
                   chr = r$chr, pos = r$pos)]
        }
        n_keep <- n_keep + coloc_acc_add(acc, sres, annot)
      }
    }

    results[[length(results) + 1L]] <-
      coloc_acc_rows(acc, disease_no_cs_status(p$o_leads, oloc), function(st) stub(p, st),
                     unit = paste(unlist(p[, intersect(c("exposure_id", "tissue", "ENSG", "outcome_id"), names(p)), with = FALSE]),
                                  collapse = " x "))
    if (i %% cfg$progress_every == 0L)
      msg("done ", i, "/", nrow(pair), " (PP4>", PP4_KEEP, ": ", n_keep, ")")
  }

  dt <- rbindlist(results, fill = TRUE)
  # Named by job, not row range: row ranges of two submissions can coincide.
  jid <- Sys.getenv("SLURM_JOB_ID", "")
  chunk_file <- file.path(TMP_DIR, if (nzchar(jid)) sprintf("chunk_job%s.tsv", jid) else
                            sprintf("chunk_%07d_%07d.tsv", ROW_START, ROW_END))
  dt[, computed_at := RUN_STARTED]
  fwrite(dt, paste0(chunk_file, ".part"), sep = "\t")
  if (!file.rename(paste0(chunk_file, ".part"), chunk_file)) stop("could not write ", chunk_file)
  msg("wrote ", nrow(dt), " rows -> ", chunk_file, " (PP4>", PP4_KEEP, " RDS: ", n_keep, ")")
  coloc_failed_report()
  invisible(chunk_file)
}
