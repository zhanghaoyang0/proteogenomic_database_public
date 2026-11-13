# Shared helpers for coloc runs (protein x disease, protein x mb, etc.).
# Window-based: loci are compared when their credible-set spans lie within CIS_HALF_BP
# (loci_near).

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

MIN_N_OVERLAP <- 50L
PP4_KEEP      <- 0.7  # PP4 keep threshold (between moderate 0.5 and strong 0.8)

source(path.expand("code/lib/common.R"))   # GWAS_BASE

finemap_path <- function(src, id) {   # vectorised over src and id
  sub <- ifelse(src == "BioFINDER1", "bf1/finemapping/protein", "bf2/finemapping/protein")
  file.path(GWAS_BASE, sub, paste0(id, ".finemapped.tsv.gz"))
}

# ---- coloc.susie (both sides fine-mapped) -------------------------------------
# Cross-cohort SNP matching is by rsid. Fine-mapping (01_finemapping.R) keeps only
# variants present in the 1kg ref (EUR.bim is 100% bare rsid), so every variant in a
# saved susie fit carries a 1kg rsid — empirically fit$snp is 100% bare rsid in all
# three cohorts (protein/FinnGen/MetaBrain). rsid matching is therefore equivalent to
# CHR:POS matching on the fits, with no clean-file position lookup needed. coloc.susie
# uses lbf magnitude, so no allele flip is needed.

# Saved-susie path per cohort: <cohort>/finemapping/.../susie_rds/<id>__<clump_lead>.susie.rds
# Cohort tag from the exposure_id prefix. The saved-susie path is a pure function of
# result-table columns, so a per-pair index of these paths never needs storing.
protein_src <- function(exposure_id) {
  ifelse(startsWith(exposure_id, "bf1__"), "BioFINDER1", "BioFINDER2")
}
protein_susie_path <- function(src, id, lead) {   # vectorised
  sub <- ifelse(src == "BioFINDER1", "bf1/finemapping/protein/susie_rds",
                "bf2/finemapping/protein/susie_rds")
  file.path(GWAS_BASE, sub, paste0(id, "__", lead, ".susie.rds"))
}
disease_finemap_path <- function(src, id) {   # vectorised over src and id
  sub <- ifelse(src == "Other", "other/finemapping", "finngen/finemapping")
  file.path(GWAS_BASE, sub, paste0(id, ".finemapped.tsv.gz"))
}
disease_susie_path <- function(src, id, lead) {   # vectorised
  sub <- ifelse(src == "Other", "other/finemapping/susie_rds", "finngen/finemapping/susie_rds")
  file.path(GWAS_BASE, sub, paste0(id, "__", lead, ".susie.rds"))
}
# MetaBrain eQTL side: one finemap locus per (tissue, ENSG); same susie/finemap layout.
metabrain_susie_path <- function(tissue, ensg, lead) {
  file.path(GWAS_BASE, "metabrain/finemapping/susie_rds",
            paste0(tissue, "_", ensg, "__", lead, ".susie.rds"))
}
metabrain_finemap_path <- function(tissue, ensg) {
  file.path(GWAS_BASE, "metabrain/finemapping",
            paste0(tissue, "_", ensg, ".finemapped.tsv.gz"))
}
# Every molecular-QTL resource uses the same finemapping layout -- one locus per
# (tissue, ENSG), files keyed "<tissue>_<ENSG>" -- so the three differ only by
# directory.
eqtl_susie_path <- function(resource, tissue, ensg, lead) {
  file.path(GWAS_BASE, resource, "finemapping/susie_rds",
            paste0(tissue, "_", ensg, "__", lead, ".susie.rds"))
}
eqtl_finemap_path <- function(resource, tissue, ensg) {
  file.path(GWAS_BASE, resource, "finemapping",
            paste0(tissue, "_", ensg, ".finemapped.tsv.gz"))
}
# Per-locus finemap rows (one per clump lead): chr/pos of the top-PIP variant (the
# anchor every result row reports) and lo/hi, the span of all the locus's credible-set
# variants. cis_only restricts to cis CS (protein side).
finemap_loci <- function(finemap_file, cis_only = FALSE) {
  if (!file.exists(finemap_file)) return(NULL)
  d <- fread(finemap_file)
  if (!nrow(d)) return(NULL)
  if (cis_only) { d <- d[is_cis %in% TRUE]; if (!nrow(d)) return(NULL) }
  setorder(d, -pip)
  d[, .(chr = chr[1L], pos = pos[1L], lo = min(pos), hi = max(pos)), by = clump_lead_snp]
}

# The candidate rule every chain uses: loci whose credible-set span lies within
# CIS_HALF_BP of locus r's span. Until 2026-10-25 the two anchors had to be within
# CIS_HALF_BP, which skipped loci whose credible sets overlap while their top-PIP
# variants are far apart (151 such protein x disease combos, 8 joint MR + coloc pairs
# missed; coloc/docs/issues.md history).
loci_near <- function(loci, r) {
  if (is.null(loci) || !nrow(loci)) return(loci)
  loci[as.character(chr) == as.character(r$chr) &
       lo <= r$hi + CIS_HALF_BP & hi >= r$lo - CIS_HALF_BP]
}

# Attach rsid colnames to a susie fit's lbf_variable so coloc.susie can match by SNP.
# Duplicate rsids (rare) get unique dummy names so they never cross-match. Keeps full
# width (sets$cs indices stay valid). Returns NULL if the fit has no credible set.
rekey_susie <- function(fit) {
  if (is.null(fit) || is.null(fit$sets$cs) || !length(fit$sets$cs)) return(NULL)
  k   <- fit$snp
  dup <- duplicated(k); k[dup] <- paste0(k[dup], "_d", which(dup))
  colnames(fit$lbf_variable) <- k
  fit$snp <- k
  fit
}

# coloc.susie on two saved fits aligned by rsid. Returns a data.table with one row
# per (CS_e, CS_o) combination; NULL when either side has no CS; or a status string
# when the test could not run: "low_overlap" (< MIN_N_OVERLAP shared SNPs) or
# "coloc_error" (coloc.susie errored). These used to be NULL too, so a failed test
# was recorded as "no CS". There is no coloc.abf fallback.
#
# A missing fit is an error, not a NULL: callers only pass loci listed in a
# finemapped file, which lists only loci with a credible set, so the fit must exist.
# Returning NULL here used to record the pair as done with no coloc.
#
# fe: the already-loaded, already-rekeyed exposure fit. The exposure side is constant
# across a caller's inner candidate loop, so callers load it once per region and pass
# it here; leaving fe NULL keeps the old self-contained behaviour.
run_coloc_susie <- function(rds_e, rds_o, fe = NULL) {
  if (!file.exists(rds_o)) stop("finemapped locus without its SuSiE fit: ", rds_o)
  if (is.null(fe)) {
    if (!file.exists(rds_e)) stop("finemapped locus without its SuSiE fit: ", rds_e)
    fe <- rekey_susie(readRDS(rds_e))
  }
  fo <- rekey_susie(readRDS(rds_o))
  if (is.null(fe) || is.null(fo)) return(NULL)
  if (length(intersect(fe$snp, fo$snp)) < MIN_N_OVERLAP) return("low_overlap")
  res <- tryCatch(coloc.susie(fe, fo), error = function(e) {
    message("coloc.susie failed (", basename(rds_o), "): ", conditionMessage(e))
    NULL
  })
  if (is.null(res)) return("coloc_error")
  if (is.null(res$summary) || !nrow(res$summary)) return(NULL)
  s <- as.data.table(res$summary)
  # idx1/idx2 say which credible set of each fit this row compared, as the SuSiE
  # component number -- look them up through the fit's sets$cs_index, not by
  # subscripting sets$cs. Carrying them lets a downstream consumer (the locus plot)
  # address the CS directly instead of re-deriving it by searching for the hit
  # variant, which breaks at a locus with several CS if the hit falls outside them.
  s[, .(n_snps = nsnps, e_hit = as.character(hit1), o_hit = as.character(hit2),
        e_idx = as.integer(idx1), o_idx = as.integer(idx2),
        PP0 = PP.H0.abf, PP1 = PP.H1.abf, PP2 = PP.H2.abf,
        PP3 = PP.H3.abf, PP4 = PP.H4.abf)]
}

# ---- which compared components the fine-mapping keeps --------------------------
# coloc.susie compares every component of both fits, but a .finemapped.tsv.gz keeps
# only the GWS credible sets (2026-10-24 audit: 29% of protein x disease ok rows used
# a component absent on at least one side). A component is "retained" when its
# credible set, as a set of rsids, is a credible set of the table under ANY clump
# window: the cross-window dedup keeps an identical set under a neighbouring lead,
# which is the same signal. cis_only (protein exposures) also requires is_cis, the
# way finemap_loci(cis_only = TRUE) anchors the run. "converged" is the fit's own
# fit$converged. primary = ok with both sides retained and converged; status-only
# and pre-idx rows get NA flags and primary FALSE.
cs_key <- function(snps) paste(sort(unique(snps)), collapse = ";")

# 01_finemapping.R names the fit of a clump lead without an rsid "<id>__NA.susie.rds",
# but the result tables carry that lead as "" (7 protein_x_gtex rows, 2026-10-24).
lead_or_na <- function(x) replace(x, !is.na(x) & x == "", NA)

# computed_at (ledger seconds) guards against a fit refit after the row was computed:
# its idx would then point at a different, or no, component. Stops rather than flag
# the wrong component -- the chain's pair step re-queues such rows; rerun it first.
# The clock is the finemapped table's mtime, the one the pair step re-queues on: 01
# renames the table last, so it is never older than its fits. (The fit's own mtime
# here could stop the collect for a row the pair step would never re-queue.)
component_flags <- function(rds, fm, idx, computed_at, cis_only = FALSE) {
  out <- data.table(retained = rep(NA, length(idx)), converged = rep(NA, length(idx)))
  has <- !is.na(idx)
  key_cache <- new.env()   # one table serves every fit (clump window) of its trait
  table_keys <- function(path) {
    if (is.null(key_cache[[path]])) {
      if (!file.exists(path)) stop("result row without its finemapped table: ", path)
      d <- fread(path, select = c("clump_lead_snp", "cs_id", "is_cis", "snp"))
      if (cis_only) d <- d[is_cis %in% TRUE]
      key_cache[[path]] <- d[, cs_key(snp), by = .(clump_lead_snp, cs_id)]$V1
    }
    key_cache[[path]]
  }
  fits <- unique(rds[has])
  miss <- fits[!file.exists(fits)]
  if (length(miss)) stop(length(miss), " result fit(s) missing, e.g. ", miss[1L])
  fm_mt <- as.numeric(file.mtime(fm))
  stale <- has & (fm_mt > as.numeric(computed_at)) %in% TRUE
  if (any(stale))
    stop(sum(stale), " row(s) on ", uniqueN(rds[stale]), " fit(s) refit after the row was ",
         "computed, e.g. ", rds[stale][1L], " -- rerun the chain (pair -> run) before collecting")
  for (f in fits) {
    i <- which(has & rds == f)
    keys <- table_keys(fm[i[1L]])
    fit <- readRDS(f)
    k <- match(idx[i], fit$sets$cs_index)
    if (anyNA(k)) stop("component not a credible set of its fit: ", f)
    set(out, i, "retained", vapply(k, function(j) cs_key(fit$snp[fit$sets$cs[[j]]]) %chin% keys, NA))
    set(out, i, "converged", isTRUE(fit$converged))
  }
  out
}

# Adds e_retained, o_retained, e_converged, o_converged and primary to a coloc result
# table. Paths are per row, built by the chain's collector.
add_component_flags <- function(dt, e_rds, e_fm, o_rds, o_fm, e_cis_only) {
  if (!"computed_at" %in% names(dt)) stop("add_component_flags needs the ledger's computed_at")
  e <- component_flags(e_rds, e_fm, dt$e_idx, dt$computed_at, cis_only = e_cis_only)
  o <- component_flags(o_rds, o_fm, dt$o_idx, dt$computed_at)
  dt[, `:=`(e_retained = e$retained, o_retained = o$retained,
            e_converged = e$converged, o_converged = o$converged)]
  dt[, primary := status == "ok" & e_retained %in% TRUE & o_retained %in% TRUE &
                    e_converged %in% TRUE & o_converged %in% TRUE]
  if (dt[status == "ok", anyNA(e_retained) || anyNA(o_retained)])
    stop("ok rows without component flags (missing e_idx/o_idx)")
  msg("component flags: ok ", sum(dt$status == "ok"), " -> primary ", sum(dt$primary),
      " (e not retained ", dt[status == "ok", sum(!e_retained)],
      ", o not retained ", dt[status == "ok", sum(!o_retained)],
      ", non-converged fit ", dt[status == "ok", sum(!(e_converged & o_converged))], ")")
  invisible(dt)
}

# Status of a pair that kept no PP4 row: "no_sig_coloc" if any locus pair was tested,
# else why none could be ("coloc_error" before "low_overlap"), else no CS on a side.
# Downstream counts only "ok" and "no_sig_coloc" as tested, so the denominators do
# not change; the label just stops calling a failed test "no CS".
pair_status <- function(ran_any, fails, no_cs_status) {
  if (ran_any) "no_sig_coloc"
  else if ("coloc_error" %in% fails) "coloc_error"
  else if ("low_overlap" %in% fails) "low_overlap"
  else no_cs_status
}

# What one unit keeps from its coloc.susie comparisons -- the rule both runners share
# (coloc_x_disease.R per pair, coloc_protein_x_qtl.R per protein). Until 2026-10-24
# each runner had its own copy; the 2026-10-23 best-PP4 change reached only one.
#   acc <- coloc_acc()                      one per unit
#   coloc_acc_add(acc, sres, annot)         after each run_coloc_susie();
#                                           annot(d, status) stamps the unit's columns
#   coloc_acc_rows(acc, no_cs_status, stub) the unit's rows
# Rows with PP4 >= PP4_KEEP are "ok". A unit tested only below it gets one
# "no_sig_coloc" row carrying its highest-PP4 comparison; otherwise stub(status).
coloc_acc <- function() {
  acc <- new.env()
  acc$rows <- list(); acc$best <- NULL; acc$ran_any <- FALSE; acc$fails <- character()
  acc
}
coloc_acc_add <- function(acc, sres, annot) {
  if (is.character(sres)) { acc$fails <- c(acc$fails, sres); return(0L) }
  if (is.null(sres)) return(0L)
  acc$ran_any <- TRUE
  top <- sres[!is.na(PP4)]; top <- top[which.max(top$PP4)]
  if (nrow(top) && (is.null(acc$best) || top$PP4 > acc$best$PP4))
    acc$best <- annot(top, "no_sig_coloc")
  sig <- sres[!is.na(PP4) & PP4 >= PP4_KEEP]
  if (!nrow(sig)) return(0L)
  acc$rows[[length(acc$rows) + 1L]] <- annot(sig, "ok")
  1L
}
# A coloc.susie error in a unit that also ran other comparisons is not in its rows (the
# unit is ok or no_sig_coloc), so it is reported here: an errored test, perhaps on the
# unit's strongest locus, must not pass for a tested negative unseen. Only coloc_error
# counts; low_overlap (too few shared variants) is a property of the data, not an error.
# None in the 2026-10 results, so the result tables carry no extra column for it.
COLOC_ERRORS <- new.env(); COLOC_ERRORS$units <- character(0); COLOC_ERRORS$n <- 0L
coloc_failed_report <- function() {
  if (COLOC_ERRORS$n > 0L)
    message("WARNING: ", COLOC_ERRORS$n, " coloc.susie comparison(s) errored in ",
            length(COLOC_ERRORS$units), " unit(s) that still kept ok or no_sig_coloc rows: ",
            paste(head(COLOC_ERRORS$units, 20), collapse = ", "))
}
coloc_acc_rows <- function(acc, no_cs_status, stub, unit = "") {
  n_err <- sum(acc$fails == "coloc_error")
  if (acc$ran_any && n_err) {
    COLOC_ERRORS$units <- c(COLOC_ERRORS$units, unit)
    COLOC_ERRORS$n <- COLOC_ERRORS$n + n_err
  }
  if (length(acc$rows)) return(rbindlist(acc$rows, fill = TRUE))
  st <- pair_status(acc$ran_any, acc$fails, no_cs_status)
  if (st == "no_sig_coloc" && !is.null(acc$best)) acc$best else stub(st)
}

# The no-CS status of an x_disease pair. The pair builder put the pair on the list
# for its disease clump leads in the exposure's cis window (o_leads, ";"-joined).
# If one of those leads has a credible set but no disease locus is near an exposure
# locus (loci_near), nothing was tested although the disease
# does have a CS there: "disease_cs_outside_window", not "no_disease_cs" (handoff
# r010, ARHGAP27 x PD). Neither is a tested status.
disease_no_cs_status <- function(o_leads, oloc) {
  leads <- strsplit(as.character(o_leads), ";", fixed = TRUE)[[1L]]
  if (!is.null(oloc) && any(leads %chin% oloc$clump_lead_snp)) "disease_cs_outside_window"
  else "no_disease_cs"
}

# ---- shared by the two x_disease pair builders --------------------------------

# Disease clump leads (CHR, BP, SNP, outcome_id, outcome_source) from named clump dirs,
# e.g. list("FinnGen meta-analysis" = <dir>, "Other" = <dir>). A header-only
# .clumped (no significant SNP) contributes nothing.
read_disease_leads <- function(clump_dirs) {
  one <- function(src_name, dir) {
    files <- list.files(dir, pattern = "\\.clumped$", full.names = TRUE)
    if (!length(files)) return(NULL)
    rbindlist(lapply(files, function(f) {
      if (file.size(f) == 0) return(NULL)
      # Stop on a read error: NULL would drop all of the outcome's pairs unseen.
      d <- tryCatch(fread(f, select = c("CHR", "BP", "SNP"), showProgress = FALSE),
                    error = function(e) stop("cannot read ", f, ": ", conditionMessage(e)))
      if (!nrow(d)) return(NULL)
      d[, `:=`(outcome_id = sub("\\.clumped$", "", basename(f)), outcome_source = src_name)]
    }), fill = TRUE)
  }
  d <- rbindlist(mapply(one, names(clump_dirs), clump_dirs, SIMPLIFY = FALSE), fill = TRUE)
  d <- d[!is.na(CHR) & !is.na(BP)]
  d[, CHR := as.integer(CHR)]
  d
}

# Write <info_dir>/coloc_pair_todo.csv: the pairs not yet in the ledger, plus those
# whose ledger rows are older than one of their fine-maps. `pair` carries the key
# columns and input_mtime (epoch seconds of the newer fine-map of the two sides).
# Un-imported chunks in tmp_dir count as done: a run that finished before its collect
# has not been lost. The ledger survives the collect's chunk cleanup, which is why
# done is read from it rather than from tmp_chunks/ alone.
write_coloc_todo <- function(pair, db_file, db_table, key_cols, tmp_dir, info_dir) {
  key_of <- function(d) do.call(paste, c(lapply(key_cols, function(k) as.character(d[[k]])), sep = KEY_SEP))
  done_keys <- db_done_keys(db_file, db_table, key_cols)
  msg("ledger done: ", length(done_keys))
  chunk_files <- list.files(tmp_dir, pattern = "^chunk_.*\\.tsv$", full.names = TRUE)
  if (length(chunk_files)) {
    done <- rbindlist(lapply(chunk_files, fread, select = c(key_cols, "computed_at")), fill = TRUE)
    # Same staleness rule as the ledger: a chunk row older than its fine-map is not done.
    mt <- pair[, c(key_cols, "input_mtime"), with = FALSE]
    for (k in key_cols) { done[, (k) := as.character(get(k))]; mt[, (k) := as.character(get(k))] }
    done <- merge(done, mt, by = key_cols, all.x = TRUE)
    done <- done[is.na(input_mtime) | (!is.na(computed_at) & computed_at >= input_mtime)]
    done_keys <- unique(c(done_keys, key_of(done)))
    msg("plus un-imported chunks: ", length(chunk_files), " file(s)")
  }
  stale_keys <- db_stale_keys(db_file, db_table, key_cols, pair[, c(key_cols, "input_mtime"), with = FALSE])
  msg("stale (a fine-map is newer than its ledger rows): ", length(stale_keys))
  todo <- pair[!key_of(pair) %in% setdiff(done_keys, stale_keys)]
  todo[, input_mtime := NULL]
  out_todo <- file.path(info_dir, "coloc_pair_todo.csv")
  fwrite(todo, out_todo)
  msg("todo: ", nrow(todo), "/", nrow(pair), " -> ", out_todo)
  invisible(todo)
}
