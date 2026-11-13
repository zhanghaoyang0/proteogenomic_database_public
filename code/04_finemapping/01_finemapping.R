#!/usr/bin/env Rscript
# Fine mapping with SuSiE (genome-wide: cis + trans).
# Output includes is_cis column; filter downstream by is_cis as needed.
# Usage: Rscript 01_finemapping.R <clump_file> <clean_dir> <ref_bfile> <out_dir> <plink_bin>

library(susieR)
library(data.table)

args       <- commandArgs(trailingOnly = TRUE)
clump_file <- args[1]
clean_dir  <- args[2]
ref_bfile  <- args[3]
out_dir    <- args[4]
plink_bin  <- args[5]

source(path.expand("code/lib/common.R"))   # REPO_ROOT, FM_WINDOW_BP, P_GWS
WINDOW   <- FM_WINDOW_BP
P_THRESH <- P_GWS
# MHC loci are not fine-mapped, for every resource (rule: ../lib/mhc.R).
source(path.expand("code/lib/mhc.R"))

trait <- sub("\\.clumped$", "", basename(clump_file))
cat("Processing:", trait, "\n")

# Read clump file (all loci, cis + trans)
clump <- fread(clump_file)
# cis windows of this assay, one per subunit (none for eQTL/disease traits: is_cis FALSE)
source(file.path(REPO_ROOT, "code/lib/cis_windows.R"))
pa_str     <- sub("^bf[12]__", "", trait)
cis_win    <- read_cis_windows(REPO_ROOT)[protein_assay == pa_str]
if (nrow(clump) == 0) {
  cat("No lead SNPs, skipping.\n")
  quit(status = 0)
}

write_status <- function(x) {
  status_file <- file.path(out_dir, paste0(trait, ".finemap_status.tsv.gz"))
  tmp_status <- paste0(status_file, ".tmp")
  fwrite(x, tmp_status, sep = "\t", compress = "gzip")
  if (!file.rename(tmp_status, status_file)) stop("could not rename status file: ", tmp_status)
  cat("Locus status:", paste(x[, .N, by = finemap_status][,
                                paste0(finemap_status, "=", N)], collapse = ", "), "\n")
}

# Read clean GWAS
# Protein cohorts are flat: <clean_dir>/<trait>.tsv.gz. The eQTL resources moved to a
# per-tissue subdir layout (<clean_dir>/<tissue>/<ENSG>.tsv.gz) and this resolver was
# never updated. GTEx and MetaBrain had already been fine-mapped under the flat layout,
# so nothing failed until SingleBrain -- the first resource fine-mapped after the
# migration -- where all 6,007 traits errored one by one while the array still reported
# COMPLETED. Try flat first, then the subdir form for "<tissue>_<ENSG...>" traits.
clean_file <- file.path(clean_dir, paste0(trait, ".tsv.gz"))
if (!file.exists(clean_file) && grepl("_ENSG", trait, fixed = TRUE)) {
  tissue <- sub("_ENSG.*$", "", trait)
  gene   <- substring(trait, nchar(tissue) + 2L)
  alt    <- file.path(clean_dir, tissue, paste0(gene, ".tsv.gz"))
  if (file.exists(alt)) clean_file <- alt
}
if (!file.exists(clean_file)) stop("Clean file not found: ", clean_file)
ss <- fread(clean_file)
ss <- ss[!is.na(BETA) & !is.na(SE) & SE > 0]

# Temp dir for the plink files: node-local disk inside a LUNARC job ($SNIC_TMP, e.g.
# /local/slurmtmp.<jobid>; SLURM_TMPDIR is not set there), else OUT_DIR/tmp on the
# shared FS. On the login node SNIC_TMP is a relative placeholder, hence the check.
# Cleared first, so files a killed run left behind are never read.
tmp_base <- Sys.getenv("SNIC_TMP")
if (!startsWith(tmp_base, "/") || !dir.exists(tmp_base)) tmp_base <- file.path(out_dir, "tmp")
tmp_dir <- file.path(tmp_base, trait)
unlink(tmp_dir, recursive = TRUE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
cat("plink temp dir:", tmp_dir, "\n")

# Load reference BIM once per trait. Each trait is its own Rscript (xargs in the sbatch),
# so "once per job" would need one long-lived R process for many traits; a few seconds
# per trait were not worth that (code review 2026-11-01, C47).
bim_all <- fread(paste0(ref_bfile, ".bim"),
                 col.names = c("CHR", "SNP_ref", "CM", "BP", "A1", "A2"))

# Fine-map one locus (independent of all others) -> list of per-CS data.tables.
# Pulled out of the loop so loci can be run in parallel via mclapply (big traits
# like T2D have 1000+ loci; within-trait parallelism is the main speedup).
# Each locus uses its own tmp files (indexed by i) and its own RDS, so forks
# never collide. Returns both credible-set rows and one explicit locus-status row.
fmap_locus <- function(i) {
  lead_snp <- clump$SNP[i]
  result <- function(status, detail = "", cs = list(), n_snps = NA_integer_, n_cs = length(cs)) {
    detail <- gsub("[\r\n\t]+", " ", detail)
    list(cs = cs,
         status = data.table(trait = trait, clump_lead_snp = lead_snp,
                             finemap_status = status, detail = detail,
                             n_snps = as.integer(n_snps), n_cs = as.integer(n_cs)))
  }
  # Every bail-out is "this locus is unusable, here is why": one place to say it,
  # so the reasons stay uniform and greppable in the array logs.
  skip <- function(why, status = "insufficient_data") {
    cat("  ", lead_snp, ": ", why, ", skipping\n", sep = "")
    result(status, why)
  }
  lead_chr <- clump$CHR[i]
  lead_bp  <- clump$BP[i]
  if (mhc_variant(lead_chr, lead_bp))
    return(skip("lead is in the MHC, not fine-mapped", "mhc"))

  # cis is decided per credible set, by its top-PIP variant (below), not by the clump
  # lead: one +-500 kb window can hold several credible sets at different positions,
  # and the top-PIP variant is the best estimate of where each signal sits. (Until
  # 2026-10-06 every set took the clump lead's status; 348 of 12,887 BF2 sets differed.)
  in_cis <- function(chr, pos) in_cis_window(pa_str, chr, pos, cis_win)

  reg_start <- max(1, lead_bp - WINDOW)
  reg_end   <- lead_bp + WINDOW

  # Summary stats in window, remove duplicates by position
  ss_win <- ss[CHR == lead_chr & POS >= reg_start & POS <= reg_end]
  ss_win[, pos_key := paste(CHR, POS, sep = ":")]
  ss_win <- ss_win[!duplicated(pos_key)]

  # SuSiE-RSS takes one n for the whole region. In a meta-analysis (MetaBrain) N varies
  # per SNP -- in cortex windows the smallest is ~4% of the largest -- and z-scores from
  # SNPs seen in only a few cohorts disagree with the LD the reference implies, which
  # left 43% of MetaBrain cognate fits unconverged (2026-10-02). Keep SNPs with
  # N >= 80% of the window maximum, as FinnGen and PolyFun do. BF2, GTEx and
  # SingleBrain have one N per trait, so this drops nothing there.
  if ("N" %in% names(ss_win) && any(is.finite(ss_win$N))) {
    n_max  <- max(ss_win$N, na.rm = TRUE)
    keep_n <- (ss_win$N >= 0.8 * n_max) %in% TRUE   # a missing N is dropped too
    if (any(!keep_n)) cat(sprintf("  dropped %d of %d SNPs with N < 80%% of max N (%d) for %s\n",
                                  sum(!keep_n), nrow(ss_win), as.integer(n_max), lead_snp))
    ss_win <- ss_win[keep_n]
  }

  if (nrow(ss_win) < 50)
    return(skip(sprintf("too few summary-stat SNPs in the window (n=%d)", nrow(ss_win))))

  # 1000G BIM in window, remove duplicates by position
  bim_win <- bim_all[CHR == lead_chr & BP >= reg_start & BP <= reg_end]
  bim_win[, pos_key := paste(CHR, BP, sep = ":")]
  bim_win <- bim_win[!duplicated(pos_key)]

  # Intersect by CHR:POS
  common  <- intersect(ss_win$pos_key, bim_win$pos_key)
  ss_win  <- ss_win[pos_key %in% common]
  bim_win <- bim_win[pos_key %in% common]

  if (nrow(bim_win) < 50)
    return(skip(sprintf("too few SNPs shared with the LD reference (n=%d)", nrow(bim_win))))

  # Step 1: plink --make-bed to get exact genotype subset with true SNP order
  snp_file       <- file.path(tmp_dir, paste0("snps_", i, ".txt"))
  tmp_ref_prefix <- file.path(tmp_dir, paste0("ref_", i))
  tmp_ld_prefix  <- file.path(tmp_dir, paste0("ld_", i))

  fwrite(data.table(SNP = unique(bim_win$SNP_ref)), snp_file, col.names = FALSE)

  cmd_makebed <- sprintf(
    "%s --bfile %s --chr %d --from-bp %d --to-bp %d --extract %s --make-bed --out %s --silent",
    plink_bin, ref_bfile, lead_chr, reg_start, reg_end, snp_file, tmp_ref_prefix
  )
  # A plink failure is an I/O or environment fault, not "no signal": stop, so the
  # worker error fails the trait (see the try-error check after mclapply) instead of
  # recording the locus as having no credible set. Never seen in 4,592 logs to date.
  if (system(cmd_makebed) != 0) stop("plink --make-bed failed for ", lead_snp)

  tmp_bim_file <- paste0(tmp_ref_prefix, ".bim")
  if (!file.exists(tmp_bim_file)) stop("plink --make-bed produced no BIM for ", lead_snp)

  # This BIM is the true SNP order used by plink
  plink_bim <- fread(tmp_bim_file,
                     col.names = c("CHR", "SNP_ref", "CM", "BP", "A1", "A2"))
  plink_bim[, pos_key := paste(CHR, BP, sep = ":")]

  if (nrow(plink_bim) < 50)
    return(skip(sprintf("too few SNPs in the plink subset (n=%d)", nrow(plink_bim))))

  # Step 2: compute signed LD matrix
  cmd_ld <- sprintf(
    "%s --bfile %s --r square spaces --out %s --silent",
    plink_bin, tmp_ref_prefix, tmp_ld_prefix
  )
  if (system(cmd_ld) != 0) stop("plink --r failed for ", lead_snp)

  ld_file <- paste0(tmp_ld_prefix, ".ld")
  if (!file.exists(ld_file)) stop("plink --r produced no LD matrix for ", lead_snp)

  R_mat <- as.matrix(fread(ld_file, header = FALSE))

  if (nrow(R_mat) != nrow(plink_bim) || ncol(R_mat) != nrow(plink_bim))
    stop(sprintf("%s: LD matrix is %dx%d but the BIM has %d SNPs",
                 lead_snp, nrow(R_mat), ncol(R_mat), nrow(plink_bim)))

  # Match ss to plink BIM order
  idx   <- match(plink_bim$pos_key, ss_win$pos_key)
  valid <- !is.na(idx)

  plink_bim  <- plink_bim[valid]
  R_mat      <- R_mat[valid, valid, drop = FALSE]
  ss_ordered <- ss_win[idx[valid]]

  if (nrow(ss_ordered) < 50)
    return(skip(sprintf("too few SNPs left after matching to the plink order (n=%d)",
                        nrow(ss_ordered))))

  # Allele alignment: keep only SNPs with matching allele pairs
  same <- ss_ordered$EA == plink_bim$A1 & ss_ordered$NEA == plink_bim$A2
  flip <- ss_ordered$EA == plink_bim$A2 & ss_ordered$NEA == plink_bim$A1
  keep <- same | flip

  if (sum(keep) < 50)
    return(skip(sprintf("too few allele-matched SNPs (n=%d)", sum(keep))))

  ss_ordered <- ss_ordered[keep]
  plink_bim  <- plink_bim[keep]
  R_mat      <- R_mat[keep, keep, drop = FALSE]

  z <- ss_ordered$BETA / ss_ordered$SE
  z[flip[keep]] <- -z[flip[keep]]

  stopifnot(length(z) == nrow(R_mat))

  n_eff <- as.integer(median(ss_ordered$N, na.rm = TRUE))
  if (is.na(n_eff) || n_eff <= 0) return(skip("no usable sample size (N)"))

  # Force symmetry
  R_mat <- (R_mat + t(R_mat)) / 2
  diag(R_mat) <- 1

  # Run SuSiE with L=10, filter by credible set purity. max_iter: the default 100 left
  # 313 of 10,602 BF2 fits unconverged (2026-10-02 rerun), mostly strong loci with
  # many independent signals; 1000 lets them finish. Non-convergence is still logged.
  fit_error <- NULL
  fit <- tryCatch(
    susie_rss(z = z, R = R_mat, n = n_eff, L = 10, max_iter = 1000, verbose = FALSE),
    error = function(e) {
      fit_error <<- conditionMessage(e)
      cat("  susie_rss error for", lead_snp, ":", fit_error, "\n")
      NULL
    }
  )

  if (!is.null(fit) && !isTRUE(fit$converged))
    cat("  WARNING: SuSiE did not converge in", fit$niter, "iterations for", lead_snp, "\n")

  if (!is.null(fit_error)) return(result("susie_error", fit_error, n_snps = length(z)))

  if (is.null(fit) || length(fit$sets$cs) == 0) {
    cat("  SuSiE found no credible sets for", lead_snp,
        "(n_snps=", length(z), ")\n")
    return(result("no_credible_set", "SuSiE returned no credible set", n_snps = length(z)))
  }

  # Save SuSiE fit object for coloc.susie() downstream
  rds_dir <- file.path(out_dir, "susie_rds")
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  rds_file <- file.path(rds_dir, paste0(trait, "__", lead_snp, ".susie.rds"))
  fit$snp <- ss_ordered$SNP  # attach SNP names so coloc can match
  tmp_rds <- paste0(rds_file, ".tmp.", Sys.getpid())
  saveRDS(fit, tmp_rds)
  if (!file.rename(tmp_rds, rds_file)) stop("could not rename ", tmp_rds, " -> ", rds_file)

  pip_vec  <- fit$pip
  cs_list  <- fit$sets$cs
  purity   <- fit$sets$purity  # data.frame: min.abs.corr, mean.abs.corr, median.abs.corr

  locus_cs <- list()
  for (cs_id in seq_along(cs_list)) {
    cs_idx <- cs_list[[cs_id]]

    # Keep only GWS credible sets
    if (min(ss_ordered$P[cs_idx], na.rm = TRUE) > P_THRESH) next

    cs_purity <- if (!is.null(purity) && cs_id <= nrow(purity))
      purity[cs_id, "min.abs.corr"] else NA_real_
    cs_top_idx <- cs_idx[which.max(pip_vec[cs_idx])]

    locus_cs[[length(locus_cs) + 1]] <- data.table(
      trait       = trait,
      clump_lead_snp = lead_snp,
      lead_snp    = ss_ordered$SNP[cs_top_idx],
      cs_id       = cs_id,
      cs_purity   = cs_purity,
      is_cis      = in_cis(ss_ordered$CHR[cs_top_idx], ss_ordered$POS[cs_top_idx]),
      chr         = ss_ordered$CHR[cs_idx],
      pos         = ss_ordered$POS[cs_idx],
      snp         = ss_ordered$SNP[cs_idx],
      ea          = ss_ordered$EA[cs_idx],
      nea         = ss_ordered$NEA[cs_idx],
      beta        = ss_ordered$BETA[cs_idx],
      se          = ss_ordered$SE[cs_idx],
      p           = ss_ordered$P[cs_idx],
      pip         = pip_vec[cs_idx],
      # A non-converged fit still yields credible sets; coloc's primary result drops
      # them (coloc/functions_coloc.R add_component_flags), the table only flags them.
      converged   = isTRUE(fit$converged)
    )
  }
  if (length(locus_cs)) result("credible_set", cs = locus_cs, n_snps = length(z))
  else result("no_credible_set", "no credible set passed the genome-wide-significance filter",
              n_snps = length(z))
}

# Run loci: parallel within-trait when FINEMAP_LOCI_CORES > 1 (default 1 = old
# sequential behavior). mclapply forks share bim_all/ss copy-on-write; per-locus
# R_mat peaks ~70MB so even 16 cores stay well within job memory.
loci_cores <- suppressWarnings(as.integer(Sys.getenv("FINEMAP_LOCI_CORES", "1")))
if (is.na(loci_cores) || loci_cores < 1) loci_cores <- 1
cat("Fine-mapping", nrow(clump), "loci with", loci_cores, "core(s)\n")
# Optional split of one trait's loci over several jobs (T2D: 1,638 loci took ~13 h in one
# 16-core task). FINEMAP_CHUNK=k/n fits loci k, k+n, k+2n, ... and saves their results to
# <out_dir>/chunks/<trait>.chunk_<k>_of_<n>.rds, then exits. FINEMAP_CHUNKS_FROM=n fits
# nothing: it puts the n chunks back in locus order, so everything below (dedup across
# windows, cs_id numbering, output) sees the same `res` as an unsplit run.
# submit_01_finemapping_split.sh drives both.
chunk       <- Sys.getenv("FINEMAP_CHUNK", "")
chunks_from <- suppressWarnings(as.integer(Sys.getenv("FINEMAP_CHUNKS_FROM", "")))
chunk_file  <- function(k, n) file.path(out_dir, "chunks", sprintf("%s.chunk_%d_of_%d.rds", trait, k, n))
idx <- seq_len(nrow(clump))
if (nzchar(chunk)) {
  kn <- suppressWarnings(as.integer(strsplit(chunk, "/", fixed = TRUE)[[1]]))
  if (length(kn) != 2L || anyNA(kn) || kn[1] < 1L || kn[1] > kn[2]) stop("FINEMAP_CHUNK must be k/n, got: ", chunk)
  idx <- idx[(idx - 1L) %% kn[2] == kn[1] - 1L]
  cat("Chunk", chunk, ":", length(idx), "of", nrow(clump), "loci\n")
}
if (!is.na(chunks_from)) {
  res <- vector("list", nrow(clump))
  for (k in seq_len(chunks_from)) {
    f <- chunk_file(k, chunks_from)
    if (!file.exists(f)) stop("missing chunk ", k, "/", chunks_from, ": ", f)
    part <- readRDS(f)
    # the staging must not have changed between the chunk runs and this one
    if (!identical(part$snps, clump$SNP[part$idx])) stop("chunk ", f, " was fitted on other loci than ", clump_file)
    res[part$idx] <- part$res
  }
  cat("Assembled", chunks_from, "chunks\n")
} else if (loci_cores > 1) {
  suppressPackageStartupMessages(library(parallel))
  res <- mclapply(idx, fmap_locus, mc.cores = loci_cores)
} else {
  res <- lapply(idx, fmap_locus)
}
# A child that errored (mclapply returns a try-error, or NULL if it was killed)
# would silently drop its loci from an otherwise "complete" output: stop instead.
bad <- which(!vapply(res, function(x) is.list(x) && !inherits(x, "try-error") &&
                               all(c("cs", "status") %in% names(x)), logical(1)))
if (length(bad)) stop(trait, ": ", length(bad), " locus job(s) failed, e.g. locus ", bad[1], ": ",
                      if (inherits(res[[bad[1]]], "try-error")) as.character(res[[bad[1]]]) else "NULL")
if (nzchar(chunk)) {
  f <- chunk_file(kn[1], kn[2]); dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp_chunk <- paste0(f, ".tmp.", Sys.getpid())
  saveRDS(list(idx = idx, snps = clump$SNP[idx], res = res), tmp_chunk)
  if (!file.rename(tmp_chunk, f)) stop("could not rename ", tmp_chunk, " -> ", f)
  unlink(tmp_dir, recursive = TRUE)
  cat("Chunk saved:", f, "\n")
  quit(status = 0)
}
all_cs <- unlist(lapply(res, `[[`, "cs"), recursive = FALSE)
locus_status <- rbindlist(lapply(res, `[[`, "status"))
if (nrow(locus_status) != nrow(clump) || anyDuplicated(locus_status[, .(trait, clump_lead_snp)]))
  stop("locus status is not exactly one row per clump lead")

# Cleanup plink temp files
unlink(tmp_dir, recursive = TRUE)

# Save output
final_cs_counts <- data.table(clump_lead_snp = character(), n_cs = integer())
# The finemapped table's columns, for both the table and its header-only form.
OUT_COLS <- c("trait", "clump_lead_snp", "lead_snp", "cs_id", "cs_purity",
              "is_cis", "chr", "pos", "snp", "ea", "nea", "beta", "se", "p", "pip", "converged")
if (length(all_cs) > 0) {
  out <- rbindlist(all_cs)
  out[, run_id := .GRP, by = .(trait, clump_lead_snp, cs_id)]
  cs_key <- out[, .(
    cs_key = paste(sort(unique(snp)), collapse = ";"),
    min_p = min(p, na.rm = TRUE),
    max_pip = max(pip, na.rm = TRUE)
  ), by = run_id]
  setorder(cs_key, cs_key, min_p, -max_pip)
  keep_run <- cs_key[, .SD[1], by = cs_key]$run_id
  out <- out[run_id %in% keep_run]
  out[, old_run_id := run_id]
  setorder(out, chr, pos, old_run_id, -pip)
  out[, cs_id := .GRP, by = old_run_id]
  out[, c("run_id", "old_run_id") := NULL]
  setcolorder(out, OUT_COLS)
  setorder(out, cs_id, -pip)
  final_cs_counts <- out[, .(n_cs = uniqueN(cs_id)), by = clump_lead_snp]
  n_cs <- uniqueN(out$cs_id)
  out_file <- file.path(out_dir, paste0(trait, ".finemapped.tsv.gz"))
  # tmp + rename: the skip test is this file, so a killed task must not leave it half-written.
  tmp_out <- paste0(out_file, ".tmp")
  fwrite(out, tmp_out, sep = "\t", compress = "gzip")
  cat("Done:", nrow(out), "variants across", n_cs, "credible sets saved.\n")
} else {
  cat("No credible sets passed filters for", trait, "\n")
  empty_cols <- OUT_COLS
  out_file <- file.path(out_dir, paste0(trait, ".finemapped.tsv.gz"))
  # fwrite() on a 0-row data.table writes header bytes but leaves gzip footer
  # unflushed -> `gzip -t` reports unexpected EOF. Write via gzfile() + close().
  tmp_out <- paste0(out_file, ".tmp")
  con <- gzfile(tmp_out, "w")
  writeLines(paste(empty_cols, collapse = "\t"), con)
  close(con)
}

# A separate sidecar keeps failures out of the credible-set variant table while
# allowing collectors to distinguish a model error from a genuine negative locus.
# The final table deduplicates identical CS across overlapping clump windows, so
# reconcile status after that step: a window whose only CS was kept under a
# neighbour must not claim to carry a CS absent from its own output rows.
locus_status[finemap_status == "credible_set", n_cs := 0L]
if (nrow(final_cs_counts)) {
  locus_status[final_cs_counts, on = "clump_lead_snp", n_cs := i.n_cs]
}
locus_status[finemap_status == "credible_set" & n_cs == 0L,
             `:=`(finemap_status = "no_credible_set",
                  detail = "credible set duplicated by another clump window")]
# Only BF proteins have a reader (04_build_signal_tables.R). An eQTL or FinnGen trait is
# one file per (tissue, gene) or trait already; a sidecar each would add ~25k files.
if (grepl("^bf[12]__", trait)) write_status(locus_status)
# The table last: it is what the freshness check reads, so it must not be current
# while the sidecar is stale (a task killed between the two was skipped from then on).
stopifnot(file.rename(tmp_out, out_file))
# The chunks are cache for this assembly only; the output above replaces them.
if (!is.na(chunks_from)) unlink(vapply(seq_len(chunks_from), function(k) chunk_file(k, chunks_from), ""))
