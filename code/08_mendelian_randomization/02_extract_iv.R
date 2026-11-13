#!/usr/bin/env Rscript

# Extract the MR instruments once per run, for every exposure in mr_pair_todo.csv.gz.
# IVs: clump lead SNPs (5e-8/1000kb/0.01r2) restricted to the cis window (follow
#      Eldjarn et al.); effect sizes pulled from clean protein sumstats.
# Each 03 chunk used to read every exposure's clean file in full (7.5M rows) to pull a
# few SNPs: ~5.8 h of an ~8 h chunk, repeated in all 80 chunks (~460 CPU-h, 2026-10-01).
# Output: info/mr_iv.tsv.gz, one row per (exposure, IV); an exposure without an IV
#         gets one row with SNP NA, so 03 can tell "no IV" from "not extracted".
#         Merged, not overwritten: this run's exposures replace their rows, all others
#         are kept. A todo-only run (2026-10-23, 18 multi-subunit exposures) used to
#         overwrite the table, losing the instruments behind every earlier MR result.
# Env:    MR_PAIR_FILE, MR_IV_FILE override the paths (tests); cores from
#         SLURM_CPUS_PER_TASK.

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

INFO_DIR  <- path.expand("./info")
source(path.expand("code/lib/common.R"))   # GWAS_BASE, CLUMP_SUBDIR, msg
source(path.expand("code/lib/mhc.R"))
MR_DIR    <- file.path(MR_ROOT, "protein2disease")
PAIR_FILE <- Sys.getenv("MR_PAIR_FILE", file.path(MR_DIR, "info", "mr_pair_todo.csv.gz"))
IV_FILE   <- Sys.getenv("MR_IV_FILE",   file.path(MR_DIR, "info", "mr_iv.tsv.gz"))
# When each exposure's rows were extracted (exposure_id, extracted_at in epoch seconds).
# The reuse rule needs a time per exposure: IV_FILE's own mtime moves on every run, so
# an exposure rebuilt before an unrelated run looked current afterwards.
TIME_FILE <- sub("\\.tsv\\.gz$", "_times.tsv.gz", IV_FILE)
CORES     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
t0        <- proc.time()[[3L]]

# Per-assay cis windows, one per subunit (gene +/-500kb); windows of a gene in the MHC
# are dropped (../lib/mhc.R), as in 02.
source(file.path(REPO_ROOT, "code/lib/cis_windows.R"))
cis_win <- read_cis_windows(REPO_ROOT)
cis_win <- cis_win[!mhc_cis_window(chr, cis_start, cis_end)]

cohort_of   <- function(src) if (src == "BioFINDER1") "bf1" else "bf2"
clump_path  <- function(src, id) file.path(GWAS_BASE, cohort_of(src),
  "clump/protein", CLUMP_SUBDIR, paste0(id, ".clumped"))
clean_path  <- function(src, id) file.path(GWAS_BASE, cohort_of(src),
  "clean/protein", paste0(id, ".tsv.gz"))

load_iv <- function(id, src) {
  # Clump-based cis-pQTL instruments (follow Eldjarn et al.): independent
  # (r2<0.01) genome-wide-sig clump leads inside the cis window. Cis-only avoids
  # trans pleiotropy; clump r2<0.01 guarantees IV independence (IVW assumption),
  # unlike fine-mapping PIP-max which can split one signal into LD-correlated CS.
  cf <- clump_path(src, id)
  # A missing .clumped FILE means this assay had no significant SNPs -- legitimate.
  # A missing clump DIRECTORY means CLUMP_SUBDIR is stale: every pair would return
  # NULL and the run would report success with zero results. Fail loudly instead.
  if (!dir.exists(dirname(cf))) stop("clump dir not found: ", dirname(cf), call. = FALSE)
  if (!file.exists(cf)) return(NULL)
  cl <- fread(cf)
  if (!nrow(cl)) return(NULL)
  pa <- sub("^bf[12]__", "", id)
  in_cis <- in_cis_window(rep(pa, nrow(cl)), cl$CHR, cl$BP, cis_win)
  cl <- cl[in_cis & !mhc_variant(cl$CHR, cl$BP)]
  if (!nrow(cl)) return(NULL)
  ss <- fread(clean_path(src, id),
              select = c("SNP", "EA", "NEA", "BETA", "SE", "P"))
  d <- ss[SNP %in% cl$SNP]
  if (!nrow(d)) return(NULL)
  d[, .(SNP = SNP, EA_e = EA, NEA_e = NEA,
        BETA_e = BETA, SE_e = SE, P_e = P)]
}

expo <- unique(fread(PAIR_FILE, select = c("exposure_id", "exposure_source")))
if (anyDuplicated(expo$exposure_id)) stop("an exposure_id has more than one exposure_source")
# Reuse an exposure's rows in IV_FILE when its clean and clump files and the cis windows
# are all older than its own extraction: a new outcome pairs with every exposure, and
# re-extracting all of them gave the same table (2026-10-30). Delete IV_FILE to force a
# full run (for a change in this script or ../lib/mhc.R, which no mtime shows).
times <- if (file.exists(TIME_FILE)) fread(TIME_FILE) else
  data.table(exposure_id = character(), extracted_at = numeric())
if (file.exists(IV_FILE)) {
  have <- data.table(exposure_id = unique(fread(IV_FILE, select = "exposure_id")$exposure_id))
  # Rows written before TIME_FILE existed date from IV_FILE's mtime at the latest.
  have[, t := times$extracted_at[match(exposure_id, times$exposure_id)]]
  have[is.na(t), t := as.numeric(file.mtime(IV_FILE))]
  mt <- function(f) if (file.exists(f)) as.numeric(file.mtime(f)) else -Inf
  win_t <- mt(file.path(REPO_ROOT, "info/bf_protein/ensg/assay_cis_windows.csv"))
  t_ex  <- have$t[match(expo$exposure_id, have$exposure_id)]
  src_t <- mapply(function(id, src) max(mt(clean_path(src, id)), mt(clump_path(src, id)), win_t),
                  expo$exposure_id, expo$exposure_source)
  reuse <- !is.na(t_ex) & src_t <= t_ex
  msg(sum(reuse), " of ", nrow(expo), " exposures reused from ", IV_FILE)
  expo <- expo[!reuse]
}
# Nothing to extract: touch IV_FILE, since 03 requires it to be newer than the todo list.
if (!nrow(expo)) {
  if (!Sys.setFileTime(IV_FILE, Sys.time())) stop("could not touch ", IV_FILE)
  msg("no exposure to extract; ", IV_FILE, " content unchanged, mtime updated")
  quit(save = "no")
}
msg(nrow(expo), " exposures, ", CORES, " core(s)")

run_start <- as.numeric(Sys.time())   # a source changed during the run counts as newer
res <- mclapply(seq_len(nrow(expo)), function(i) {
  id <- expo$exposure_id[i]
  iv <- load_iv(id, expo$exposure_source[i])
  if (is.null(iv)) iv <- data.table(SNP = NA_character_)
  cbind(exposure_id = id, iv)
}, mc.cores = CORES)
# A forked child that errored returns a try-error; it must not become "no IV".
bad <- which(!vapply(res, is.data.table, logical(1)))
if (length(bad)) stop(length(bad), " exposure(s) failed, e.g. ", expo$exposure_id[bad[1]], ": ",
                      as.character(res[[bad[1]]]))
iv <- rbindlist(res, fill = TRUE)
if (!setequal(unique(iv$exposure_id), expo$exposure_id)) stop("IV table does not cover every exposure")

if (file.exists(IV_FILE)) {
  old <- fread(IV_FILE, colClasses = c(SNP = "character", EA_e = "character", NEA_e = "character"),
              na.strings = c("", "NA"))   # fwrite writes the no-IV row's NA SNP as ""
  # A run whose exposures all lack an instrument has only exposure_id and SNP.
  miss <- setdiff(names(old), names(iv))
  if (length(miss) && all(is.na(iv$SNP))) iv[, (miss) := lapply(miss, function(cc) old[[cc]][NA_integer_])]
  if (!setequal(names(old), names(iv))) stop("columns of ", IV_FILE, " differ from this run's")
  old <- old[!exposure_id %in% expo$exposure_id]
  msg("kept ", uniqueN(old$exposure_id), " exposures from ", IV_FILE, ", replaced ", nrow(expo))
  iv <- rbind(old, iv, use.names = TRUE)
}

tmp <- paste0(IV_FILE, ".tmp.gz")
fwrite(iv, tmp, sep = "\t", compress = "gzip")
if (!file.rename(tmp, IV_FILE)) stop("could not rename ", tmp, " -> ", IV_FILE)
# Times after the table, so a crash in between re-extracts rather than trusts.
times <- rbind(times[!exposure_id %in% expo$exposure_id],
               data.table(exposure_id = expo$exposure_id, extracted_at = run_start))
times <- times[exposure_id %in% iv$exposure_id]
tmp <- paste0(TIME_FILE, ".tmp.gz")
fwrite(times, tmp, sep = "\t", compress = "gzip")
if (!file.rename(tmp, TIME_FILE)) stop("could not rename ", tmp, " -> ", TIME_FILE)
msg("wrote ", iv[!is.na(SNP), .N], " IVs for ", iv[!is.na(SNP), uniqueN(exposure_id)], " of ",
    uniqueN(iv$exposure_id), " exposures -> ", IV_FILE)
