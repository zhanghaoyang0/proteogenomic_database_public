#!/usr/bin/env Rscript
# Credible-set deduplication: which credible sets of a trait are the same signal.
#
# Overlapping +/-500 kb windows fine-map the same signal more than once. CFHR4 in
# CSF is the clearest case: 10 windows produced 28 credible sets, five of which
# were led by rs375046 -- one per window.
#
# Rule (Nature Communications operational criterion): two credible sets of one
# trait coming from DIFFERENT windows are non-independent if ANY variant pair
# across them is in LD (r2 > 0.1); retain the set whose lead variant is more
# significant. Applied greedily in ascending lead-P order, not as transitive
# connected components -- with A-B and B-C in LD but A-C not, components would
# collapse all three, while greedy keeps A, drops B, then judges C against A.
#
# Sets from the SAME window are never compared: SuSiE already declared them
# distinct effects of one model.
#
# No distance limit, on purpose (code review 2026-11-01, C33). Of the Olink sets merged
# into another, 64 have leads more than 1 Mb apart; most sit in known long-range LD
# (chr3 47-52 Mb at 3p21.31, chr20 across the centromere; in SomaScan also the 8p23
# inversion and chr11 across the centromere), some with r2 near 0.8 (ENTPD6, 4.9 Mb).
# A cap would count those as separate signals. Only 27 far merges have r2 < 0.2, 0.65%
# of the 4,163 Olink signals, which bounds what chance LD could have removed.
#
# Nothing is merged and nothing is deleted. Every credible set keeps its own
# variants, PIPs and locus plot; this only labels which one represents a signal.
# Output is a side table, so *.finemapped.tsv.gz is untouched.
#
# Output: info/bf_protein/cs_dedup_map.tsv.gz -- same place and kind as
# finemap_lead_snps.csv, small enough to version, and read by the Shiny app.
#
# Incremental: a trait is recomputed only when its .finemapped.tsv.gz changed since the
# last run. The record is tmp/finemap/cs_dedup_state.tsv (repo root) (trait, file mtime, r2, output), not the mtime
# of the output: that file is in git, so a checkout resets it. No state, a different
# r2 or output file, or FORCE=1 recomputes every trait (~1.5 h for ~5,400). A different ref_bfile
# needs FORCE=1.
#
# Usage: Rscript 03_dedup_credible_sets.R [out.tsv.gz] [ref_bfile] [plink_bin] [r2=0.1] [dir ...]
suppressMessages(library(data.table))
a <- commandArgs(TRUE)
here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
# "cs_dedup_map", not "cs_signal_map": its signal_group_id is a DEDUPLICATED
# GROUP of credible sets, a different partition from the discovery signals in
# info/bf_protein/signals.tsv.gz. Two files both offering a "signal_id" invited a
# merge() that would join silently wrong, so the file and the column are renamed
# together -- a stale reader now fails with "file not found" instead.
DEF_OUT <- normalizePath(file.path(here, "..", "..", "info", "bf_protein",
                                   "cs_dedup_map.tsv.gz"), mustWork = FALSE)
OUT   <- if (length(a) >= 1 && nzchar(a[1])) a[1] else DEF_OUT
REF   <- if (length(a) >= 2 && nzchar(a[2])) a[2] else
           LD_REF
PLINK <- if (length(a) >= 3 && nzchar(a[3])) a[3] else Sys.which("plink")
if (!nzchar(PLINK)) stop("plink not found: module load PLINK/1.9b_6.21-x86_64 first")
R2  <- if (length(a) >= 4) as.numeric(a[4]) else 0.1
DIRS <- if (length(a) >= 5) a[-(1:4)] else
  file.path(GWAS_BASE, "bf2", "finemapping/protein")

files <- unlist(lapply(DIRS, list.files, pattern = "[.]finemapped[.]tsv[.]gz$", full.names = TRUE))
cat(sprintf("traits: %d   r2 threshold: %.2f\n", length(files), R2))
trait_of <- function(f) sub("[.]finemapped[.]tsv[.]gz$", "", basename(f))
if (anyDuplicated(trait_of(files))) stop("a trait is fine-mapped in more than one of: ", paste(DIRS, collapse = " "))
# mtime as a string: a numeric mtime loses its fractional seconds through fwrite/fread.
cur <- data.table(trait = trait_of(files), mtime = format(file.mtime(files), "%Y-%m-%d %H:%M:%OS6"), r2 = R2,
                  out = normalizePath(OUT, mustWork = FALSE))

STATE <- Sys.getenv("CS_DEDUP_STATE", file.path(here, "..", "..", "tmp", "finemap", "cs_dedup_state.tsv"))
reuse <- character(0)
if (Sys.getenv("FORCE") != "1" && file.exists(STATE) && file.exists(OUT)) {
  st <- fread(STATE, colClasses = c(trait = "character", mtime = "character"))
  hit <- merge(cur, st, by = "trait", suffixes = c("", ".st"))
  reuse <- hit[mtime == mtime.st & r2 == r2.st & out == out.st, trait]
}
prev <- if (length(reuse)) fread(OUT, colClasses = c(trait = "character"))[trait %in% reuse] else NULL
files <- files[!trait_of(files) %in% reuse]
cat(sprintf("reused from the last run: %d   to compute: %d\n", length(reuse), length(files)))
tmp <- tempfile("csdedup"); dir.create(tmp)

one <- function(f) {
  fm <- fread(cmd = paste("zcat", shQuote(f)))
  if (!nrow(fm)) return(NULL)   # fine-mapped, no credible set
  trait <- sub("[.]finemapped[.]tsv[.]gz$", "", basename(f))
  cs_snps <- split(as.character(fm$snp), fm$cs_id)
  meta <- unique(fm[, .(cs_id, window = clump_lead_snp, lead_snp)])
  meta <- merge(meta, fm[snp == lead_snp, .(lead_p = min(p)), by = cs_id], by = "cs_id")
  meta[, size := lengths(cs_snps)[as.character(cs_id)]]
  setorder(meta, lead_p)

  # A trait whose sets all come from one window has nothing to compare, and one
  # set has nothing to compare against: skip plink entirely for those.
  if (nrow(meta) == 1L || uniqueN(meta$window) == 1L) {
    meta[, `:=`(signal_id = seq_len(.N), is_representative = TRUE,
                duplicate_of = NA_integer_, max_r2 = NA_real_)]
  } else {
    # tmp is shared across traits: clear the last trait's plink files, or a
    # failed run here would silently reuse its LD.
    unlink(list.files(tmp, pattern = "^(sub|ld)[.]", full.names = TRUE))
    uni <- unique(as.character(fm$snp))
    fwrite(data.table(uni), file.path(tmp,"s.txt"), col.names = FALSE)
    system(sprintf("%s --bfile %s --extract %s --make-bed --out %s --silent",
                   PLINK, REF, file.path(tmp,"s.txt"), file.path(tmp,"sub")), ignore.stderr = TRUE)
    bimf <- file.path(tmp,"sub.bim")
    # Without LD every set would come out representative, so a plink failure must stop.
    if (!file.exists(bimf)) stop(trait, ": plink --make-bed wrote no ", bimf)
    system(sprintf("%s --bfile %s --r2 square spaces --out %s --silent",
                   PLINK, file.path(tmp,"sub"), file.path(tmp,"ld")), ignore.stderr = TRUE)
    if (!file.exists(file.path(tmp,"ld.ld"))) stop(trait, ": plink --r2 wrote no ld.ld")
    bim <- fread(bimf, select = 2, col.names = "SNP")
    LD <- as.matrix(fread(file.path(tmp,"ld.ld"), header = FALSE))
    dimnames(LD) <- list(bim$SNP, bim$SNP)
    vin <- function(c) intersect(cs_snps[[as.character(c)]], rownames(LD))
    ids <- meta$cs_id
    st <- setNames(rep(NA, length(ids)), ids); du <- setNames(rep(NA_integer_, length(ids)), ids)
    dr <- setNames(rep(NA_real_, length(ids)), ids); sg <- setNames(rep(NA_integer_, length(ids)), ids)
    win <- setNames(meta$window, as.character(meta$cs_id)); n <- 0L
    for (cid in ids) {
      k <- as.character(cid); if (!is.na(st[k])) next
      n <- n + 1L; st[k] <- TRUE; sg[k] <- n
      ix <- vin(cid); if (!length(ix)) next
      for (oth in ids) {
        o <- as.character(oth)
        if (!is.na(st[o]) || win[o] == win[k]) next
        iy <- vin(oth); if (!length(iy)) next
        r <- suppressWarnings(max(LD[ix, iy, drop = FALSE], na.rm = TRUE))
        if (is.finite(r) && r > R2) { st[o] <- FALSE; du[o] <- cid; dr[o] <- r; sg[o] <- n }
      }
    }
    meta[, `:=`(signal_id = sg[as.character(cs_id)], is_representative = st[as.character(cs_id)],
                duplicate_of = du[as.character(cs_id)], max_r2 = round(dr[as.character(cs_id)], 3))]
    unlink(list.files(tmp, full.names = TRUE))
  }
  cbind(trait = trait, meta)
}

res <- vector("list", length(files)); t0 <- Sys.time()
for (i in seq_along(files)) {
  res[[i]] <- one(files[i])
  if (i %% 500 == 0) cat(sprintf("  %d/%d  (%.1f min)\n", i, length(files),
                                 as.numeric(difftime(Sys.time(), t0, units="mins"))))
}
out <- rbindlist(res, fill = TRUE)
if (nrow(out)) {
  setcolorder(out, c("trait","signal_id","cs_id","window","lead_snp","lead_p","size",
                     "is_representative","duplicate_of","max_r2"))
  # Written-out names differ from the working ones on purpose: `signal_id` here is
  # a dedup GROUP, and `window` is the clump lead the window was centred on. Both
  # names clashed with the discovery-signal meaning used by signals.tsv.gz.
  setnames(out, c("signal_id", "window", "duplicate_of"),
                c("signal_group_id", "clump_lead_snp", "duplicate_of_cs_id"))
}
out <- rbindlist(list(prev, out), use.names = TRUE, fill = TRUE)
setorder(out, trait, signal_group_id, -is_representative, lead_p)
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(out, paste0(OUT, ".tmp.gz"), sep = "\t", compress = "gzip")
if (!file.rename(paste0(OUT, ".tmp.gz"), OUT)) stop("could not write ", OUT)
dir.create(dirname(STATE), showWarnings = FALSE, recursive = TRUE)
fwrite(cur, STATE, sep = "\t")   # after OUT: a run killed before here recomputes next time
unlink(tmp, recursive = TRUE)

s <- out[, .(n_cs = .N, n_signal = uniqueN(signal_group_id)), by = trait]
cat(sprintf("\ncredible sets %d -> distinct signals %d  (%.1f%% labelled duplicate)\n",
            nrow(out), sum(s$n_signal), 100*mean(!out$is_representative)))
cat(sprintf("traits with any duplicate: %d / %d (%.1f%%)\n",
            sum(s$n_cs > s$n_signal), nrow(s), 100*mean(s$n_cs > s$n_signal)))
cat("wrote ->", OUT, "\n")
