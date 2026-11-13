#!/usr/bin/env Rscript
# Build the two app-facing tables that make a discovery locus the first-class
# object and a credible set its evidence:
#
#   info/bf_protein/signals.tsv.gz        one row per discovery locus (clump lead)
#   info/bf_protein/credible_sets.tsv.gz  one row per credible set, joined to its locus
#
# Why: the app used to key everything on cs_id, so a clump lead whose SuSiE window
# produced no credible set had no row, no plot and no explanation -- 1,700 of the
# 12,290 loci. It also reported the CS lead (max PIP) in one panel and the clump
# lead (min P) in another without saying they were the same signal.
#
# Inputs (all existing, nothing is recomputed here):
#   info/bf_protein/clump_lead_pqtl.csv.gz  authoritative locus list, one row per
#                                           (cohort, protein_assay, SNP)
#   gwas/bf2/finemapping/protein/<trait>.finemapped.tsv.gz       CS members + PIP
#   info/bf_protein/cs_dedup_map.tsv.gz     cross-window dedup (see 04_dedup_*.R)
#
# The clump table is the authority on which loci exist. The identity
#   (loci with >=1 CS) + (loci with none) == (distinct clump leads)
# must hold exactly; we stop() if it does not rather than silently intersecting.

suppressPackageStartupMessages(library(data.table))
source(path.expand("code/lib/common.R"))   # FM_WINDOW_BP

INFO      <- path.expand("./info/bf_protein")
GWAS_ROOT <- path.expand(Sys.getenv("BF_GWAS_ROOT", GWAS_BASE))
WINDOW    <- as.integer(Sys.getenv("BF_FM_WINDOW", FM_WINDOW_BP))  # 01_finemapping.R uses FM_WINDOW_BP
COHORTS   <- "bf2"   # BF1 is not fine-mapped in the live pipeline (docs/project_status.md)

CLUMP  <- file.path(INFO, "clump_lead_pqtl.csv.gz")
DEDUP  <- file.path(INFO, "cs_dedup_map.tsv.gz")
OUT_S  <- file.path(INFO, "signals.tsv.gz")
OUT_CS <- file.path(INFO, "credible_sets.tsv.gz")

assert_file <- function(p) if (!file.exists(p)) stop("Missing input: ", p)
assert_file(CLUMP); assert_file(DEDUP)

# ---- 1. loci from the clump table -------------------------------------------
cl <- fread(CLUMP)
need <- c("cohort", "protein_assay", "SNP", "CHR", "BP", "P", "is_cis", "in_mhc",
          "gene", "ensembl_gene_id")
if (!all(need %in% names(cl))) stop("clump_lead_pqtl.csv.gz missing: ",
                                    paste(setdiff(need, names(cl)), collapse = ", "))
sig <- cl[, .(trait = paste0(cohort, "__", protein_assay),
              cohort, assay = protein_assay,
              signal_id = SNP, clump_lead_snp = SNP,
              chr = CHR, pos = BP, pval = P, is_cis = as.logical(is_cis),
              in_mhc = as.logical(in_mhc),
              # The assay's own gene, so the locus plot can highlight it in the
              # gene track. Carried as ENSG too: 144 Olink assays are named for a
              # legacy symbol (U_PAR, TNF_R2) that no longer matches GENCODE, so
              # matching the track by symbol would miss them.
              gene = gene, ensembl_gene_id = ensembl_gene_id)]
if (anyDuplicated(sig[, .(trait, signal_id)]))
  stop("clump table has duplicate (trait, lead) rows")
# signal_id ends up in filenames, so it must stay filename-safe.
bad <- sig[!grepl("^[A-Za-z0-9._:-]+$", signal_id)]
if (nrow(bad)) stop(nrow(bad), " clump leads are not filename-safe, e.g. ", bad$signal_id[1])

sig[, `:=`(window_start = pmax(1L, pos - WINDOW), window_end = pos + WINDOW)]

# ---- 2. credible sets from the finemapped files ------------------------------
fm_files <- unlist(lapply(COHORTS, function(co)
  list.files(file.path(GWAS_ROOT, co, "finemapping/protein"),
             pattern = "\\.finemapped\\.tsv\\.gz$", full.names = TRUE)))
if (!length(fm_files)) stop("No .finemapped.tsv.gz under ", GWAS_ROOT)
message("reading ", length(fm_files), " finemapped files ...")

fm <- rbindlist(lapply(fm_files, function(f) {
  x <- tryCatch(fread(f, select = c("trait", "clump_lead_snp", "lead_snp",
                                    "cs_id", "cs_purity", "is_cis", "chr", "pos", "snp", "pip")),
                error = function(e) stop("unreadable: ", f, " -- ", conditionMessage(e)))
  if (!nrow(x)) return(NULL)
  x[, .(cs_size = .N,
        top_pip_snp = snp[which.max(pip)],
        top_pip = max(pip),
        lead_snp = lead_snp[1L],
        cs_purity = cs_purity[1L],
        is_cis = as.logical(is_cis[1L]),    # per CS, set by 01_finemapping.R
        top_chr = chr[which.max(pip)], top_pos = pos[which.max(pip)]),
    by = .(trait, clump_lead_snp, cs_id)]
}), use.names = TRUE)
if (anyNA(fm$is_cis)) stop(sum(is.na(fm$is_cis)), " credible sets have no is_cis")
# 01 decides is_cis from the cis windows (top-PIP variant inside one of the assay's
# windows), but its freshness check cannot see the windows file: it is tracked in git, so
# its mtime moves on every checkout. Recompute here and stop on a fit that predates a
# window change.
source(file.path(REPO_ROOT, "code/lib/cis_windows.R"))
fm[, cis_now := in_cis_window(sub("^bf[12]__", "", trait), top_chr, top_pos, read_cis_windows(REPO_ROOT))]
stale_cis <- fm[cis_now != is_cis, unique(trait)]
if (length(stale_cis))
  stop(length(stale_cis), " traits were fine-mapped with older cis windows (is_cis differs), e.g. ",
       paste(head(stale_cis, 3), collapse = ", "), " -- refit them with FORCE=1")
fm[, c("top_chr", "top_pos", "cis_now") := NULL]
# A file with rows but no CS is impossible; a file with zero rows is the
# "fine-mapped, no credible set" case and legitimately contributes nothing here.
finemapped_traits <- sub("\\.finemapped\\.tsv\\.gz$", "", basename(fm_files))

# Per-locus sidecars distinguish SuSiE failures from genuine runs with no CS.
status_files <- unlist(lapply(COHORTS, function(co)
  list.files(file.path(GWAS_ROOT, co, "finemapping/protein"),
             pattern = "\\.finemap_status\\.tsv\\.gz$", full.names = TRUE)))
fm_status <- rbindlist(lapply(status_files, function(f) fread(
  f, select = c("trait", "clump_lead_snp", "finemap_status", "detail"))), fill = TRUE)
if (nrow(fm_status) && anyDuplicated(fm_status[, .(trait, clump_lead_snp)]))
  stop("duplicate (trait, clump_lead_snp) rows in fine-mapping status sidecars")
# 01 writes this detail when the window's credible set is the one another window found.
DUP_DETAIL <- "credible set duplicated by another clump window"
if (nrow(fm_status)) fm_status[detail %chin% DUP_DETAIL, finemap_status := "cs_in_other_window"]
status_lookup <- if (nrow(fm_status)) {
  setNames(fm_status$finemap_status,
           paste(fm_status$trait, fm_status$clump_lead_snp, sep = "\r"))
} else setNames(character(), character())

# ---- 3. dedup facts ----------------------------------------------------------
dd <- fread(DEDUP)
need <- c("trait", "signal_group_id", "cs_id", "clump_lead_snp",
          "is_representative", "duplicate_of_cs_id")
if (!all(need %in% names(dd))) stop("cs_dedup_map.tsv.gz missing: ",
                                    paste(setdiff(need, names(dd)), collapse = ", "))
setnames(dd, "clump_lead_snp", "window")   # local alias; the finemapped file uses both

cs <- merge(fm, dd[, .(trait, cs_id, window, dd_lead = lead_snp, signal_group_id,
                       is_representative, duplicate_of_cs_id, max_r2)],
            by = c("trait", "cs_id"), all.x = TRUE)
if (anyNA(cs$signal_group_id))
  stop(sum(is.na(cs$signal_group_id)), " credible sets are absent from ", basename(DEDUP),
       " -- re-run 03_dedup_credible_sets.R")
# Row by row: cs_id is renumbered on every fine-mapping run, so a dedup map older
# than a partial rerun can join each key to a different credible set.
n_bad <- cs[window != clump_lead_snp | dd_lead != lead_snp, .N]
if (n_bad) stop(n_bad, " credible sets disagree with ", basename(DEDUP),
                " on window or lead SNP -- re-run 03_dedup_credible_sets.R")
cs[, c("window", "dd_lead") := NULL]
setnames(cs, "clump_lead_snp", "signal_id")

# ---- 4. per-locus rollup -----------------------------------------------------
# A locus with credible sets is cis when any of them is (cis is decided per credible
# set in 01_finemapping.R), so it agrees with coloc, which uses the cis sets. Taking
# only the strongest set left CBLN2 with 4 cis sets and no cis locus. A locus without
# credible sets keeps the clump lead's status.
roll <- cs[, .(n_cs = .N, n_cs_rep = sum(is_representative %in% TRUE), cs_is_cis = any(is_cis)),
           by = .(trait, signal_id)]
sig <- merge(sig, roll, by = c("trait", "signal_id"), all.x = TRUE)
sig[is.na(n_cs), `:=`(n_cs = 0L, n_cs_rep = 0L)]
sig[n_cs > 0L, is_cis := cs_is_cis]
sig[, cs_is_cis := NULL]
# Sharing is a property of a credible set, not of a discovery locus: a locus whose
# CS was grouped with a neighbour's is NOT "the same signal" as that neighbour, it
# just has one credible set that fine-maps to the same effect. The UI says
# "2 CS - 1 shared" off this count and puts the "shared with <locus>" badge on the
# credible set itself; it never folds or relabels the discovery signal.
sig[, n_cs_shared := n_cs - n_cs_rep]
sig[, status_key := paste(trait, signal_id, sep = "\r")]
sig[, recorded_status := unname(status_lookup[status_key])]

# MHC loci are skipped by 01_finemapping.R (../lib/mhc.R), so they get their own
# status rather than reading as "SuSiE ran and found nothing".
sig[, finemap_status := fcase(n_cs > 0L, "credible_set",
                              in_mhc, "mhc",
                              trait %in% finemapped_traits, "no_credible_set",
                              default = "not_finemapped")]
sig[n_cs == 0L & !in_mhc & recorded_status %chin%
      c("susie_error", "insufficient_data", "no_credible_set", "cs_in_other_window"),
    finemap_status := recorded_status]
# A trait with a sidecar lists every locus it fine-mapped; a locus missing from it was
# added by a re-clump after the fit, and is not "no credible set" until it is refit.
side_traits <- unique(fm_status$trait)
n_new <- sig[n_cs == 0L & !in_mhc & trait %chin% side_traits & is.na(recorded_status), .N]
if (n_new) stop(n_new, " loci are missing from their trait's status sidecar (added after the fit) -- refit those traits")
if (sig[n_cs == 0L & recorded_status == "credible_set", .N])
  stop("status sidecar says credible_set for a locus absent from the fine-mapping table")
sig[, c("status_key", "recorded_status") := NULL]
if (sig[in_mhc & n_cs > 0L, .N]) stop("MHC loci carry credible sets -- fine-mapped before the MHC rule")
# Every assay with a cis lead is meant to be fine-mapped (submit_01_finemapping_gwwide.sh),
# so a cis trait left "not_finemapped" in a fine-mapped cohort is a failed or
# skipped run, not a category.
fm_cohorts <- unique(sub("__.*", "", finemapped_traits))
miss <- sig[is_cis %in% TRUE & cohort %in% fm_cohorts & finemap_status == "not_finemapped", unique(trait)]
if (length(miss)) stop(length(miss), " cis traits have no .finemapped file, e.g. ",
                       paste(head(miss, 3), collapse = ", "), " -- re-run submit_01_finemapping_gwwide.sh")

# INTERNAL METADATA ONLY -- never shown to a user and never used to hide a row.
# TRUE when the locus contributes at least one credible set that is not shared
# with another window. Useful for counting; the picker still lists all 12,290
# discovery signals. Deliberately NOT connected components: a locus can carry one
# group's shared CS and another group's representative, and component merging
# would swallow the signal that exists only there (CFHR4 in CSF would collapse
# 10 loci holding 9 representative CS into one).
sig[, is_primary := (n_cs == 0L) | (n_cs_rep > 0L)]

# Per credible set, the locus holding the representative of its group -- the
# "Shared with <locus>" badge. Empty when this CS is itself the representative.
rep_of_group <- unique(cs[is_representative %in% TRUE, .(trait, signal_group_id, rep_signal = signal_id)])
cs <- merge(cs, rep_of_group, by = c("trait", "signal_group_id"), all.x = TRUE)
cs[, shared_with_signal := fifelse(is_representative %in% TRUE | is.na(rep_signal) |
                                     rep_signal == signal_id, "", rep_signal)]
cs[, rep_signal := NULL]

# ---- 5. display index and plot basenames -------------------------------------
setorder(sig, trait, chr, pos)
sig[, signal_index := seq_len(.N), by = trait]          # display only, never a key
sig[, locus_plot := paste0(trait, "__", signal_id, "__locus.png")]

setorder(cs, trait, signal_id, cs_id)
cs[, cs_index := seq_len(.N), by = .(trait, signal_id)]  # display only
cs[, plot_path := paste0(trait, "__", signal_id, "__cs", cs_id, ".png")]

# ---- 6. invariants -----------------------------------------------------------
if (nrow(cs[!paste(trait, signal_id) %in% sig[, paste(trait, signal_id)]]))
  stop("credible sets reference a locus absent from the clump table")
if (sum(sig$n_cs) != nrow(cs))
  stop("sum(n_cs) = ", sum(sig$n_cs), " != nrow(credible_sets) = ", nrow(cs))
if (uniqueN(cs[, paste(trait, cs_id)]) != nrow(cs))
  stop("cs_id is not unique within a trait")
span <- sig[window_start > 1L, unique(window_end - window_start)]
if (length(span) != 1L || span != 2L * WINDOW)
  stop("window width is not 2 x ", WINDOW)

setcolorder(sig, c("trait", "cohort", "assay", "signal_id", "signal_index",
                   "clump_lead_snp", "chr", "pos", "pval", "is_cis", "in_mhc",
                   "gene", "ensembl_gene_id",
                   "window_start", "window_end", "n_cs", "n_cs_rep", "n_cs_shared",
                   "finemap_status", "locus_plot", "is_primary"))
setcolorder(cs, c("trait", "signal_id", "cs_id", "cs_index", "lead_snp", "top_pip_snp",
                  "top_pip", "cs_size", "cs_purity", "is_cis", "shared_with_signal",
                  "signal_group_id", "is_representative", "duplicate_of_cs_id",
                  "max_r2", "plot_path"))

fwrite(sig, OUT_S,  sep = "\t", compress = "gzip")
fwrite(cs,  OUT_CS, sep = "\t", compress = "gzip")

st <- sig[, .N, by = finemap_status][order(-N)]
message(sprintf("signals: %d loci (%d traits) | %s -> %s",
                nrow(sig), uniqueN(sig$trait),
                paste(sprintf("%s %d", st$finemap_status, st$N), collapse = ", "), OUT_S))
message(sprintf("credible sets: %d observed -> %d representative, %d shared -> %s",
                nrow(cs), sum(cs$is_representative %in% TRUE),
                sum(nzchar(cs$shared_with_signal)), OUT_CS))
