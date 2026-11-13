# Credible-set -> distinct-signal map, shared by every script that counts signals.
#
# Overlapping +/-500 kb fine-mapping windows resolve one signal more than once, so a
# raw credible-set count over-counts by ~15%. code/04_finemapping/03_dedup_credible_sets.R
# writes the map that says which credible set represents a signal; this is the one
# place the report joins it, so no script can quietly count raw rows again.
#
# Usage (d must carry cohort, protein_assay, cs_id):
#   source("code/09_statistical_analysis/lib/cs_signals.R")
#   fm <- distinct_signals(fm)   # drop duplicate credible sets
#   fm <- add_cs_signal(fm)      # keep all rows, add is_representative
#
# repo_root defaults to the working directory (the repo root), matching lib/eqtl_coloc_scope.R.
#
# Threshold: the r2 > 0.1 map only (no r2 0.5 sensitivity run; author decision 2026-11-01).
#
# Reporting threshold: fine-mapping keeps only credible sets holding a variant at
# P < 5e-8, which is the cis threshold. A trans signal is reported only when its lead
# variant reaches the study-wide 5e-8 / (number of analysed assays), as in Western
# 2024 (Nat Genet) and Sun 2023 (Nature). Between the two thresholds BF2 trans signals
# replicated in BF1 at ~35% against ~80% for cis signals of the same P (2026-10-24).
# distinct_signals() returns reported signals only; add_cs_signal() keeps every row and
# adds is_reported. Lookups that ask whether a signal is seen elsewhere (sharing with
# the other fluid, replication) still use every credible set, as replication does.
N_ANALYSED_ASSAYS <- 3538L   # 1,406 CSF + 2,132 plasma, 00_kept_assays.csv
TRANS_P <- 5e-8 / N_ANALYSED_ASSAYS
KEPT_ASSAYS_CSV <- "output/00_kept_assays.csv"
#
# Requires data.table, loaded by the caller.

CS_DEDUP_MAP <- function(repo_root = ".") file.path(repo_root, "info", "bf_protein", "cs_dedup_map.tsv.gz")

add_cs_signal <- function(d, repo_root = ".") {
  f <- CS_DEDUP_MAP(repo_root)
  if (!file.exists(f)) {
    stop("missing credible-set signal map: ", f,
         "\n  run: Rscript code/04_finemapping/03_dedup_credible_sets.R")
  }
  for (col in c("cohort", "protein_assay", "cs_id")) {
    if (!col %in% names(d)) stop("add_cs_signal() needs a ", col, " column")
  }
  if ("is_reported" %in% names(d)) stop("add_cs_signal(): d already has is_reported")
  # The trans threshold divides by the analysed assay count; a dedup rerun that
  # changes it must change the constant too.
  kept_f <- file.path(repo_root, KEPT_ASSAYS_CSV)
  if (!file.exists(kept_f)) stop("missing ", kept_f, " (needed for the trans threshold)")
  n_kept <- nrow(fread(kept_f, select = "protein_assay"))
  if (n_kept != N_ANALYSED_ASSAYS) {
    stop("00_kept_assays.csv has ", n_kept, " assays but N_ANALYSED_ASSAYS is ",
         N_ANALYSED_ASSAYS, " -- update it in lib/cs_signals.R")
  }
  # The map's grouping id is signal_group_id, not signal_id: a "signal" in
  # info/bf_protein/signals.tsv.gz is a DISCOVERY LOCUS, a different partition
  # (12,290 loci vs 12,119 dedup groups; 872 groups span more than one locus).
  # This helper counts dedup groups and must keep doing so -- do not repoint it
  # at signals.tsv.gz, it would silently change a published count.
  map <- fread(f, select = c("trait", "cs_id", "signal_group_id",
                             "is_representative", "duplicate_of_cs_id"))
  d <- copy(d)
  d[, trait_key := paste(cohort, protein_assay, sep = "__")]
  out <- merge(d, map, by.x = c("trait_key", "cs_id"), by.y = c("trait", "cs_id"),
               all.x = TRUE)
  # A credible set the map has never seen means the map is stale relative to the
  # fine-mapping output. Silently treating it as unique would restore the very
  # over-count this file exists to remove, so stop instead.
  n_bad <- sum(is.na(out$is_representative))
  if (n_bad > 0) {
    stop(n_bad, " of ", nrow(out), " credible sets are absent from ", f,
         "\n  rerun code/04_finemapping/03_dedup_credible_sets.R")
  }
  # cs_id is renumbered on every fine-mapping run, so a map built before a rerun can
  # still match on (trait, cs_id) while pointing at a different set. Its set leads
  # must agree with the current lead table, which the same run writes.
  cur <- fread(file.path(repo_root, "info", "bf_protein", "finemap_lead_snps.csv"),
               select = c("cohort", "protein_assay", "cs_id", "snp"))
  cur[, trait := paste(cohort, protein_assay, sep = "__")]
  leads <- merge(cur[, .(trait, cs_id, snp)],
                 fread(f, select = c("trait", "cs_id", "lead_snp")), by = c("trait", "cs_id"))
  n_moved <- sum(!((leads$snp == leads$lead_snp) %in% TRUE))   # a missing lead counts as moved
  if (n_moved > 0) {
    stop(n_moved, " credible sets in ", f, " have a different lead SNP than ",
         "finemap_lead_snps.csv -- the map predates the last fine-mapping run",
         "\n  rerun code/04_finemapping/03_dedup_credible_sets.R")
  }
  rep <- reported_sets(repo_root)
  out <- merge(out, rep[, .(trait_key = trait, cs_id, is_reported)],
               by = c("trait_key", "cs_id"), all.x = TRUE)
  if (anyNA(out$is_reported)) stop("credible sets without a cis/trans call in finemap_lead_snps.csv")
  out[, trait_key := NULL]
  out[]
}

# One row per credible set: is it a reported signal (cis, or trans with lead P < TRANS_P)?
reported_sets <- function(repo_root = ".") {
  map <- fread(CS_DEDUP_MAP(repo_root), select = c("trait", "cs_id", "lead_snp", "lead_p"))
  cur <- fread(file.path(repo_root, "info", "bf_protein", "finemap_lead_snps.csv"),
               select = c("cohort", "protein_assay", "cs_id", "is_cis"))
  cur <- unique(cur[, .(trait = paste(cohort, protein_assay, sep = "__"), cs_id,
                        is_cis = as.logical(is_cis))])
  r <- merge(map, cur, by = c("trait", "cs_id"), all.x = TRUE)
  if (anyNA(r$is_cis)) stop(sum(is.na(r$is_cis)), " credible sets in the dedup map have no is_cis")
  r[, is_reported := is_cis | lead_p < TRANS_P]
  r[]
}

# (trait, cs_lead_snp) of the reported distinct signals, for tables keyed on the lead
# variant (the replication tables).
reported_leads <- function(repo_root = ".") {
  map <- fread(CS_DEDUP_MAP(repo_root), select = c("trait", "cs_id", "is_representative"))
  r <- merge(reported_sets(repo_root), map, by = c("trait", "cs_id"))
  unique(r[is_representative == TRUE & is_reported == TRUE, .(trait, cs_lead_snp = lead_snp)])
}

distinct_signals <- function(d, repo_root = ".") {
  out <- add_cs_signal(d, repo_root)
  n_rep <- sum(out$is_representative)
  n_out <- sum(out$is_representative & out$is_reported)
  message("credible sets: ", nrow(out), " -> ", n_rep, " distinct (dropped ",
          nrow(out) - n_rep, " duplicates; map ", basename(CS_DEDUP_MAP(repo_root)),
          ") -> ", n_out, " reported (dropped ", n_rep - n_out,
          " trans signals with lead P >= ", signif(TRANS_P, 3), ")")
  out[is_representative == TRUE & is_reported == TRUE]
}
