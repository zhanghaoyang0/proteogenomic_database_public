#!/usr/bin/env Rscript
library(data.table)

source(path.expand("code/lib/common.R"))   # REPO_ROOT, GWAS_BASE
repo_root = REPO_ROOT

args = commandArgs(trailingOnly = TRUE)
replication_dir = if (length(args) >= 1L) path.expand(args[1]) else file.path(repo_root, "info/bf_protein/replication")
extract_dir = if (length(args) >= 2L) path.expand(args[2]) else file.path(GWAS_BASE, "bf1/replication/extract")

candidate_file = file.path(replication_dir, "bf1_replication_candidates.tsv.gz")
request_file = file.path(replication_dir, "bf1_replication_extract_requests.tsv.gz")
candidates = fread(candidate_file)
requests = fread(request_file)

status_files = list.files(extract_dir, pattern = "^extract_[0-9]+_[0-9]+\\.status\\.tsv\\.gz$", full.names = TRUE)
data_files = list.files(extract_dir, pattern = "^extract_[0-9]+_[0-9]+\\.tsv\\.gz$", full.names = TRUE)
if (!length(status_files)) stop("No extraction status files in ", extract_dir)
status = rbindlist(lapply(status_files, fread))
if (anyDuplicated(status$request_id)) stop("Duplicated request IDs in extraction status")
if (!setequal(status$request_id, requests$request_id)) {
  stop("Extraction is incomplete: completed ", uniqueN(status$request_id),
       " of ", nrow(requests), " requests")
}
# request_id is a row number, so extracts made from an earlier request or candidate file
# would pass the completeness check above: each id must still name the same (cohort,
# assay) and the same number of requested SNPs as this run's files.
n_now = rbind(
  candidates[, .(cohort = "bf2", n_now = uniqueN(SNP)), by = .(assay = bf2_assay)],
  candidates[bf1_assay_available == TRUE, .(cohort = "bf1", n_now = uniqueN(SNP)), by = .(assay = bf1_assay)])
chk = merge(status[, .(request_id, s_cohort = cohort, s_assay = assay, requested_snps)],
            requests[, .(request_id, cohort, assay)], by = "request_id")
chk = merge(chk, n_now, by = c("cohort", "assay"), all.x = TRUE)
n_bad = chk[s_cohort != cohort | s_assay != assay | is.na(n_now) | requested_snps != n_now, .N]
if (n_bad) stop(n_bad, " extracts do not match this run's request and candidate files; rerun 02")
extracted = rbindlist(lapply(data_files, fread), use.names = TRUE, fill = TRUE)

bf2 = extracted[cohort == "bf2"]
setnames(bf2, c("assay", "EA", "NEA", "EA_FREQ", "BETA", "SE", "P", "N"),
         c("bf2_assay", "bf2_ea", "bf2_nea", "bf2_eaf", "bf2_beta", "bf2_se", "bf2_p", "bf2_n"))
bf2 = bf2[, .(bf2_assay, SNP, bf2_ea, bf2_nea, bf2_eaf, bf2_beta, bf2_se, bf2_p, bf2_n)]

bf1 = extracted[cohort == "bf1"]
setnames(bf1, c("assay", "EA", "NEA", "EA_FREQ", "BETA", "SE", "P", "N"),
         c("bf1_assay", "bf1_ea", "bf1_nea", "bf1_eaf", "bf1_beta_raw", "bf1_se", "bf1_p", "bf1_n"))
bf1 = bf1[, .(bf1_assay, SNP, bf1_ea, bf1_nea, bf1_eaf, bf1_beta_raw, bf1_se, bf1_p, bf1_n)]

# Left joins: a duplicated (assay, SNP) row in an extract would multiply candidates.
n_cand = nrow(candidates)
result = merge(candidates, bf2, by = c("bf2_assay", "SNP"), all.x = TRUE)
if (nrow(result) != n_cand) stop("BF2 extract has duplicate (assay, SNP) rows")
result = merge(result, bf1, by = c("bf1_assay", "SNP"), all.x = TRUE)
if (nrow(result) != n_cand) stop("BF1 extract has duplicate (assay, SNP) rows")

comp = function(x) chartr("ACGT", "TGCA", toupper(x))
result[, allele_alignment := fcase(
  is.na(bf1_ea) | is.na(bf2_ea), NA_character_,
  bf1_ea == bf2_ea & bf1_nea == bf2_nea, "direct",
  bf1_ea == bf2_nea & bf1_nea == bf2_ea, "swap",
  comp(bf1_ea) == bf2_ea & comp(bf1_nea) == bf2_nea, "strand_flip",
  comp(bf1_ea) == bf2_nea & comp(bf1_nea) == bf2_ea, "strand_flip_swap",
  default = "mismatch"
)]
result[, bf1_beta := fifelse(
  allele_alignment %chin% c("swap", "strand_flip_swap"),
  -bf1_beta_raw,
  bf1_beta_raw
)]
result[, bf1_eaf_aligned := fifelse(
  allele_alignment %chin% c("swap", "strand_flip_swap"),
  1 - bf1_eaf,
  bf1_eaf
)]

result[, measurement_non_testable_reason := fcase(
  !bf1_assay_available, "protein_not_measured_same_fluid",
  is.na(bf2_beta), "bf2_sentinel_absent_from_clean",
  is.na(bf1_beta_raw), "bf1_sentinel_absent_from_clean",
  allele_alignment == "mismatch", "allele_mismatch",
  !is.finite(bf1_beta) | !is.finite(bf1_se) | bf1_se <= 0 | !is.finite(bf1_p), "invalid_bf1_statistics",
  !is.finite(bf2_beta) | !is.finite(bf2_se) | bf2_se <= 0 | !is.finite(bf2_p), "invalid_bf2_statistics",
  default = NA_character_
)]
result[, measurement_testable := is.na(measurement_non_testable_reason)]

# Collapse the lead + credible-set-proxy rows down to one row per signal: the
# credible-set lead when it is testable in BF1, otherwise the highest-PIP member of
# the same credible set that is. member_rank 0 marks the lead, so ordering by
# (testable first, then rank) puts the preferred row on top; a signal with no
# testable variant at all keeps its lead row and its non-testable reason.
# .SD[1L] keeps whatever row sorts first, so check what it relies on beforehand:
# one row per (signal, member) and a lead row for every signal.
n_signals = uniqueN(candidates$replication_id)
if (anyDuplicated(result[, .(replication_id, member_rank)]))
  stop("Duplicated (replication_id, member_rank) rows before the signal collapse")
if (result[member_rank == 0L, uniqueN(replication_id)] != n_signals)
  stop("Some signals have no lead (member_rank 0) row")
setorder(result, replication_id, -measurement_testable, member_rank)
result = result[, .SD[1L], by = replication_id]
result[, proxy_used := snp_role == "proxy"]
# One BF1 association must count for one signal only: a proxy may not be shared by two
# signals of an assay, nor be another signal's tested variant (284 members are shared
# between signals in the candidate pool; none is picked twice today).
px = result[proxy_used == TRUE, .(bf2_assay, SNP)]
if (anyDuplicated(px) || nrow(result[proxy_used == FALSE][px, on = .(bf2_assay, SNP), nomatch = NULL]))
  stop("A proxy variant stands in for more than one signal of an assay")
message("Proxy fallback used for ", sum(result$proxy_used), " of ", n_signals, " signals (",
        result[strict_replication_eligible == TRUE, sum(proxy_used)], " in the strict family)")

result[, strict_non_testable_reason := fcase(
  !strict_replication_eligible & exploratory_cross_assay_eligible,
    "no_exact_assay_cross_assay_only",
  !strict_replication_eligible, "no_exact_assay",
  !measurement_testable, measurement_non_testable_reason,
  default = NA_character_
)]
result[, replication_testable := strict_replication_eligible & measurement_testable]
result[, exploratory_testable := exploratory_cross_assay_eligible & measurement_testable]
result[, direction_concordant := fifelse(
  replication_testable, sign(bf1_beta) == sign(bf2_beta), NA
)]
result[, nominal_replication := fifelse(
  replication_testable, direction_concordant & bf1_p < 0.05, NA
)]
result[, corrected_replication := fifelse(
  replication_testable,
  direction_concordant & bf1_p < bonferroni_alpha,
  NA
)]
result[, bf1_fdr := NA_real_]
result[replication_testable == TRUE, bf1_fdr := p.adjust(bf1_p, method = "BH", n = correction_family_n[1])]
result[, fdr_replication := fifelse(
  replication_testable, direction_concordant & bf1_fdr < 0.05, NA
)]
result[, exploratory_nominal_support := fifelse(
  exploratory_testable,
  sign(bf1_beta) == sign(bf2_beta) & bf1_p < 0.05,
  NA
)]
result[, heterogeneity_z := fifelse(
  replication_testable,
  (bf2_beta - bf1_beta) / sqrt(bf2_se^2 + bf1_se^2),
  NA_real_
)]
result[, heterogeneity_p := 2 * pnorm(-abs(heterogeneity_z))]

setorder(result, fluid, bf2_assay, bf2_lead_p, SNP)
preferred = c(
  "replication_id", "fluid", "gene", "ensembl_gene_id", "bf2_assay", "bf1_assay",
  "assay_match_type", "SNP", "is_cis", "bf2_ea", "bf2_nea", "bf2_beta", "bf2_se",
  "bf2_p", "bf2_n", "bf1_ea", "bf1_nea", "bf1_beta", "bf1_se", "bf1_p", "bf1_n",
  "allele_alignment", "direction_concordant", "replication_testable",
  "nominal_replication", "corrected_replication", "bf1_fdr", "fdr_replication",
  "heterogeneity_p", "strict_non_testable_reason",
  "measurement_non_testable_reason", "exploratory_testable",
  "exploratory_nominal_support"
)
setcolorder(result, c(preferred, setdiff(names(result), preferred)))

summarize = function(x, label_cols = character()) {
  x[, .(
    bf2_signals = .N,
    strict_eligible = sum(strict_replication_eligible),
    testable = sum(replication_testable),
    testable_via_proxy = sum(replication_testable & proxy_used),
    direction_concordant = sum(direction_concordant, na.rm = TRUE),
    nominal_replication = sum(nominal_replication, na.rm = TRUE),
    bonferroni_replication = sum(corrected_replication, na.rm = TRUE),
    fdr_replication = sum(fdr_replication, na.rm = TRUE),
    testable_fraction = if (sum(strict_replication_eligible) > 0L) {
      sum(replication_testable) / sum(strict_replication_eligible)
    } else NA_real_,
    direction_concordance = mean(direction_concordant, na.rm = TRUE),
    nominal_rate_testable = mean(nominal_replication, na.rm = TRUE),
    bonferroni_rate_testable = mean(corrected_replication, na.rm = TRUE),
    nominal_rate_planned = if (sum(strict_replication_eligible) > 0L) {
      sum(nominal_replication, na.rm = TRUE) / sum(strict_replication_eligible)
    } else NA_real_,
    bonferroni_rate_planned = if (sum(strict_replication_eligible) > 0L) {
      sum(corrected_replication, na.rm = TRUE) / sum(strict_replication_eligible)
    } else NA_real_,
    # Testable rows only: other rows can carry an unaligned or invalid BF1 beta.
    beta_pearson = if (sum(replication_testable) > 1L)
      cor(bf2_beta[replication_testable], bf1_beta[replication_testable]) else NA_real_,
    beta_spearman = if (sum(replication_testable) > 1L)
      cor(bf2_beta[replication_testable], bf1_beta[replication_testable], method = "spearman") else NA_real_
  ), by = label_cols]
}

summary = rbindlist(list(
  cbind(stratum = "overall", summarize(result)),
  cbind(stratum = "fluid", summarize(result, "fluid")),
  cbind(stratum = "cis_trans", summarize(result, "is_cis")),
  cbind(stratum = "fluid_cis_trans", summarize(result, c("fluid", "is_cis"))),
  cbind(stratum = "match_type", summarize(result, "assay_match_type"))
), use.names = TRUE, fill = TRUE)

fwrite(result, file.path(replication_dir, "bf1_targeted_replication.tsv.gz"), sep = "\t")
fwrite(summary, file.path(replication_dir, "bf1_targeted_replication_summary.tsv"), sep = "\t")
fwrite(result[, .N, by = strict_non_testable_reason][order(-N)],
       file.path(replication_dir, "bf1_targeted_replication_non_testable.tsv"), sep = "\t")
print(summary[stratum == "overall"])
message("Wrote BF1 targeted replication results for ", nrow(result), " BF2 signals")
