#!/usr/bin/env Rscript
# One verdict per credible set, computed once.
#
# 03 writes the raw replication table: one row per submitted sentinel, with the
# ingredients of a verdict spread over four columns (assay_match_type,
# replication_testable, exploratory_testable, the two betas and bf1_p). Every
# consumer then had to combine them itself, and three did -- the manuscript's
# 01a, 01b and 06 each rebuilt the same fcase. Add the Shiny app as a fourth and
# the repo carries four definitions of "replicated". This script is the single
# owner: it applies the rule, writes the verdict as a column, and asserts that
# the funnel it produces is still the published one.
#
# The rule is the manuscript's, not the prespecified Bonferroni family that 03
# reports. Exact-assay and cross-assay matches both count as testable, and
# replication is a concordant direction plus BF1 P < 0.05 with no multiplicity
# correction: the sentinels are pre-specified from BF2 rather than rediscovered,
# so each is a single planned test. 03's corrected_replication is carried
# through unchanged so the strict analysis stays inspectable.
#
# Usage: Rscript 04_credible_set_status.R [replication_dir] [info_dir]

library(data.table)

source(path.expand("code/lib/common.R"))   # REPO_ROOT, GWAS_BASE
repo_root = REPO_ROOT

args = commandArgs(trailingOnly = TRUE)
info_dir = if (length(args) >= 2L) path.expand(args[2]) else file.path(repo_root, "info/bf_protein")
replication_dir = if (length(args) >= 1L) path.expand(args[1]) else file.path(info_dir, "replication")

repl_file = file.path(replication_dir, "bf1_targeted_replication.tsv.gz")
dedup_file = file.path(info_dir, "cs_dedup_map.tsv.gz")
out_file = file.path(replication_dir, "cs_replication.tsv.gz")
for (f in c(repl_file, dedup_file)) {
  if (!file.exists(f)) stop("missing input: ", f)
}

repl = fread(repl_file, select = c(
  "cohort", "bf2_assay", "SNP", "cs_id", "assay_match_type", "replication_testable",
  "exploratory_testable", "bf1_assay", "bf1_beta", "bf1_se", "bf1_p", "bf1_n",
  "bf2_beta", "corrected_replication", "proxy_used", "snp_role",
  "measurement_non_testable_reason", "nominal_replication", "exploratory_nominal_support"
))
# The app and the info tables key on trait, not on (cohort, assay).
repl[, trait := paste(cohort, bf2_assay, sep = "__")]

# ---- The verdict -------------------------------------------------------------
repl[, testable := (assay_match_type == "exact_assay" & replication_testable == TRUE) |
                   (assay_match_type == "same_gene_fluid" & exploratory_testable == TRUE)]
# Concordance only where the test was made: elsewhere a BF1 beta can be unaligned
# (allele mismatch) or missing.
repl[, direction_concordant := fifelse(testable, sign(bf1_beta) == sign(bf2_beta), NA)]
# "Replicated" is 03's verdict, not a second copy of its rule (same direction and BF1
# P < 0.05): nominal_replication for exact-assay rows, exploratory_nominal_support for
# cross-assay rows.
repl[, replicated := fifelse(assay_match_type == "exact_assay", nominal_replication,
                             exploratory_nominal_support) %in% TRUE]
# One label per 03 reason, so S1 and the app name the actual one.
repl[, bf1_status := fcase(
  testable & replicated, "replicated",
  testable, "tested, not replicated",
  assay_match_type == "not_measured_same_fluid", "not testable: protein not measured in BF1",
  measurement_non_testable_reason == "bf1_sentinel_absent_from_clean", "not testable: sentinel absent from BF1 sumstats",
  measurement_non_testable_reason == "bf2_sentinel_absent_from_clean", "not testable: sentinel absent from BF2 sumstats",
  measurement_non_testable_reason == "allele_mismatch", "not testable: alleles do not match in BF1",
  measurement_non_testable_reason %chin% c("invalid_bf1_statistics", "invalid_bf2_statistics"),
    "not testable: invalid summary statistics",
  default = NA_character_
)]
if (anyNA(repl$bf1_status))
  stop(repl[is.na(bf1_status), .N], " sentinels have no bf1_status; a new 03 reason needs a label here")

# ---- Expand to credible-set grain --------------------------------------------
# 03 submitted one row per sentinel VARIANT, and a sentinel can be the lead of
# two credible sets for the same protein; only one of the two carries the row's
# cs_id. Joining on the row's own cs_id alone would leave the other set
# unlabelled although its lead was tested, so the verdict is spread across every
# credible set that shares the lead.
#
# But `SNP` is NOT reliably that lead: when the sentinel is absent from the BF1
# sumstats, 03 substitutes the highest-PIP credible-set member and writes the
# PROXY into SNP. So the lead is recovered from the dedup map through cs_id
# first, and only then used to expand. Joining straight on SNP silently dropped
# every proxy-rescued sentinel.
map = fread(dedup_file, select = c("trait", "cs_id", "lead_snp", "is_representative"))
map = map[startsWith(trait, "bf2__")]

verdict = repl[, .(trait, cs_id, bf1_status, bf1_assay,
                   bf1_assay_match = assay_match_type, bf1_beta, bf1_se, bf1_p,
                   bf1_n, direction_concordant, corrected_replication,
                   proxy_used, tested_variant = SNP)]
# cs_id -> the sentinel's own lead variant, whatever variant was looked up.
# Fine-mapping renumbers cs_id on every run, so a rerun of it after replication 01-03
# would attach these verdicts to other sets: every row must find its set, and a row
# that tested the lead itself must find that same lead there.
n_verdict = nrow(verdict)
verdict = merge(verdict, map[, .(trait, cs_id, cs_lead_snp = lead_snp)],
                by = c("trait", "cs_id"))
if (nrow(verdict) != n_verdict)
  stop(n_verdict - nrow(verdict), " sentinels have no credible set under their cs_id; rerun replication 01-03")
n_moved = verdict[proxy_used == FALSE & tested_variant != cs_lead_snp, .N]
if (n_moved) stop(n_moved, " sentinels' cs_id now names another lead; rerun replication 01-03")
verdict[, cs_id := NULL]
if (verdict[, uniqueN(bf1_status), by = .(trait, cs_lead_snp)][, any(V1 > 1L)])
  stop("one credible-set lead has conflicting BF1 verdicts")
verdict = unique(verdict, by = c("trait", "cs_lead_snp"))

cs = merge(map, verdict, by.x = c("trait", "lead_snp"), by.y = c("trait", "cs_lead_snp"))
setnames(cs, "lead_snp", "cs_lead_snp")
setcolorder(cs, c("trait", "cs_id", "cs_lead_snp", "is_representative", "bf1_status"))
setorder(cs, trait, cs_id)

# A credible set with no row here was never submitted -- it holds no genome-wide
# significant variant, so it is absent rather than untestable. Consumers assign
# that level on the left join; writing it here would claim a test that was never
# requested.

# ---- The published funnel must not move --------------------------------------
# Table 1 of the manuscript, on the deduplicated sentinel set. These numbers are
# in print, so a drift here is a silent correction to a published result.
sent = unique(cs[is_representative == TRUE],
              by = c("trait", "cs_lead_snp"))
funnel = c(
  submitted = nrow(sent),
  measured_in_bf1 = sent[, sum(bf1_assay_match != "not_measured_same_fluid")],
  testable = sent[, sum(startsWith(bf1_status, "replicated") |
                        startsWith(bf1_status, "tested"))],
  replicated = sent[, sum(bf1_status == "replicated")]
)
#
# Updated 2026-09-22, when 01a_extract_cs_members.R was finally run and the proxy
# rescue it feeds became active for the first time. The denominators did not move
# -- submitted and measured_in_bf1 are identical -- and the non-proxy subset still
# reproduces the published figures exactly:
#
#   non-proxy       3,540 sentinels   2,807 testable   2,434 replicated  <- as printed
#   proxy-rescued     122 sentinels     122 testable      89 replicated
#   total           3,662 sentinels   2,929 testable   2,523 replicated
#
# So Table 1, Fig. 1c,d and S1 need reprinting, not rechecking: 86.7% becomes
# 86.1% because a proxy is a weaker test than its own sentinel, not because
# anything about the original analysis changed.
#
# Updated 2026-10-02 after the genotype-QC rerun (dosage missingness; BF2 7.70M
# variants, BF2 genome-wide fine-mapping, SuSiE max_iter 1000): 3,662 -> 3,985
# submitted, 2,929 -> 3,417 testable, 2,523 -> 2,876 replicated (86.1% -> 84.2%).
# Updated 2026-10-04 after dropping the high-LD exclusion from clump (MHC loci are
# counted but not fine-mapped, so not submitted): 3,985 -> 4,032 submitted,
# 3,417 -> 3,459 testable, 2,876 -> 2,904 replicated (84.2% -> 84.0%).
# The manuscript reads these numbers live (since 2026-10-16).
# Updated 2026-10-12 after the genetic-key / BF1 sid / sex-check rerun: 4,032 -> 4,033
# submitted, 3,577 -> 3,570 measured, 3,459 -> 3,455 testable, 2,904 -> 2,920
# replicated (84.0% -> 84.5%). Manuscript edits: issues/resolved/geno_key_rerun_20261009/README.md.
# Updated 2026-11-01 when the proxy pool moved to the reported credible set (code review
# 2026-11-01, C52): 3,455 -> 3,456 testable (CTSB CSF, cross-assay, not replicated),
# 2,920 replicated (84.5%).
# Updated 2026-11-04 after the BF1 redo (imputed from the raw array, 1,497 GWAS samples;
# BF1 protein sids from the genetic key, issues/bf1_protein_sid_spelling_20261102, so
# most BF1 assays gained several hundred samples; three new BF1 panels; request r029):
# 3,570 -> 3,588 measured, 3,456 -> 3,536 testable, 2,920 -> 3,110 replicated
# (84.5% -> 88.0%).
expected = c(submitted = 4033L, measured_in_bf1 = 3588L, testable = 3536L, replicated = 3110L)
if (!identical(as.integer(funnel), as.integer(expected))) {
  stop("replication funnel moved -- these numbers are published in Table 1\n",
       "  got:      ", paste(names(funnel), funnel, sep = "=", collapse = " "), "\n",
       "  expected: ", paste(names(expected), expected, sep = "=", collapse = " "), "\n",
       "  If the change is intended, update `expected` here and the manuscript together.")
}

fwrite(cs, out_file, sep = "\t")
print(cs[, .N, by = bf1_status][order(-N)])
message("Wrote ", out_file, ": ", nrow(cs), " credible sets with a BF1 verdict (",
        sum(cs$is_representative), " representative), from ", nrow(repl), " submitted sentinels")
message("Funnel (Table 1): ", paste(names(funnel), funnel, sep = " = ", collapse = ", "))
