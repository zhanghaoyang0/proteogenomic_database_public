#!/usr/bin/env Rscript
# Summarize BF2 protein-disease evidence across FinnGen disease categories.
# Analysis only: this script reads the secure MR and coloc result tables and writes
# tidy summary tables; 05c_plot.R draws from those tables and never touches
# the secure inputs.
# Evidence layers:
#   1. MR: global BH-FDR < 0.05 among successful protein -> disease tests
#   2. Coloc: maximum PP4 > 0.8 per assay-disease pair
#   3. Joint: both criteria for the same assay-disease pair
# The figure counts unique genes so repeated assays do not inflate categories.

suppressPackageStartupMessages({
  library(data.table)
})
source("code/09_statistical_analysis/lib/common.R")   # assay_platform, assay_fluid, CHR_LEN_HG19
source("code/09_statistical_analysis/lib/eqtl_coloc_scope.R")   # load_assay_annot

out_summary <- "output/05a_disease_category_summary.csv"
out_pairs <- "output/05a_disease_joint_pairs.csv"
# Full per-row pair table (856k rows) that 05b, 05d, 05e, 06 and 09d read. In output/
# but gitignored on its own line: too large for git, and rerunning 05a rebuilds it.
out_pair_evidence <- "output/05a_disease_pair_evidence.csv.gz"

mr_file <- path.expand(
  file.path(MR_ROOT, "protein2disease/result/mr_protein2disease.tsv.gz")
)
coloc_file <- path.expand(
  file.path(COLOC_ROOT, "protein_x_disease/result/coloc_protein_x_disease.tsv")
)
trait_file <- "info/trait_info.csv"

MR_FDR <- 0.05

trait <- fread(
  trait_file,
  select = c("trait_id", "trait_name", "source", "category")
)
trait <- unique(
  trait[source == "FinnGen meta-analysis",
        .(outcome_id = trait_id, disease = trait_name, disease_category = category)],
  by = "outcome_id"
)

# BF2 Olink / NuLISA assays of the deduped analyzed set, one per protein and fluid
# (SomaLogic dropped; lib/eqtl_coloc_scope.R).
annot <- load_assay_annot()

mr <- fread(
  mr_file,
  select = c(
    "exposure_id", "outcome_id", "status", "n_ivs",
    "ivw_b", "ivw_se", "ivw_p", "steiger_p", "steiger_correct"
  )
)
mr <- merge(mr, annot, by = "exposure_id", all = FALSE)
mr <- merge(mr, trait, by = "outcome_id", all = FALSE)
mr[, mr_ok := status == "ok" & is.finite(ivw_p)]
mr[, mr_fdr := NA_real_]
mr[mr_ok == TRUE, mr_fdr := p.adjust(ivw_p, method = "BH")]
mr[, mr_sig := mr_ok & mr_fdr < MR_FDR]

coloc <- fread(
  coloc_file,
  select = c(
    "exposure_id", "outcome_id", "status", "PP4",
    "coloc_method", "lead_snp", "o_lead", "chr", "pos", "primary"
  )
)
coloc <- merge(coloc, annot, by = "exposure_id", all = FALSE)
coloc <- merge(coloc, trait, by = "outcome_id", all = FALSE)
# The coloc result is sparse: status == "ok" denotes a retained PP4 >= 0.7
# result, not merely an attempted test. no_sig_coloc also records a completed
# attempt, so include both statuses in the tested denominator. Only reported rows
# (coloc_reported, lib/common.R) can colocalize or supply the pair's PP4.
coloc[, coloc_tested_row := status %chin% c("ok", "no_sig_coloc")]
coloc[, coloc_primary_row := coloc_reported(status, primary, PP4)]
coloc[, coloc_strong_row := coloc_primary_row & PP4 > COLOC_PP4]
coloc_eval <- coloc[, .(
  coloc_tested = any(coloc_tested_row),
  coloc_sig = any(coloc_strong_row)
), by = .(exposure_id, outcome_id)]
coloc_best <- coloc[coloc_primary_row == TRUE][
  order(exposure_id, outcome_id, -PP4)
][, .SD[1], by = .(exposure_id, outcome_id)]
coloc_pair <- merge(
  coloc_eval,
  coloc_best[, .(
    exposure_id, outcome_id, PP4, coloc_method, lead_snp, o_lead, chr, pos
  )],
  by = c("exposure_id", "outcome_id"),
  all.x = TRUE
)

pair <- merge(
  mr,
  coloc_pair[, .(
    exposure_id, outcome_id, coloc_tested, coloc_sig, PP4,
    coloc_method, lead_snp, o_lead, chr, pos
  )],
  by = c("exposure_id", "outcome_id"),
  all.x = TRUE
)
pair[is.na(coloc_tested), `:=`(coloc_tested = FALSE, coloc_sig = FALSE)]
pair[, joint_sig := mr_sig & coloc_sig]

# Two complementary summaries per disease category, distinguished by `analysis`:
#
#   "All testable"   every protein for which the test could be run in that fluid. This
#                    describes each fluid's discovery yield. It is NOT a fair
#                    CSF-vs-plasma comparison: the two fluids instrument different
#                    protein sets, so the denominators hold different proteins.
#   "Head-to-head"   restricted to proteins testable in BOTH fluids within that
#                    category, so the denominator is the same protein set on both sides
#                    and the disease GWAS is identical. A difference that survives here
#                    is not driven by protein coverage or instrument availability.
#
# Unlike the eQTL colocalization in 06, the denominator here is not saturated: MR and
# coloc are attempted for every cis-instrumented protein against every disease, so
# "tested" is a broad set (e.g. neurological: 7 of 158 CSF proteins reach joint
# evidence), and restricting it to the shared proteins keeps it broad.
summarize_evidence <- function(d, evidence, tested_col, sig_col, paired = FALSE) {
  if (paired) {
    shared_keys <- d[get(tested_col) == TRUE, .(nf = uniqueN(fluid)),
                     by = .(disease_category, gene)][nf == 2L,
                                                     .(disease_category, gene)]
    d <- merge(d, shared_keys, by = c("disease_category", "gene"))
  }
  tested <- d[get(tested_col) == TRUE, .(
    n_tested_pairs = uniqueN(paste(exposure_id, outcome_id)),
    n_tested_genes = uniqueN(gene)
  ), by = .(fluid, disease_category)]
  sig <- d[get(sig_col) == TRUE, .(
    n_pairs = uniqueN(paste(exposure_id, outcome_id)),
    n_unique_proteins = uniqueN(gene),
    n_diseases = uniqueN(outcome_id)
  ), by = .(fluid, disease_category)]
  out <- merge(tested, sig, by = c("fluid", "disease_category"), all = TRUE)
  for (col in c(
    "n_tested_pairs", "n_tested_genes", "n_pairs",
    "n_unique_proteins", "n_diseases"
  )) {
    set(out, which(is.na(out[[col]])), col, 0L)
  }
  out[, `:=`(
    evidence = evidence,
    analysis = if (paired) "Head-to-head" else "All testable",
    pair_percent = fifelse(
      n_tested_pairs > 0, 100 * n_pairs / n_tested_pairs, NA_real_
    )
  )]
  out
}

pair[, both_tested := mr_ok & coloc_tested]

# Pair-level export: one row per assay x disease, every evidence field this script
# derives. 06_supplementary_tables.R turns this into the protein-disease
# supplementary table; exporting from here rather than re-deriving it there keeps
# one definition of MR_FDR, COLOC_PP4 and the joint call.
dir.create(dirname(out_pair_evidence), recursive = TRUE, showWarnings = FALSE)
fwrite(pair, out_pair_evidence)
message("Saved: ", out_pair_evidence, " (", nrow(pair), " assay x disease rows)")
summary_dt <- rbindlist(list(
  summarize_evidence(pair, "MR", "mr_ok", "mr_sig"),
  summarize_evidence(pair, "Coloc", "coloc_tested", "coloc_sig"),
  summarize_evidence(pair, "MR + coloc", "both_tested", "joint_sig"),
  summarize_evidence(pair, "MR", "mr_ok", "mr_sig", paired = TRUE),
  summarize_evidence(pair, "Coloc", "coloc_tested", "coloc_sig", paired = TRUE),
  summarize_evidence(pair, "MR + coloc", "both_tested", "joint_sig", paired = TRUE)
), use.names = TRUE)

setorder(summary_dt, analysis, fluid, disease_category, evidence)
dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)
fwrite(summary_dt, out_summary)

# One row per fluid-gene-disease, preferring the strongest MR evidence and then
# the highest colocalization posterior across repeated assays.
joint <- pair[joint_sig == TRUE]
setorder(joint, fluid, gene, outcome_id, mr_fdr, -PP4)
joint <- joint[, .SD[1], by = .(fluid, gene, outcome_id)]
joint_out <- joint[, .(
  fluid, gene, outcome_id, disease, disease_category,
  exposure_id, protein_assay, platform,
  n_ivs, mr_beta = ivw_b, mr_se = ivw_se, mr_p = ivw_p, mr_fdr,
  coloc_pp4 = PP4, coloc_method, protein_lead = lead_snp,
  disease_lead = o_lead, chr, pos
)]
setorder(joint_out, fluid, disease_category, mr_fdr, -coloc_pp4)
dir.create(dirname(out_pairs), recursive = TRUE, showWarnings = FALSE)
fwrite(joint_out, out_pairs)

# ---- per-gene evidence, for the fluid-specificity panel ----------------------
# One row per fluid x gene: was the gene testable at all (both methods attempted on
# at least one disease), and did it reach joint evidence. A gene-level table is what
# the specificity panel needs -- the category summary above is already aggregated,
# so the denominator for "proteins of this pQTL class that were tested" cannot be
# recovered from it.

# ---- CSF associations a plasma-only design could not have made ----------------
# Published pQTL-MR is almost entirely plasma, so a joint hit that exists only on the
# CSF side is, operationally, out of reach of that literature. The three reasons are
# kept apart because they are not the same claim: no plasma assay is a measurement
# limit; no plasma cis instrument means plasma has no genome-wide significant cis
# signal, which is not the same as no plasma effect (in the paired sample 9 of the 25
# such proteins have a same-direction plasma effect at P < 0.05, CD58 among them;
# 2026-10-27); a plasma test that simply misses is the weakest and may be power alone.
plasma_genes <- annot[fluid == "Plasma", unique(gene)]
gd <- unique(pair[joint_sig == TRUE, .(gene, outcome_id, disease, disease_category)])
arm <- merge(gd[, .(fluid = c("CSF", "Plasma")), by = names(gd)],
             pair[, .(joint = any(joint_sig), tested = any(mr_ok & coloc_tested)),
                  by = .(gene, outcome_id, fluid)],
             by = c("gene", "outcome_id", "fluid"), all.x = TRUE)
arm[is.na(joint), joint := FALSE][is.na(tested), tested := FALSE]
w <- dcast(arm, gene + outcome_id + disease + disease_category ~ fluid,
           value.var = c("joint", "tested"))
csf_only <- w[joint_CSF == TRUE & joint_Plasma == FALSE]
csf_only[, why_plasma := fcase(
  !(gene %in% plasma_genes), "not on plasma panel",
  tested_Plasma == FALSE, "no plasma cis instrument",
  default = "plasma tested, not significant")]
csf_only <- merge(
  csf_only[, .(gene, outcome_id, disease, disease_category, why_plasma)],
  joint_out[fluid == "CSF", .(gene, outcome_id, n_ivs, mr_beta, mr_fdr, coloc_pp4,
                              protein_lead)],
  by = c("gene", "outcome_id"))
setorder(csf_only, why_plasma, mr_fdr)
out_csf <- file.path(dirname(out_summary), "05a_disease_csf_only_pairs.csv")
fwrite(csf_only, out_csf)

# ---- Joint pairs by fluid-specificity class -------------------------------------
# How a supported gene-disease link is distributed across fluids, on two sets. The
# first is every joint-supported link; the second keeps only links whose protein was
# testable in BOTH fluids for that disease, so "CSF-only" there means the plasma test
# was run and missed, not that it was never possible. `pair` already carries the
# testability flags, so nothing is re-read from the secure tables.
testable_both <- pair[both_tested == TRUE, .(nf = uniqueN(fluid)),
                      by = .(gene, outcome_id)][nf == 2L]
jw <- dcast(joint_out, gene + outcome_id + disease + disease_category ~ fluid,
            value.var = "coloc_pp4", fun.aggregate = length, fill = 0L)
jw[, cls := fcase(CSF > 0 & Plasma > 0, "both",
                  CSF > 0, "CSF-only", default = "plasma-only")]
jw[, testable_both := paste(gene, outcome_id) %chin%
     testable_both[, paste(gene, outcome_id)]]

by_cls <- rbind(
  jw[, .(set = "all joint-supported pairs", n = .N), by = cls],
  jw[testable_both == TRUE, .(set = "pair testable in both fluids", n = .N), by = cls]
)
by_cls[, pct := round(100 * n / sum(n), 1), by = set]
by_cat <- jw[testable_both == TRUE, .(n = .N), by = .(disease_category, cls)]
by_cat[, pct := round(100 * n / sum(n), 1), by = disease_category]

out_xtab <- file.path(dirname(out_summary), "05a_joint_pairs_by_specificity.csv")
fwrite(rbind(by_cls[, .(table = "by_specificity", cls, set, disease_category = NA_character_, n, pct)],
             by_cat[order(disease_category, cls)][
               , .(table = "by_category", cls, set = NA_character_, disease_category, n, pct)]),
       out_xtab)
message("Joint pairs by fluid specificity:")
print(by_cls)
message("Saved: ", out_xtab)

message("MR tests: ", mr[mr_ok == TRUE, .N],
        "; significant at global BH-FDR < 0.05: ", mr[mr_sig == TRUE, .N])
message("Attempted coloc pairs: ", nrow(coloc_pair),
        "; PP4 > ", COLOC_PP4, ": ", coloc_pair[coloc_sig == TRUE, .N])
message("Joint fluid-gene-disease pairs: ", nrow(joint_out))
message("Saved: ", out_summary)
message("Saved: ", out_pairs)
message("Saved: ", out_csf, " (CSF-only pairs: ", nrow(csf_only), ", genes: ",
        uniqueN(csf_only$gene), ")")
print(csf_only[, .(pairs = .N, genes = uniqueN(gene)), by = why_plasma])
