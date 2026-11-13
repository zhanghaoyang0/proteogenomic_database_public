#!/usr/bin/env Rscript
# The six apparently unreported CSF-supported protein-disease pairs, flagged in
# Supplementary Table 4 (no_prior_report; their forest panel was dropped 2026-11-13),
# and the record of why these pairs and not others.
# It computes the evidence; 05a_disease_data.R decides what is supported.
#
# (Until 2026-10-29 this file also wrote the three PIP loci of an earlier panel,
# 05b_case_tracks.csv / 05b_case_studies.csv; the locus panels are now 05f / 05g.)
#
# Usage:
#   Rscript code/09_statistical_analysis/05b_case_data.R [out_csv]

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/common.R")   # assay_platform, assay_fluid
source("code/09_statistical_analysis/lib/eqtl_coloc_scope.R")   # load_assay_annot

args <- commandArgs(trailingOnly = TRUE)

# CSF-vs-plasma MR estimates for the three CSF-only pairs that a literature search
# found no prior protein-level report for (panel c).
#
# What this script can and cannot establish. It computes the evidence, not the
# novelty: whether a pair has been reported before was decided by hand, by literature
# search, and is recorded here only as the reason the pair was picked. The wording the
# figure and the manuscript use is deliberately narrow -- "apparently previously
# unreported genetically supported protein-disease association" -- because "no one has
# ever linked this gene to this disease" and "no one has shown the protein level is
# genetically linked to this disease" are different claims, and only the second is
# ours. SPON2 has glaucoma expression literature back to 2006; TMPRSS5 sits in a
# neuropathy locus reported in 2026. Neither is a new gene-disease link.
#
# Why these six and not the other 133 CSF-only pairs:
#   SPON2   - no plasma cis instrument, so plasma genetics could not have tested it at
#             all. 2 IVs, and the second hit is the wide definition of the same
#             disease rather than an unrelated phenotype.
#   CGREF1  - plasma was tested and did not reach joint evidence. 2 IVs, one disease.
#             The cleanest form of the compartment claim.
#   TMPRSS5 - plasma tested, not significant, at a disease locus already known. Kept
#             on purpose as the honest case: new protein-level layer, not a new gene.
#   ICAM5   - another disease class, and the strongest CSF statistic here. Its plasma
#             arm agrees in direction and is nominally significant, failing only the
#             global BH threshold, so it argues CSF precision rather than a plasma
#             blind spot. Single IV in the ICAM1/ICAM4/ICAM5 cluster, so protein
#             attribution needs care; rs892188 also carries diverticular disease.
#   NRP1    - same class as ICAM5 (plasma agrees nominally), kept for the sleep/
#             neurological angle. Single IV, which also carries seborrhoeic keratosis.
#   ITIH4   - a sixth disease class (musculoskeletal) and, like SPON2, a pair plasma
#             genetics could not have tested: no plasma cis instrument. Replaced
#             FAS -> CLL, which was picked first but whose CSF estimate (beta -0.48)
#             was three times any other and alone set the x axis, flattening the
#             other five rows. Single IV, rs4687657, which also carries hernia, so
#             this is a locus with more than one outcome and attribution to the
#             protein needs the same care as ICAM5 and NRP1.
#
# An audit of all 139 CSF-only pairs is what these labels rest on. Among the 39 where
# plasma was actually tested, 22 have a plasma estimate in the same direction that is
# nominally significant and fails only the global BH threshold; only 13 are genuinely
# null. "Plasma tested, not significant" is therefore usually power, not biology, and
# the panel is titled and captioned so it does not claim otherwise.
#
# The MR FDR is read from 05a_disease_data.R's pair table, not recomputed: the numbers
# printed on the figure are the ones the disease-category counts rest on, over the same
# BH family (FinnGen outcomes; code review 2026-10-30, R9).
#
# Output:
#   05b_novel_pairs.csv
out_novel_csv <- if (length(args) >= 1) args[1] else
  "output/05b_novel_pairs.csv"

pair_file <- "output/05a_disease_pair_evidence.csv.gz"
coloc_file <- path.expand(
  file.path(COLOC_ROOT, "protein_x_disease/result/coloc_protein_x_disease.tsv")
)
csf_only_file <- "output/05a_disease_csf_only_pairs.csv"

cases <- data.table(
  case_order = 1:6,
  gene = c("SPON2", "CGREF1", "TMPRSS5", "ICAM5", "NRP1", "ITIH4"),
  outcome_id = c(
    "H7_GLAUCOMA_POAG_meta_fin_ukbb_mvp",
    "T2D_meta_fin_ukbb_mvp",
    "G6_POLYOTHUNS_meta_fin_ukbb_mvp",
    "I9_HYPTENSESS_meta_fin_ukbb_mvp",
    "G6_SLEEPAPNO_meta_fin_ukbb_mvp",
    "finngen_R12_M13_ARTHROSIS_meta_fin_ukbb"
  ),
  # FinnGen names carry their definition suffix, which is noise on an axis label.
  disease_short = c(
    "Glaucoma",
    "Type 2 diabetes",
    "Polyneuropathy",
    "Hypertension",
    "Sleep apnoea",
    "Arthrosis"
  )
)

# ---- same exposure basis as 05a_disease_data.R -------------------------------
annot <- load_assay_annot()

mr <- fread(pair_file, select = c("exposure_id", "outcome_id", "status", "n_ivs",
                                 "ivw_b", "ivw_se", "ivw_p", "mr_ok", "mr_fdr", "mr_sig"))
mr <- merge(mr, annot, by = "exposure_id", all = FALSE)

coloc <- fread(coloc_file, select = c("exposure_id", "outcome_id", "status", "PP4",
                                      "lead_snp", "primary"))
coloc <- merge(coloc, annot[, .(exposure_id)], by = "exposure_id", all = FALSE)
# A missing PP4 means two things: "no_sig_coloc" was tested and fell below the
# pipeline's PP4_KEEP (code/07_colocalization/functions_coloc.R), anything else was not
# tested. 05c labels them differently, so the status travels with the row. An ok row
# that is not primary (code/07_colocalization/README.md) is tested but not shown.
coloc[, primary_ok := coloc_reported(status, primary, PP4)]
coloc_status <- coloc[, .(coloc_status = if (any(primary_ok)) "ok"
                                         else if ("ok" %in% status) "not_primary"
                                         else status[1]),
                      by = .(exposure_id, outcome_id)]
coloc <- coloc[primary_ok == TRUE]
setorder(coloc, exposure_id, outcome_id, -PP4)
coloc <- coloc[, .SD[1], by = .(exposure_id, outcome_id)]

d <- merge(mr[mr_ok == TRUE], cases, by = c("gene", "outcome_id"), all = FALSE)
n_d <- nrow(d)
d <- merge(d, coloc[, .(exposure_id, outcome_id, pp4 = PP4, lead_snp)],
           by = c("exposure_id", "outcome_id"), all.x = TRUE)
d <- merge(d, coloc_status, by = c("exposure_id", "outcome_id"), all.x = TRUE)
stopifnot(nrow(d) == n_d)
# One assay per gene x fluid survives the dedup, so a second row here means the
# dedup did not do its job and the figure would silently draw one of two estimates.
if (anyDuplicated(d[, .(case_order, fluid)])) {
  stop("More than one assay per case x fluid after dedup")
}
# Every case is a CSF-supported pair; a case without a CSF row (assay or MR result
# gone) would drop out of the figure silently (code review 2026-10-30, R14).
if (!setequal(d[fluid == "CSF", case_order], cases$case_order)) {
  stop("Cases without a CSF MR row: ",
       paste(cases[!case_order %in% d[fluid == "CSF", case_order], gene], collapse = ", "))
}

# ---- why the plasma arm is missing where it is missing -----------------------
# Taken from 02b rather than re-derived: the three reasons are not equivalent and the
# figure must not blur "plasma has no instrument" into "plasma was tested and missed".
csf_only <- fread(csf_only_file, select = c("gene", "outcome_id", "why_plasma"))
d <- merge(d, unique(csf_only), by = c("gene", "outcome_id"), all.x = TRUE)
if (anyNA(d[fluid == "CSF", why_plasma])) {
  stop("A case is not in the CSF-only table: ",
       paste(d[fluid == "CSF" & is.na(why_plasma), gene], collapse = ", "))
}

out <- d[, .(
  case_order, gene, disease_short, outcome_id, fluid, platform, protein_assay,
  why_plasma, n_ivs,
  mr_beta = ivw_b, mr_se = ivw_se, mr_p = ivw_p, mr_fdr,
  mr_sig,
  pp4, coloc_status, lead_snp
)]
setorder(out, case_order, -fluid)  # CSF first within a case

fwrite(out, out_novel_csv)
message("Saved: ", out_novel_csv, " (", nrow(out), " rows)")
print(out[, .(case_order, gene, fluid, n_ivs, mr_beta = round(mr_beta, 3),
              mr_fdr = signif(mr_fdr, 3), pp4 = round(pp4, 3), mr_sig)])
