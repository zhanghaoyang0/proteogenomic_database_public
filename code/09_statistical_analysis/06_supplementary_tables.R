#!/usr/bin/env Rscript
# Supplementary data tables: the underlying per-row results behind every figure,
# in three files that share one protein key (gene / ensembl_gene_id / fluid).
#
#   S1  pQTL fine-mapping    one row per credible-set SNP
#   S2  protein x disease    one row per protein x disease x fluid
#   S3  protein x tissue     one row per protein x tissue x fluid x pQTL credible set
#
# These are assembled, not recomputed. S2 comes from the pair table 05a writes
# (output/05a_disease_pair_evidence.csv.gz) and S3 from the in-scope coloc table 04 writes
# (04_protein_tissue_pair_evidence.csv.gz), so the MR FDR threshold, the PP4
# threshold, the joint call and the cis+cognate scope have exactly one definition
# each and cannot drift between the figures and the supplement. S1 is read directly
# from the per-assay fine-mapping files, which already hold only credible-set
# members, and annotated with the cross-fluid sharing 03 computes.
#
# Every table is restricted to the deduped analyzed set (00_kept_assays.csv):
# SomaLogic excluded, GWAS N > 1000, one highest-N assay per protein per fluid. That
# is the set the manuscript counts, so a protein appears at most once per fluid.
#
# Usage: Rscript code/09_statistical_analysis/06_supplementary_tables.R

suppressPackageStartupMessages({
  library(data.table)
})

OUT_DIR  <- "output"
KEPT     <- file.path(OUT_DIR, "00_kept_assays.csv")
SPEC     <- file.path(OUT_DIR, "03_fluid_specificity_protein.csv")
CS_SHARE <- file.path(OUT_DIR, "03_fluid_specificity_credible_sets.csv")
PAIR     <- "output/05a_disease_pair_evidence.csv.gz"
TISSUE   <- file.path(OUT_DIR, "04_protein_tissue_pair_evidence.csv.gz")
CS_LEAD  <- "info/bf_protein/finemap_lead_snps.csv"
REPL_TSV <- "info/bf_protein/replication/cs_replication.tsv.gz"
source("code/09_statistical_analysis/lib/common.R")   # GWAS_BASE, read_finemap_cs, check_finemap_cs, add_novelty, add_western
FM_DIR   <- file.path(GWAS_BASE, "bf2/finemapping/protein")

out_s1 <- file.path(OUT_DIR, "S1_pqtl_finemapping.csv.gz")
out_s2 <- file.path(OUT_DIR, "S2_protein_disease.csv.gz")
# S2 is the evidence table: the rows clearing MR FDR < 0.05 or PP4 > 0.8. The full
# table is every MR test that ran -- 856k rows, ~33 MB gzipped, of which 850k are
# null -- and adds only the negative record, so it is written on request rather than
# every run. Instrument availability, the other thing the full table carried, is a
# property of protein x fluid rather than of protein x disease and is now a set of
# columns on S1, where it appears once per protein instead of 554 times.
#   Rscript code/09_statistical_analysis/06_supplementary_tables.R --full-s2
WRITE_FULL_S2 <- "--full-s2" %chin% commandArgs(trailingOnly = TRUE)
out_s2_full <- file.path(OUT_DIR, "S2_protein_disease_all_tested.csv.gz")
out_s3 <- file.path(OUT_DIR, "S3_protein_tissue.csv.gz")

for (f in c(KEPT, SPEC, CS_SHARE, PAIR, TISSUE, CS_LEAD, REPL_TSV)) {
  if (!file.exists(f)) {
    stop("missing input: ", f,
         "\n  run 00, 03, 04 and 05a first (see code/README.md)")
  }
}
if (!dir.exists(FM_DIR)) stop("missing fine-mapping directory: ", FM_DIR)

source("code/09_statistical_analysis/lib/cs_signals.R")

kept <- fread(KEPT)
# gwas_n is the GWAS median per-variant N (00b), not 00's N, which counts every
# measured sample, genotyped or not.
gwas_qc <- fread(file.path(OUT_DIR, "00b_gwas_qc.csv"), select = c("protein_assay", "gwas_n"))
kept <- merge(kept, gwas_qc, by = "protein_assay")
if (anyNA(kept$gwas_n) || nrow(kept) != nrow(fread(KEPT))) stop("00b_gwas_qc.csv does not cover 00_kept_assays.csv; run 00b")
kept_key <- unique(kept[, .(cohort, protein_assay, fluid, platform,
                            ensembl_gene_id, gene, gwas_n)])

# Protein-level fluid class from 09: shared / CSF-only / plasma-only / neither,
# defined over the 1,202 proteins assayed in both fluids. Proteins assayed in one
# fluid only are absent from that file and are labelled here, because "not assayed
# in the other fluid" and "assayed but no pQTL there" are different statements.
spec <- fread(SPEC)[, .(ensembl_gene_id, protein_fluid_class = class)]
cis_both_genes <- fread(SPEC)[cis_both == TRUE, ensembl_gene_id]

add_class <- function(d) {
  d <- merge(d, spec, by = "ensembl_gene_id", all.x = TRUE)
  d[is.na(protein_fluid_class), protein_fluid_class := "assayed in one fluid only"]
  d
}

# Instrument availability, per protein x fluid. It is a property of the protein and
# its fluid, not of the protein-disease pair, so it belongs on S1 once rather than on
# S2 once per disease. Read here because S1 is written first.
if (!file.exists(PAIR)) {
  stop("missing input: ", PAIR, "\n  run 05a_disease_data.R first")
}
pair_dt <- fread(PAIR)
avail <- pair_dt[, .(
  mr_instrumented        = any(mr_ok),
  n_diseases_mr_tested   = sum(mr_ok),
  n_diseases_coloc_tested = sum(coloc_tested),
  n_diseases_joint       = sum(joint_sig)
), by = protein_assay]

# ---------------------------------------------------------------------- S1 ----
# The per-assay .finemapped.tsv.gz files hold only credible-set members, so reading
# them whole gives the SNP-level table directly; there is nothing to filter out.
s1 <- rbindlist(Map(read_finemap_cs, kept_key$cohort, kept_key$protein_assay,
                    MoreArgs = list(fm_dir = FM_DIR)),
                use.names = TRUE, fill = TRUE)
if (nrow(s1) == 0) stop("no credible sets read from ", FM_DIR)
check_finemap_cs(s1, kept_key, FM_DIR)
# converged: per SuSiE fit (max_iter 1000). Non-converged fits keep their credible sets in
# S1, flagged; coloc does not use them (functions_coloc.R add_component_flags).
if (!"converged" %in% names(s1) || anyNA(s1$converged))
  stop("fine-mapping files without a converged flag for every row")
s1 <- merge(s1, kept_key, by = c("cohort", "protein_assay"))

# cs_key is 03's key (fluid | ensembl_gene_id | cs_id), so the sharing flag joins
# without recomputing which credible sets overlap between fluids.
s1[, cs_key := paste(fluid, ensembl_gene_id, cs_id, sep = "|")]
cs_share <- fread(CS_SHARE)[, .(cs_key, cs_shared_with_other_fluid = shared)]
s1 <- merge(s1, cs_share, by = "cs_key", all.x = TRUE)

# Credible-set summary fields (size, purity, max PIP, whether the set clears
# genome-wide significance) live one level up, in the sentinel table.
# is_cis comes with the fine-mapping files: 01_finemapping.R sets it per credible set
# (top-PIP variant in the cis window), and the sentinel table carries the same value.
# (Until 2026-10-06 the files held the clump locus's status instead, and S1 exported
# both flags.)
cs_lead <- fread(CS_LEAD)[cohort == "bf2", .(
  cohort, protein_assay, cs_id,
  cs_lead_snp = snp, cs_size, cs_max_pip = max_pip, cs_hc_p5e8 = hc_p5e8
)]
cs_lead <- unique(cs_lead, by = c("cohort", "protein_assay", "cs_id"))
s1 <- merge(s1, cs_lead, by = c("cohort", "protein_assay", "cs_id"), all.x = TRUE)
# S1 releases every credible set, including the ones overlapping +/-500 kb windows
# resolved more than once, but flags which set represents the distinct signal so a
# reader can reproduce the counts in the main text.
s1 <- add_cs_signal(s1)
# Reported = counted in the paper: the representative set of a distinct signal, and for
# trans its lead reaches the study-wide threshold (lib/cs_signals.R).
s1[, signal_reported := is_representative & is_reported]
s1[, is_cs_lead := snp == cs_lead_snp]

# BioFINDER-1 replication status, per credible set. Joined on the credible-set lead
# variant rather than on cs_id: the replication step submits one row per protein x
# sentinel variant, so two credible sets of the same protein that share a lead are one
# submitted test, and joining on the lead attaches that one result to both instead of
# leaving the second looking untested. What is left unjoined afterwards is exactly the
# credible sets whose lead (top-PIP variant) is not genome-wide significant (hc_p5e8),
# which were never eligible; asserted below. The four outcomes are kept distinct because "did not replicate" and "could
# not be tested" are different statements and pooling them would flatter the rate.
#
# bf1_status is read, not recomputed: code/05_replication/04_credible_set_status.R
# owns the rule and asserts that the funnel it produces is still the published
# Table 1. It used to be rebuilt here, in 01a and in 01b, which is three chances
# for the three to drift apart.
repl <- unique(fread(REPL_TSV, select = c(
  "trait", "cs_lead_snp", "bf1_status", "bf1_assay", "bf1_assay_match",
  "bf1_beta", "bf1_se", "bf1_p"
)), by = c("trait", "cs_lead_snp"))
s1[, trait := paste(cohort, protein_assay, sep = "__")]
s1 <- merge(s1, repl, by = c("trait", "cs_lead_snp"), all.x = TRUE)
s1[, trait := NULL]
if (s1[is.na(bf1_status) & cs_hc_p5e8 %in% TRUE, .N])
  stop("credible sets with a genome-wide significant lead have no BF1 status")
s1[is.na(bf1_status), bf1_status := "not submitted: lead variant not genome-wide significant"]
# No BF1 assay is NA in S1 whatever the reason: the replication table writes "" for a set
# it tested without one, while a set it never received gets NA from the merge (r024).
s1[bf1_assay == "", bf1_assay := NA_character_]
s1 <- merge(s1, avail, by = "protein_assay", all.x = TRUE)
# A protein absent from the MR table was never instrumented against disease at all.
s1[is.na(mr_instrumented), `:=`(mr_instrumented = FALSE, n_diseases_mr_tested = 0L,
                                n_diseases_coloc_tested = 0L, n_diseases_joint = 0L)]
s1 <- add_class(s1)

s1_out <- s1[, .(
  cohort, fluid, platform, protein_assay, gene, ensembl_gene_id, gwas_n,
  cs_id, cs_key, is_cis,
  cs_size, cs_purity, cs_max_pip, cs_hc_p5e8, signal_reported,
  clump_lead_snp, cs_lead_snp, is_cs_lead,
  snp, chr, pos, ea, nea, beta, se, p, pip, susie_converged = converged,
  cs_shared_with_other_fluid, protein_fluid_class,
  bf1_status, bf1_assay, bf1_assay_match, bf1_beta, bf1_se, bf1_p,
  mr_instrumented, n_diseases_mr_tested, n_diseases_coloc_tested, n_diseases_joint
)]
setorder(s1_out, fluid, gene, chr, pos)
# Subnormal P values (below ~2.2e-308) lose digits on every fread/fwrite round trip, and
# 07b re-reads and rewrites S1, so an unchanged rerun showed S1 as modified. Floor them at
# 1e-307: double.xmin itself prints with 15 digits as a value that parses back subnormal.
s1_out[p > 0 & p < 1e-307, p := 1e-307]
s1_nov <- add_novelty(s1_out, file.path(OUT_DIR, "07_novelty_credible_sets.csv.gz"))
if (is.null(s1_nov)) {
  warning("S1 written without the novelty columns (credible sets changed since 07b ran); ",
          "run 07b_novelty_scan.R next", call. = FALSE)
} else s1_out <- s1_nov
s1_w <- add_western(s1_out, file.path(OUT_DIR, "07c_western_replication.csv"))
if (is.null(s1_w)) {
  warning("S1 written without the Western replication columns (credible sets changed since 07c ran); ",
          "run 07c_western_replication.R next", call. = FALSE)
} else s1_out <- s1_w
fwrite(s1_out, out_s1)
message("Saved: ", out_s1, " (", nrow(s1_out), " credible-set SNPs, ",
        uniqueN(s1_out$cs_key), " credible sets, ",
        uniqueN(s1_out[, paste(fluid, ensembl_gene_id)]), " protein x fluid)")
message("  credible sets by fluid and cis/trans (is_cis, the manuscript definition):")
print(unique(s1_out[, .(cs_key, fluid, is_cis)])[, .N, by = .(fluid, is_cis)][
  order(fluid, -is_cis)])
message("  BioFINDER-1 replication status, per credible set:")
print(unique(s1_out[, .(cs_key, fluid, bf1_status)])[, .N, by = .(fluid, bf1_status)][
  order(fluid, -N)])

# ---------------------------------------------------------------------- S2 ----
s2 <- copy(pair_dt)
s2 <- merge(s2, kept_key[, .(protein_assay, ensembl_gene_id, gwas_n)],
            by = "protein_assay", all.x = TRUE)

# Steiger directionality: the MR table carries it for every pair (via 05a); S2 reports
# it for the MR-significant pairs only, the set Methods describe testing.
s2[!(mr_sig %in% TRUE), c("steiger_p", "steiger_correct") := NA]
s2[, steiger_pass := steiger_correct & steiger_p < 0.05]
if (s2[mr_sig %in% TRUE & is.na(steiger_pass), .N] > 0L)
  stop("MR-significant pairs without Steiger; the MR table predates 03's Steiger columns (2026-10-25)")

# Link class: for a gene x disease that reaches joint evidence, whether both fluids
# reach it or only one. Rows without joint evidence get NA rather than a class,
# because "not supported here" is not the same as "supported only elsewhere".
jl <- s2[joint_sig == TRUE, .(nf = uniqueN(fluid), fl = paste(sort(unique(fluid)),
                                                              collapse = "+")),
         by = .(gene, outcome_id)]
jl[, link_class := fcase(nf == 2L, "both fluids",
                         fl == "CSF", "CSF-only",
                         default = "plasma-only")]
s2 <- merge(s2, jl[, .(gene, outcome_id, link_class)],
            by = c("gene", "outcome_id"), all.x = TRUE)

# The head-to-head restriction, recorded per row rather than only in the summary, so
# a reader can reproduce the restricted panel of 05_supp_disease.png a from this file alone. Same
# definition as summarize_evidence(paired = TRUE) in 03a: tested in both fluids
# within the disease category.
shared_keys <- s2[both_tested == TRUE, .(nf = uniqueN(fluid)),
                  by = .(disease_category, gene)][nf == 2L,
                                                  .(disease_category, gene,
                                                    in_head_to_head = TRUE)]
s2 <- merge(s2, shared_keys, by = c("disease_category", "gene"), all.x = TRUE)
s2[is.na(in_head_to_head), in_head_to_head := FALSE]
s2 <- add_class(s2)

s2_out <- s2[, .(
  fluid, platform, protein_assay, gene, ensembl_gene_id, gwas_n, exposure_id,
  outcome_id, disease, disease_category,
  mr_status = status, mr_testable = mr_ok, n_ivs,
  mr_beta = ivw_b, mr_se = ivw_se, mr_p = ivw_p, mr_fdr, mr_sig,
  steiger_p, steiger_correct, steiger_pass,
  coloc_tested, coloc_pp4 = PP4, coloc_method, coloc_sig,
  protein_lead_snp = lead_snp, disease_lead_snp = o_lead, chr, pos,
  tested_by_both_methods = both_tested, joint_evidence = joint_sig,
  link_class, in_head_to_head, protein_fluid_class
)]
setorder(s2_out, -joint_evidence, disease_category, gene, fluid, outcome_id)
s2_sup <- s2_out[mr_sig == TRUE | coloc_sig == TRUE]
fwrite(s2_sup, out_s2)
message("Saved: ", out_s2, " (", nrow(s2_sup),
        " rows clearing MR FDR < 0.05 or PP4 > ", COLOC_PP4, ": ",
        s2_sup[mr_sig & coloc_sig, .N], " both, ",
        s2_sup[mr_sig & !coloc_sig, .N], " MR only, ",
        s2_sup[!mr_sig & coloc_sig, .N], " coloc only)")
if (WRITE_FULL_S2) {
  fwrite(s2_out, out_s2_full)
  message("Saved: ", out_s2_full, " (", nrow(s2_out), " tested rows)")
}

# ---------------------------------------------------------------------- S3 ----
s3 <- fread(TISSUE)
s3 <- merge(s3, kept_key[, .(protein_assay, ensembl_gene_id, gwas_n)],
            by = "protein_assay", all.x = TRUE)
# 06 keeps status-only rows because they record an attempted test; the same two
# statuses define the tested denominator here.
s3[, coloc_tested := status %chin% c("ok", "no_sig_coloc")]
s3[, coloc_sig := coloc_retained == TRUE & PP4 > COLOC_PP4]
# The proteins carrying a cis-pQTL in both fluids (03, cis_both): the fixed denominator
# of the head-to-head counts (04_tissue_dumbbell.csv, drawn by 10d). Flagged per row so that panel is reproducible here.
s3 <- add_class(s3)
s3[, in_head_to_head := ensembl_gene_id %in% cis_both_genes]

s3_out <- s3[, .(
  fluid, platform, protein_assay, gene, ensembl_gene_id, gwas_n, exposure_id,
  resource, tissue, tissue_class, eqtl_ensg = ENSG, eqtl_gene = GENE,
  pqtl_cs_id = cs_id, is_cis, pqtl_lead_snp = lead_snp, eqtl_lead_snp = o_lead,
  chr, pos, coloc_status = status, coloc_primary = primary, coloc_tested, coloc_pp4 = PP4,
  coloc_method, coloc_sig, in_head_to_head, protein_fluid_class
)]
setorder(s3_out, -coloc_sig, resource, tissue_class, tissue, gene, fluid)
fwrite(s3_out, out_s3)
message("Saved: ", out_s3, " (", nrow(s3_out), " protein x tissue x fluid rows, ",
        s3_out[coloc_sig == TRUE, .N], " colocalizing at PP4 > ", COLOC_PP4, ")")
message("06 done")
