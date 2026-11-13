#!/usr/bin/env Rscript
# Candidate list for the novel-pQTL validation (docs/pipeline.md):
# one row per distinct CSF credible-set signal, flagged is_candidate when it is outside
# the MHC and previously unmatched (a CSF reference, Cruchaga SomaScan or Western Olink,
# assayed the protein and neither reported anything matching):
#   cis    the protein's primary association, A or B in 07b's novelty_primary_cis, and not
#          genome-wide significant in Western's full summary statistics (07c); secondary
#          cis signals are not candidates
#   trans  tier A in 07b's locus-level call
# assoc_class carries the call: A, B, "B, detected in Western", known, unresolved, C, and
# "secondary" for cis signals that are not their protein's primary association.
#
# Every CSF signal is written, not only the candidates; that served the local
# regressions (09b, archived 2026-10-31). 09d and 09e read only the candidates.
#
# CSF only: plasma tier A is mostly trans signals UKB-PPP (34x our N) did not call,
# which the novelty plan reads as likely false positives, not discoveries.
#
# This file is summary-level (no individual data) and is committed, so the Mac
# picks it up with `git pull`.
#
# Usage: Rscript code/09_statistical_analysis/09a_novel_candidates.R

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/cs_signals.R")

OUT_DIR <- "output"
S1      <- file.path(OUT_DIR, "S1_pqtl_finemapping.csv.gz")
MHC     <- "code/lib/mhc_hg19.range"
OUT     <- file.path(OUT_DIR, "09a_novel_candidates.csv")

for (f in c(S1, MHC)) if (!file.exists(f)) stop("missing input: ", f)

s1 <- fread(S1)
if (!"novelty" %in% names(s1)) stop("S1 has no novelty column -- run 07b_novelty_scan.R")

lead <- s1[is_cs_lead == TRUE & fluid == "CSF"]
if (anyDuplicated(lead$cs_key)) stop("more than one lead per credible set")

# S1 releases every credible set; overlapping fine-mapping windows resolve one
# signal more than once, so keep the representative set only.
n0 <- nrow(lead)
lead <- distinct_signals(lead)
message("CSF credible-set leads: ", n0, " -> ", nrow(lead), " distinct signals")

mhc <- fread(MHC, header = FALSE, col.names = c("chr", "start", "end", "label"))
lead[, in_mhc := chr == mhc$chr[1] & pos >= mhc$start[1] & pos <= mhc$end[1]]

# protein_assay is "<protein>__<panel tag>" and is itself the column name in the
# GWAS .phe; the tag names that file (bf2_<tag>.phe) and the mean_protein_<tag>
# covariate the local scripts need.
lead[, panel_tag := sub("^.*__", "", protein_assay)]
if (lead[panel_tag == protein_assay, .N] > 0)
  stop("protein_assay not of the form <protein>__<panel>")

if (!all(c("novelty_primary_cis", "western_cs_min_p") %in% names(s1)))
  stop("S1 lacks novelty_primary_cis / western_cs_min_p -- run 07b and 07c")
lead[, assoc_class := fcase(
  !is_cis, sub(":.*", "", novelty),
  is.na(novelty_primary_cis) | novelty_primary_cis == "", "secondary",   # fread reads S1 NA as ""
  startsWith(novelty_primary_cis, "B") & (western_cs_min_p < 5e-8) %in% TRUE,
    "B, detected in Western",
  default = sub(":.*", "", novelty_primary_cis))]
lead[, is_candidate := !in_mhc & ((is_cis & assoc_class %in% c("A", "B")) |
                                  (!is_cis & startsWith(novelty, "A:")))]

out <- lead[, .(cs_key, protein_assay, panel_tag, platform, gene,
                ensembl_gene_id, is_cis, snp, chr, pos, ea, nea, beta, se, p,
                bf1_status, novelty, assoc_class, in_mhc, is_candidate)]
setorder(out, -is_candidate, gene, chr, pos)
fwrite(out, OUT)

message("\nsignals by cis/trans x novelty:")
print(out[, .N, by = .(is_cis, assoc_class)][order(is_cis, assoc_class)])
message("\ncandidates: ", out[is_candidate == TRUE, .N], " signals on ",
        out[is_candidate == TRUE, uniqueN(protein_assay)], " proteins (",
        out[is_candidate == TRUE & is_cis == TRUE, .N], " cis, ",
        out[is_candidate == TRUE & is_cis == FALSE, .N], " trans); excluded for MHC: ",
        out[assoc_class %in% c("A", "B") & in_mhc == TRUE, .N])
message("wrote ", OUT)
