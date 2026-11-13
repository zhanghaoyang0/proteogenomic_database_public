#!/usr/bin/env Rscript
# Do the MR instruments and the colocalizing locus mark the same signal, for the
# jointly supported protein -> disease pairs (05a's joint_sig)? Joint support is a
# pair-level call: MR combines every cis instrument, coloc counts the pair once any cis
# locus has PP4 > COLOC_PP4. This script checks the two meet at the same locus.
#
# A coloc row's lead_snp is the protein clump lead that defines its SuSiE locus
# (lead +/- FM_WINDOW_BP); the MR instruments are cis clump leads of the same clumping.
# Per pair, over its reported coloc rows with PP4 > COLOC_PP4:
#   instrument     a colocalizing locus's lead is an MR instrument
#   in_window      not, but an instrument lies inside a colocalizing locus's window
#   distinct       neither: MR and coloc rest on different cis signals
# Instruments are mr_iv.tsv.gz harmonized against the outcome with the MR runner's own
# harmonize() (code/08_mendelian_randomization/functions_mr.R); a pair whose count differs from the MR
# table's n_ivs stops the script.
#
# Input:  output/05a_disease_pair_evidence.csv.gz (run 05a first), the MR instrument table,
#         the FinnGen clean sumstats and the coloc result under ~/db/shared_data,
#         info/bf_protein/clump_lead_pqtl.csv.gz (lead positions)
# Output: output/05e_joint_iv_match.csv, one row per joint exposure x outcome
# Usage:  Rscript code/09_statistical_analysis/05e_joint_iv_match.R  (repo root;
#         one outcome file per core, cores from SLURM_CPUS_PER_TASK)

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})
source("code/08_mendelian_randomization/functions_mr.R")                       # harmonize
source("code/09_statistical_analysis/lib/common.R")     # COLOC_PP4, coloc_reported, FM_WINDOW_BP

PAIR     <- "output/05a_disease_pair_evidence.csv.gz"
OUT      <- "output/05e_joint_iv_match.csv"
IV_FILE  <- file.path(MR_ROOT, "protein2disease/info/mr_iv.tsv.gz")
OC_DIR   <- file.path(GWAS_BASE, "finngen/clean")
COLOC    <- file.path(COLOC_ROOT, "protein_x_disease/result/coloc_protein_x_disease.tsv")
LEAD     <- "info/bf_protein/clump_lead_pqtl.csv.gz"
CORES    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

for (f in c(PAIR, IV_FILE, COLOC, LEAD)) if (!file.exists(f)) stop("Missing input: ", f)

pair <- fread(PAIR)
need <- c("exposure_id", "outcome_id", "protein_assay", "gene", "fluid", "disease", "n_ivs", "joint_sig")
if (!all(need %in% names(pair))) stop(PAIR, " lacks ", paste(setdiff(need, names(pair)), collapse = ", "))
pair <- unique(pair[joint_sig %in% TRUE, ..need][, joint_sig := NULL])
if (anyDuplicated(pair, by = c("exposure_id", "outcome_id"))) stop("a joint pair appears twice")
message("joint pairs: ", nrow(pair), " (", uniqueN(pair$outcome_id), " outcomes)")

# Colocalizing loci: the rows 05a counts as strong (reported, PP4 > COLOC_PP4).
coloc <- fread(COLOC, select = c("exposure_id", "outcome_id", "status", "primary", "PP4", "lead_snp"))
coloc <- unique(coloc[coloc_reported(status, primary, PP4) & PP4 > COLOC_PP4,
                      .(exposure_id, outcome_id, lead_snp)])
coloc <- coloc[pair[, .(exposure_id, outcome_id)], on = c("exposure_id", "outcome_id"), nomatch = NULL]
if (uniqueN(coloc[, .(exposure_id, outcome_id)]) != nrow(pair))
  stop("joint pairs without a strong coloc row: ", nrow(pair) - uniqueN(coloc[, .(exposure_id, outcome_id)]))

lead <- unique(fread(LEAD, select = c("cohort", "protein_assay", "SNP", "CHR", "BP"))[cohort == "bf2"],
               by = c("protein_assay", "SNP"))

iv_all <- fread(IV_FILE, colClasses = c(SNP = "character", EA_e = "character", NEA_e = "character"),
                na.strings = c("", "NA"))   # fwrite writes the no-IV row's NA SNP as ""
iv_all <- split(iv_all[!is.na(SNP) & exposure_id %in% pair$exposure_id], by = "exposure_id", keep.by = FALSE)

by_outcome <- split(pair, by = "outcome_id")
res <- mclapply(names(by_outcome), function(oid) {
  p <- by_outcome[[oid]]
  f <- file.path(OC_DIR, paste0(oid, ".tsv.gz"))
  if (!file.exists(f)) stop("outcome file not found: ", f)
  snps <- unique(unlist(lapply(p$exposure_id, function(e) iv_all[[e]]$SNP)))
  oc <- fread(f, select = c("SNP", "EA", "NEA", "EA_FREQ", "BETA", "SE", "P", "N"))[SNP %in% snps]
  rbindlist(lapply(seq_len(nrow(p)), function(i) {
    hm <- harmonize(iv_all[[p$exposure_id[i]]], oc)
    if (is.null(hm) || nrow(hm) != p$n_ivs[i])
      stop(p$exposure_id[i], " x ", oid, ": ", if (is.null(hm)) 0L else nrow(hm),
           " harmonized IVs, MR used ", p$n_ivs[i])
    ivs <- lead[data.table(protein_assay = p$protein_assay[i], SNP = hm$SNP), on = c("protein_assay", "SNP")]
    cl  <- lead[data.table(protein_assay = p$protein_assay[i],
                           SNP = coloc[data.table(exposure_id = p$exposure_id[i], outcome_id = oid),
                                       lead_snp, on = c("exposure_id", "outcome_id")]),
                on = c("protein_assay", "SNP")]
    if (anyNA(ivs$BP) || anyNA(cl$BP)) stop(p$exposure_id[i], " x ", oid, ": a lead has no position in ", LEAD)
    d <- CJ(iv = seq_len(nrow(ivs)), lc = seq_len(nrow(cl)))
    d[, kb := fifelse(ivs$CHR[iv] == cl$CHR[lc], abs(ivs$BP[iv] - cl$BP[lc]) / 1e3, Inf)]
    top_iv <- hm$SNP[which.max((hm$BETA_e / hm$SE_o)^2)]   # largest IVW weight
    data.table(
      exposure_id = p$exposure_id[i], outcome_id = oid,
      iv_match = fcase(any(cl$SNP %chin% hm$SNP), "instrument",
                       min(d$kb) <= FM_WINDOW_BP / 1e3, "in_window",
                       default = "distinct"),
      top_iv_is_coloc_lead = top_iv %chin% cl$SNP,
      min_iv_lead_kb = min(d$kb),
      n_coloc_loci = nrow(cl),
      ivs = paste(hm$SNP, collapse = "|"), coloc_leads = paste(cl$SNP, collapse = "|"))
  }))
}, mc.cores = CORES)
bad <- which(!vapply(res, is.data.table, logical(1)))
if (length(bad)) stop(length(bad), " outcome(s) failed, e.g. ", names(by_outcome)[bad[1]], ": ",
                      as.character(res[[bad[1]]]))
m <- rbindlist(res)
if (nrow(m) != nrow(pair)) stop("rows ", nrow(m), " != joint pairs ", nrow(pair))
m <- merge(pair, m, by = c("exposure_id", "outcome_id"))
if (nrow(m) != nrow(pair)) stop("merge back onto the joint pairs lost rows")

setcolorder(m, c("fluid", "gene", "protein_assay", "disease"))
setorder(m, fluid, gene, disease)
fwrite(m, OUT)
if (!file.exists(OUT)) stop("Failed to write ", OUT)
print(m[, .N, by = .(fluid, multi_iv = n_ivs > 1, iv_match)][order(fluid, multi_iv, iv_match)])
message(sprintf("joint pairs %d | instrument %d | in_window %d | distinct %d | multi-IV with the top-weight IV colocalizing %d of %d",
                nrow(m), m[iv_match == "instrument", .N], m[iv_match == "in_window", .N],
                m[iv_match == "distinct", .N], m[n_ivs > 1 & top_iv_is_coloc_lead == TRUE, .N], m[n_ivs > 1, .N]))
message("Wrote ", OUT)
