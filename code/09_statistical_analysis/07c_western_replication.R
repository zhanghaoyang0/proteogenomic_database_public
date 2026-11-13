#!/usr/bin/env Rscript
# External replication of the CSF pQTL signals in Western et al. 2026 (Sci Transl Med,
# PMID 42054495): CSF Olink Explore HT, N = 2,524, US and Spanish cohorts, no overlap with
# BioFINDER. Full summary statistics, one GWAS Catalog study per analyte (07a lists them).
#
# Same rule as the BioFINDER-1 replication (Methods): the sentinel is the highest-PIP
# variant of the credible set; if Western lacks it, the next-highest-PIP variant of the
# set is tested. Replicated = same direction after allele harmonisation and P < 0.05.
# Variants are matched on rsID (Western files are GRCh38, ours hg19).
#
# Protein matching by UniProt: an analyte with exactly our accession set, else one sharing
# a subunit; among several, LOD pass, then CV pass, then the lowest Olink ID. GBP1 and
# MAP2K1 are three analytes each whose GWAS Catalog studies are named by gene only, so
# their study cannot be told apart: such proteins are "not testable: ambiguous".
#
# Western is also a 07b novelty reference. For a protein it reported a cis-pQTL for, a
# replicated cis signal is close to given; the informative subsets are the novel signals
# and those BioFINDER-1 could not test.
#
# Two passes. Pass 1 finds the studies not yet in LOOKUP_DIR, writes the download list and
# the rsIDs to keep to tmp/proteogenomic_database/ (repo root), and stops; run 07c_western_download.sh, then this again.
#   output/07c_western_replication.csv   one row per distinct CSF signal
# and folds the call into S1 (add_western, lib/common.R), as 07b does for novelty.
#
# Usage (repo root): Rscript code/09_statistical_analysis/07c_western_replication.R

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/cs_signals.R")   # distinct_signals
source("code/09_statistical_analysis/lib/common.R")       # add_western

OUT_DIR    <- "output"
TMP_DIR    <- "tmp/proteogenomic_database"
S1         <- file.path(OUT_DIR, "S1_pqtl_finemapping.csv.gz")
KEPT       <- file.path(OUT_DIR, "00_kept_assays.csv")
NOV7       <- file.path(OUT_DIR, "07_novelty_credible_sets.csv.gz")
WO_DIR     <- path.expand("~/db/shared_data/pqtl_ref/western_csf_olink")
WO_S4      <- file.path(WO_DIR, "table_s4_olink_analytes.csv")
WO_GCST    <- file.path(WO_DIR, "gwas_catalog_studies.csv")
LOOKUP_DIR <- file.path(WO_DIR, "lookup")
OUT        <- file.path(OUT_DIR, "07c_western_replication.csv")
P_REP      <- 0.05

for (f in c(S1, KEPT, NOV7, WO_S4, WO_GCST)) if (!file.exists(f)) stop("missing input: ", f)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- our distinct CSF signals and their credible-set variants ---------------------
s1 <- fread(S1)
kept <- fread(KEPT)
sig <- s1[is_cs_lead == TRUE & fluid == "CSF" &
            paste(cohort, protein_assay) %in% kept[, paste(cohort, protein_assay)],
          .(cs_key, cohort, protein_assay, cs_id, gene, is_cis, bf1_status, novelty)]
sig <- distinct_signals(sig)
cs <- s1[cs_key %in% sig$cs_key, .(cs_key, snp, ea, nea, beta, pip)]
setorder(cs, cs_key, -pip)
cs[, pip_rank := seq_len(.N), by = cs_key]

# ---- protein -> Western study --------------------------------------------------------
up <- unique(fread(NOV7, select = c("cs_key", "uniprot_id")))
n0 <- nrow(sig)
sig <- merge(sig, up, by = "cs_key", all.x = TRUE)
if (nrow(sig) != n0 || anyNA(sig$uniprot_id)) stop("CSF signals without exactly one UniProt set")
sig[, uniprot_id := norm_ids(uniprot_id)]

gc <- fread(WO_GCST)
wo <- fread(WO_S4)[Assay %in% gc$gene]   # the analysed panel, as in 07a
wo[, `:=`(uniprot_id = norm_ids(UniProt),
          ambiguous = Assay %in% gc[duplicated(gene), gene],
          qc_rank = (LOD != "PASS") + (CV_filter != "PASS"))]
setorder(wo, qc_rank, OlinkID)
pick <- function(u) {
  hit <- wo[uniprot_id == u]
  if (!nrow(hit)) hit <- wo[vapply(strsplit(uniprot_id, "_"),
                                   function(w) any(w %in% strsplit(u, "_")[[1]]), TRUE)]
  if (!nrow(hit)) return(data.table(w_olink_id = NA_character_, w_gene = NA_character_, ambiguous = FALSE))
  hit[1, .(w_olink_id = OlinkID, w_gene = Assay, ambiguous)]
}
map <- unique(sig[, .(uniprot_id)])[, pick(uniprot_id), by = uniprot_id]
map <- merge(map, unique(gc[!gene %in% gc[duplicated(gene), gene], .(w_gene = gene, accession, ftp_dir)]),
             by = "w_gene", all.x = TRUE)
sig <- merge(sig, map, by = "uniprot_id", all.x = TRUE)
if (nrow(sig) != n0) stop("study mapping changed the signal count")
message(sprintf("CSF signals: %d; protein analysed in Western: %d; ambiguous study: %d",
                nrow(sig), sig[!is.na(w_olink_id), .N], sig[ambiguous == TRUE, .N]))

# ---- pass 1: anything left to download? --------------------------------------------
acc <- unique(sig[!is.na(accession), .(accession, ftp_dir)])
acc[, file := file.path(LOOKUP_DIR, paste0(accession, ".tsv.gz"))]
todo <- acc[!file.exists(file)]
if (nrow(todo)) {
  fwrite(todo[, .(accession, url = sprintf("%s/%s.tsv.gz", ftp_dir, accession))],
         file.path(TMP_DIR, "07c_western_download.tsv"), sep = "\t", col.names = FALSE)
  writeLines(unique(cs[cs_key %in% sig[!is.na(accession), cs_key] & startsWith(snp, "rs"), snp]),
             file.path(TMP_DIR, "07c_western_rsids.txt"))
  message(nrow(todo), " of ", nrow(acc), " Western studies not yet in ", LOOKUP_DIR,
          "; lists written to ", TMP_DIR, ". Next:\n",
          "  sbatch code/09_statistical_analysis/07c_western_download.sh\n",
          "then run this script again.")
  quit(save = "no", status = 0)
}

# ---- pass 2: replication -----------------------------------------------------------
w <- rbindlist(lapply(seq_len(nrow(acc)), function(i)
  fread(acc$file[i], colClasses = list(character = "rs_id"))[, accession := acc$accession[i]]))
w <- w[, .(accession, snp = rs_id, w_ea = effect_allele, w_nea = other_allele,
           w_beta = beta, w_p = 10^-neg_log_10_p_value)]

# An rsID can label several alleles at one site; harmonising keeps the matching one.
tv <- merge(cs, sig[!is.na(accession), .(cs_key, accession)], by = "cs_key")
tv <- merge(tv, w, by = c("accession", "snp"), allow.cartesian = TRUE)
tv[, w_beta_h := fcase(w_ea == ea & w_nea == nea, w_beta,
                       w_ea == nea & w_nea == ea, -w_beta,
                       default = NA_real_)]
tv <- tv[!is.na(w_beta_h)]                     # allele mismatch: try the next variant
if (anyDuplicated(tv, by = c("cs_key", "snp"))) stop("a variant matches two Western rows with the same alleles")
tv <- tv[order(cs_key, pip_rank)][, .SD[1], by = cs_key]

# Smallest Western P over the whole credible set, alleles aside: 07b's novelty call reads
# P < 5e-8 here as "detected by Western but not listed as a sentinel".
wmin <- merge(merge(cs[, .(cs_key, snp)], sig[!is.na(accession), .(cs_key, accession)], by = "cs_key"),
              w[, .(accession, snp, w_p)], by = c("accession", "snp"), allow.cartesian = TRUE)[
  , .(western_cs_min_p = min(w_p)), by = cs_key]
res <- merge(sig, tv[, .(cs_key, tested_snp = snp, tested_pip_rank = pip_rank, beta_bf2 = beta,
                         beta_western = w_beta_h, p_western = w_p)], by = "cs_key", all.x = TRUE)
res <- merge(res, wmin, by = "cs_key", all.x = TRUE)
if (nrow(res) != n0) stop("replication merge changed the signal count")
res[, western_status := fcase(
  is.na(w_olink_id), "not testable: protein not analysed in Western",
  ambiguous == TRUE, "not testable: ambiguous Western study",
  is.na(tested_snp), "not testable: no credible-set variant in Western",
  sign(beta_western) == sign(beta_bf2) & p_western < P_REP, "replicated",
  default = "tested, not replicated")]

setorder(res, -is_cis, gene, cs_key)
fwrite(res[, .(cs_key, protein_assay, gene, uniprot_id, is_cis, novelty, bf1_status,
               western_gene = w_gene, western_olink_id = w_olink_id, accession, tested_snp,
               tested_pip_rank, beta_bf2, beta_western, p_western, western_status,
               western_cs_min_p)], OUT)
if (!file.exists(OUT)) stop("not written: ", OUT)
# na.strings: as in 07b, so an unchanged S1 is rewritten byte-identical
s1_new <- add_western(fread(S1, na.strings = c("", "NA")), OUT)
if (is.null(s1_new)) stop("07c credible sets missing from S1: rerun 07c after 06")
fwrite(s1_new, S1)

res[, `:=`(tier = sub(":.*", "", novelty), bf1_testable = !startsWith(bf1_status, "not"))]
print(dcast(res, is_cis + tier ~ western_status, value.var = "cs_key", fun.aggregate = length))
print(dcast(res, is_cis + bf1_testable ~ western_status, value.var = "cs_key", fun.aggregate = length))
message("07c_western_replication done -> ", OUT)
