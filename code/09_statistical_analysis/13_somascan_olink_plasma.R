#!/usr/bin/env Rscript
# Secondary, within-BioFINDER-2 plasma Olink--SomaScan pQTL comparison.
# This does not enter the primary Olink/NuLISA discovery set.
#
# Usage:
#   Rscript code/09_statistical_analysis/13_somascan_olink_plasma.R [out_dir]
#
# Outputs:
#   13_somascan_olink_summary.csv          compact target, signal and direction summary
#   13_somascan_olink_targets.csv          one-to-one UniProt target mapping and assay N
#   13_somascan_olink_signals.csv.gz       Olink lead-level comparison details
# Fig. 2c is drawn from the signals table by 01c_plot.R.

suppressPackageStartupMessages({
  library(data.table)
})

DEFAULT_OUT <- "output"
args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args)) args[[1L]] else DEFAULT_OUT

ANNOT <- "info/bf_protein/assay_annot.csv"
CLUMP <- "info/bf_protein/clump_lead_pqtl.csv.gz"
FINEMAP <- "info/bf_protein/finemap_lead_snps.csv"
P_THRESHOLD <- 5e-8
N_PRIMARY_ASSAYS <- 3538L
TRANS_THRESHOLD <- P_THRESHOLD / N_PRIMARY_ASSAYS
LOCUS_WINDOW_BP <- 1e6

for (f in c(ANNOT, CLUMP, FINEMAP)) if (!file.exists(f)) stop("missing input: ", f)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

platform <- function(tag) fcase(
  grepl("somalogic", tolower(tag)), "SomaScan",
  grepl("nulisa", tolower(tag)), "NuLISA",
  default = "Olink"
)
fluid <- function(tag) fcase(
  grepl("csf", tolower(tag)), "CSF",
  grepl("plasma", tolower(tag)), "Plasma",
  default = NA_character_
)

# Only a single, unambiguous UniProt accession is eligible. A repeated target on
# either platform is excluded rather than arbitrarily selecting one of its assays.
single_uniprot <- function(x) {
  id <- trimws(as.character(x))
  id[is.na(id) | id == "" | grepl("[_;|]", id)] <- NA_character_
  id
}

annot <- unique(fread(ANNOT), by = c("cohort", "protein_assay"))
annot[, `:=`(platform = platform(tag), fluid = fluid(tag),
             uniprot_id = single_uniprot(uniprot_id), N = as.numeric(N))]
annot <- annot[cohort == "bf2" & fluid == "Plasma" &
                 platform %in% c("Olink", "SomaScan") & !is.na(uniprot_id)]
if (!nrow(annot)) stop("no eligible BF2 plasma Olink/SomaScan assays in ", ANNOT)

counts <- annot[, .N, by = .(platform, uniprot_id)]
eligible_ids <- counts[, .(one_assay = all(N == 1L), platforms = uniqueN(platform)),
                       by = uniprot_id][one_assay & platforms == 2L, uniprot_id]
target <- annot[uniprot_id %in% eligible_ids,
                .(uniprot_id, platform, protein_assay, gene, ensembl_gene_id, assay_N = N)]
target <- dcast(target, uniprot_id ~ platform,
                value.var = c("protein_assay", "gene", "ensembl_gene_id", "assay_N"))
target_cols <- c("protein_assay_Olink", "protein_assay_SomaScan",
                 "gene_Olink", "gene_SomaScan",
                 "ensembl_gene_id_Olink", "ensembl_gene_id_SomaScan",
                 "assay_N_Olink", "assay_N_SomaScan")
if (!setequal(setdiff(names(target), "uniprot_id"), target_cols))
  stop("unexpected one-to-one target columns after platform pivot")
setcolorder(target, c("uniprot_id", target_cols))
if (!nrow(target)) stop("no one-to-one UniProt targets between Olink and SomaScan")
if (anyNA(target$protein_assay_Olink) | anyNA(target$protein_assay_SomaScan))
  stop("one-to-one target mapping is incomplete")

lead <- fread(CLUMP, select = c("cohort", "protein_assay", "SNP", "CHR", "BP", "P", "is_cis"))
lead <- lead[cohort == "bf2" & protein_assay %in%
               c(target$protein_assay_Olink, target$protein_assay_SomaScan) &
               (is_cis | P <= TRANS_THRESHOLD)]
if (!nrow(lead)) stop("no genome-wide significant leads for matched targets")
if (anyNA(lead$SNP) | anyNA(lead$CHR) | anyNA(lead$BP) | anyNA(lead$is_cis))
  stop("clumped lead table has missing signal keys")
if (anyDuplicated(lead, by = c("protein_assay", "SNP", "is_cis")))
  stop("duplicate clumped lead key within an assay")

# Effect directions are only compared for an exact sentinel. Different lead variants
# may tag the same locus with unknown LD phase, so their signed effects are not compared.
fm <- fread(FINEMAP, select = c("cohort", "protein_assay", "snp", "ea", "nea", "beta", "max_pip"))
fm <- fm[cohort == "bf2" & protein_assay %in% lead$protein_assay & snp %in% lead$SNP]
setorder(fm, protein_assay, snp, -max_pip)
fm <- fm[, .SD[1L], by = .(protein_assay, snp)]
setnames(fm, c("snp", "ea", "nea", "beta"),
         c("SNP", "effect_allele", "other_allele", "beta"))

olink <- merge(target[, .(uniprot_id, olink_assay = protein_assay_Olink)],
               lead[protein_assay %in% target$protein_assay_Olink],
               by.x = "olink_assay", by.y = "protein_assay", allow.cartesian = TRUE)
soma <- merge(target[, .(uniprot_id, soma_assay = protein_assay_SomaScan)],
              lead[protein_assay %in% target$protein_assay_SomaScan],
              by.x = "soma_assay", by.y = "protein_assay", allow.cartesian = TRUE)
if (!nrow(olink)) stop("no Olink leads among one-to-one targets")

setnames(olink, c("SNP", "CHR", "BP", "P", "is_cis"),
         c("olink_snp", "olink_chr", "olink_bp", "olink_p", "is_cis"))
setnames(soma, c("SNP", "CHR", "BP", "P", "is_cis"),
         c("soma_snp", "soma_chr", "soma_bp", "soma_p", "soma_is_cis"))

pair <- merge(olink, soma, by = "uniprot_id", allow.cartesian = TRUE)
pair <- pair[is_cis == soma_is_cis]
pair[, `:=`(same_lead = olink_snp == soma_snp,
            same_locus = olink_chr == soma_chr & abs(olink_bp - soma_bp) <= LOCUS_WINDOW_BP)]

signal <- copy(olink)
signal[, signal_key := paste(uniprot_id, olink_assay, olink_snp, is_cis, sep = "|")]
if (nrow(pair)) {
  overlap <- pair[, .(soma_signal = .N, soma_same_locus = any(same_locus),
                       soma_same_lead = any(same_lead)),
                  by = .(uniprot_id, olink_assay, olink_snp, is_cis)]
  signal <- merge(signal, overlap,
                  by = c("uniprot_id", "olink_assay", "olink_snp", "is_cis"), all.x = TRUE)
}
signal[is.na(soma_signal), `:=`(soma_signal = 0L, soma_same_locus = FALSE, soma_same_lead = FALSE)]
signal[, `:=`(comparison = fcase(
  soma_same_lead, "exact lead",
  soma_same_locus, "same locus",
  soma_signal > 0L, "different locus",
  default = "no SomaScan signal"
), signal_type = fifelse(is_cis, "cis", "trans"))]

exact <- pair[same_lead == TRUE]
if (nrow(exact)) {
  exact <- merge(exact, fm[, .(protein_assay, SNP, olink_ea = effect_allele,
                               olink_nea = other_allele, olink_beta = beta)],
                 by.x = c("olink_assay", "olink_snp"),
                 by.y = c("protein_assay", "SNP"), all.x = TRUE)
  exact <- merge(exact, fm[, .(protein_assay, SNP, soma_ea = effect_allele,
                               soma_nea = other_allele, soma_beta = beta)],
                 by.x = c("soma_assay", "soma_snp"),
                 by.y = c("protein_assay", "SNP"), all.x = TRUE)
  exact[, allele_relation := fcase(
    olink_ea == soma_ea & olink_nea == soma_nea, "same",
    olink_ea == soma_nea & olink_nea == soma_ea, "swapped",
    default = NA_character_
  )]
  exact <- exact[!is.na(allele_relation) & is.finite(olink_beta) & is.finite(soma_beta)]
  exact[, soma_beta_aligned := fifelse(allele_relation == "same", soma_beta, -soma_beta)]
  exact[, concordant_direction := sign(olink_beta) == sign(soma_beta_aligned)]
  direction <- exact[, .(uniprot_id, olink_assay, olink_snp,
                          is_cis, olink_beta, soma_beta, soma_beta_aligned,
                          allele_relation, concordant_direction)]
  signal <- merge(signal, direction,
                  by = c("uniprot_id", "olink_assay", "olink_snp", "is_cis"),
                  all.x = TRUE)
} else {
  signal[, `:=`(olink_beta = NA_real_, soma_beta = NA_real_,
                soma_beta_aligned = NA_real_, allele_relation = NA_character_,
                concordant_direction = NA)]
}
if (nrow(signal) != uniqueN(signal$signal_key)) stop("comparison changed Olink signal rows")

summary <- rbindlist(list(
  data.table(section = "targets", metric = "one-to-one UniProt targets", signal_type = "all",
             numerator = nrow(target), denominator = nrow(target), percent = 100),
  data.table(section = "availability", metric = "median assay N", signal_type = "Olink",
             numerator = median(target$assay_N_Olink), denominator = NA_real_, percent = NA_real_),
  data.table(section = "availability", metric = "median assay N", signal_type = "SomaScan",
             numerator = median(target$assay_N_SomaScan), denominator = NA_real_, percent = NA_real_),
  signal[, .(numerator = .N, denominator = .N, percent = 100),
         by = .(signal_type)][, `:=`(section = "signals", metric = "Olink leads")],
  signal[, .(numerator = sum(soma_signal > 0L), denominator = .N,
             percent = 100 * mean(soma_signal > 0L)),
         by = signal_type][, `:=`(section = "signal overlap", metric = "SomaScan signal of same cis/trans type")],
  signal[, .(numerator = sum(soma_same_locus), denominator = .N,
             percent = 100 * mean(soma_same_locus)),
         by = signal_type][, `:=`(section = "locus overlap", metric = "SomaScan lead within 1 Mb")],
  signal[, .(numerator = sum(soma_same_lead), denominator = .N,
             percent = 100 * mean(soma_same_lead)),
         by = signal_type][, `:=`(section = "lead overlap", metric = "identical sentinel rsID")],
  signal[!is.na(concordant_direction),
         .(numerator = sum(concordant_direction), denominator = .N,
           percent = 100 * mean(concordant_direction)),
         by = signal_type][, `:=`(section = "effect direction", metric = "concordant exact-lead effect")]
), use.names = TRUE, fill = TRUE)
setcolorder(summary, c("section", "metric", "signal_type", "numerator", "denominator", "percent"))

scatter <- signal[!is.na(concordant_direction)]

fwrite(summary, file.path(OUT_DIR, "13_somascan_olink_summary.csv"))
fwrite(target, file.path(OUT_DIR, "13_somascan_olink_targets.csv"))
fwrite(signal, file.path(OUT_DIR, "13_somascan_olink_signals.csv.gz"))

stopifnot(nrow(target) > 0, nrow(signal) > 0, nrow(summary) > 0,
          all(summary[section != "availability", denominator > 0]))
message("Matched one-to-one targets: ", nrow(target),
        "; Olink leads: ", nrow(signal),
        "; exact-lead effects harmonized: ", nrow(scatter))
message("Wrote 13_somascan_olink_{summary,targets,signals}.csv[.gz] to ", OUT_DIR)
