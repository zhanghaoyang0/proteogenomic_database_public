#!/usr/bin/env Rscript
# Data for 01_pqtl_architecture.png and 09_supp_replication.png a. Analysis only: this script reads the
# raw pQTL and replication inputs and writes tidy per-panel tables; 01c_plot.R
# draws from those tables and never touches the raw inputs. The hotspot tables are
# prepared by 01b_hotspot_data.R.
# Outputs:
#   01_pqtl_pos_summary.csv        - lead / credible-set counts per assay (per platform)
#                                    (`set` column: "all assays" per assay, and
#                                    "deduped assays" on the deduped protein set)
#   01a_fig1_map.csv               - one row per plotted credible set (the map)
#   01a_fig1_chr_prop.csv          - CSF/plasma proportion per variant chromosome (the strip)
#   01a_fig1_replication.csv       - BF2 vs BF1 sentinel effect sizes (09_supp_replication.png a, lib/repl_panel.R)
# Usage:
#   Rscript code/09_statistical_analysis/01a_pqtl_data.R [out_dir]

suppressPackageStartupMessages(library(data.table))

DEFAULT_OUT <- "output"

args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[1] else DEFAULT_OUT
summary_csv <- file.path(OUT_DIR, "01_pqtl_pos_summary.csv")

COHORT <- "bf2"
INFO_DIR <- "info/bf_protein"
ENSG_CSV <- file.path(INFO_DIR, "ensg", "assay_ensg.csv")
CLUMP_CSV <- file.path(INFO_DIR, "clump_lead_pqtl.csv.gz")
FINEMAP_CSV <- file.path(INFO_DIR, "finemap_lead_snps.csv")
REPL_TSV <- file.path(INFO_DIR, "replication", "bf1_targeted_replication.tsv.gz")
CS_REPL_TSV <- file.path(INFO_DIR, "replication", "cs_replication.tsv.gz")

source("code/09_statistical_analysis/lib/cs_signals.R")
source("code/09_statistical_analysis/lib/common.R")   # assay_platform, assay_fluid, CHR_LEN_HG19

classify_panel <- function(tag) list(platform = assay_platform(tag), fluid = assay_fluid(tag))

annot <- fread(
  ENSG_CSV,
  select = c("cohort", "tag", "protein_assay", "gene_chr", "gene_start")
)
annot <- unique(annot[cohort == COHORT], by = c("cohort", "protein_assay"))
annot[, c("platform", "fluid") := classify_panel(tag)]
annot <- annot[!is.na(fluid)]

prepare_signals <- function(d, variant_chr, variant_pos) {
  d <- merge(d, annot, by = c("cohort", "protein_assay"), all.x = TRUE)
  setnames(d, c(variant_chr, variant_pos), c("variant_chr", "variant_pos"))
  d[, `:=`(
    platform = factor(platform, levels = c("Olink", "NuLISA", "SomaLogic")),
    fluid = factor(fluid, levels = c("CSF", "Plasma")),
    cis_trans = factor(ifelse(is_cis, "cis", "trans"), levels = c("cis", "trans"))
  )]
  d
}

lead <- fread(
  CLUMP_CSV,
  select = c("cohort", "protein_assay", "SNP", "CHR", "BP", "is_cis")
)
lead <- prepare_signals(lead[cohort == COHORT], "CHR", "BP")

finemap <- fread(
  FINEMAP_CSV,
  select = c(
    "cohort", "protein_assay", "chr", "pos", "snp", "cs_id",
    "cs_size", "cs_purity", "max_pip", "is_cis"
  )
)
# Overlapping +/-500 kb windows fine-map one signal more than once, so the raw row
# count is not a signal count (finemap/03_dedup_credible_sets.R).
finemap <- distinct_signals(finemap[cohort == COHORT])
finemap <- prepare_signals(finemap, "chr", "pos")

chr_len <- CHR_LEN_HG19
chromosomes <- names(chr_len)
autosomes <- as.character(1:22)
chr_offset <- c(0, cumsum(chr_len[-length(chr_len)]))
names(chr_offset) <- chromosomes

to_linear <- function(chr, pos) {
  chr_offset[as.character(chr)] + as.numeric(pos) / 1e6
}

add_plot_coordinates <- function(d) {
  d[, plot_ok := as.character(variant_chr) %in% autosomes &
    as.character(gene_chr) %in% autosomes &
    !is.na(variant_pos) & !is.na(gene_start)]
  d[, `:=`(
    variant_linear = to_linear(variant_chr, variant_pos),
    gene_linear = to_linear(gene_chr, gene_start)
  )]
  d
}

lead <- add_plot_coordinates(lead)
finemap <- add_plot_coordinates(finemap)


# --- Per-panel tables ----------------------------------------------------------
# The figure uses the deduped analyzed set: one highest-N assay per protein per
# fluid (SomaLogic excluded, N>1000; 00_dedup_assays.R), so each protein appears
# once per fluid. The per-platform summary above stays per-assay on purpose.
kept_assays <- fread(file.path(OUT_DIR, "00_kept_assays.csv"))
finemap_fig <- finemap[paste(cohort, protein_assay) %in%
                       kept_assays[, paste(cohort, protein_assay)]]

# Panel a: one row per plotted credible set.
map_d <- finemap_fig[plot_ok == TRUE & as.character(platform) != "SomaLogic",
                     .(protein_assay, fluid, cis_trans, variant_chr, variant_pos,
                       gene_chr, gene_start, variant_linear, gene_linear)]
fwrite(map_d, file.path(OUT_DIR, "01a_fig1_map.csv"))
message("Saved: 01a_fig1_map.csv (", nrow(map_d), " credible sets)")

# Panel b: CSF/plasma proportion within each variant chromosome.
prop_d <- map_d[as.character(variant_chr) %in% autosomes,
                .(n = .N), by = .(variant_chr, fluid)]
prop_d[, prop := n / sum(n), by = variant_chr]
fwrite(prop_d, file.path(OUT_DIR, "01a_fig1_chr_prop.csv"))

# Panel d: BF2 vs BF1 sentinel effect sizes on two nested test sets.
#   exact assay - same assay in both cohorts (essentially all of plasma:
#                 ExploreHT_o7plasma on both sides).
#   cross-assay - same gene, same fluid, different Olink assay version (BF2 Explore
#                 3072 *_o2csf vs BF1 ExploreHT/Target-96). Only 16 CSF sentinels
#                 are testable on exact assays, so CSF replication is only visible
#                 once these are included; they replicate the genetic signal rather
#                 than an identical analyte.
# Replication is called at a concordant effect direction plus BF1 P < 0.05. These are
# 2,664 sentinels pre-specified in BF2, not a genome-wide scan, so a Bonferroni
# correction over the whole set is the wrong null: it removed almost all the trans
# signals (CSF trans 150 -> 72, plasma trans 336 -> 170), which are small effects in
# the smaller cohort, and so scored low power as failure to replicate.
r_all <- fread(REPL_TSV)
# The replication step submitted one row per credible set, so the window-duplicate
# signals inflate this funnel too. The dedup is on the sentinel's LEAD variant, not
# on the row's own cs_id: a submitted sentinel can stand for two credible sets that
# share a lead, and only one of the two carries the row's cs_id, so a cs_id filter
# would drop sentinels whose signal is in fact represented.
#
# That lead cannot be read off `SNP`. When the sentinel is absent from the BF1
# sumstats, 03 substitutes the highest-PIP credible-set member and writes the PROXY
# there, so matching `SNP` against the lead set silently threw away all 122
# proxy-rescued sentinels -- out of the numerator AND the denominator, which is why
# the funnel went on printing the pre-proxy numbers and looked reassuringly stable.
# The lead is recovered through cs_id from cs_replication.tsv.gz, which is also
# where is_representative now comes from, so this script and 04 partition the
# sentinels the same way by construction.
cs_rep <- fread(CS_REPL_TSV, select = c("trait", "cs_id", "cs_lead_snp",
                                        "is_representative"))
r_all[, trait_key := paste(cohort, bf2_assay, sep = "__")]
n_submitted_raw <- nrow(r_all)
r_all <- merge(r_all, unique(cs_rep[, .(trait_key = trait, cs_id, cs_lead_snp)]),
               by = c("trait_key", "cs_id"))
rep_leads <- unique(cs_rep[is_representative == TRUE,
                           .(trait_key = trait, cs_lead_snp)])
# Reported signals only: trans sentinels below the study-wide threshold leave the
# funnel (lib/cs_signals.R).
rep_leads <- merge(rep_leads, reported_leads()[, .(trait_key = trait, cs_lead_snp)],
                   by = c("trait_key", "cs_lead_snp"))
r_all <- merge(r_all, rep_leads, by = c("trait_key", "cs_lead_snp"))
# The analyzed assays only, as in S1 and the figure panels: replication was submitted
# for every BF2 assay, including ones 00_dedup_assays.R dropped for a second assay of
# the same protein and fluid (NPPB plasma, kept as NT_proBNP; code review 2026-10-30, R12).
n_dedup <- nrow(r_all)
r_all <- r_all[paste(cohort, bf2_assay) %in% kept_assays[, paste(cohort, protein_assay)]]
message("submitted sentinels: ", n_submitted_raw, " -> ", n_dedup,
        " after credible-set dedup -> ", nrow(r_all), " on analyzed assays")

# ---- Replication funnel: how much of the discovery set the test actually covers ----
# Table 1 of the manuscript. Two things remove a sentinel, and they are different
# statements about the cohort, so they are counted separately rather than pooled into
# one "not testable" number: the protein is not on the BioFINDER-1 panel for that
# fluid, or it is but the sentinel variant did not survive BioFINDER-1 QC.
r_all[, fluid_lab := fifelse(fluid == "csf", "CSF", "Plasma")]
r_all[, cis_lab := fifelse(is_cis == TRUE, "cis", "trans")]
# The verdict is read, not recomputed. code/05_replication/04_credible_set_status.R
# owns the rule (cross-assay counts as testable; replication is a concordant
# direction plus BF1 P < 0.05) and stops if the funnel below moves. This script,
# 01b and 06 each used to rebuild it, so the three could drift apart silently.
# Only the flags come from there -- the BF2 effect sizes and frequencies the
# panels and the attrition test need are still read from the raw table above.
status <- unique(fread(CS_REPL_TSV, select = c("trait", "cs_lead_snp", "bf1_status",
                                               "direction_concordant")),
                 by = c("trait", "cs_lead_snp"))
# The raw table carries a column of the same name, so rename before merging.
setnames(status, "direction_concordant", "cs_concordant")
# On cs_lead_snp, which was recovered above -- not on SNP, which may be a proxy.
r_all <- merge(r_all, status, by.x = c("trait_key", "cs_lead_snp"),
               by.y = c("trait", "cs_lead_snp"), all.x = TRUE)
if (anyNA(r_all$bf1_status)) {
  stop(sum(is.na(r_all$bf1_status)), " submitted sentinels have no verdict in ",
       basename(CS_REPL_TSV), " -- rerun code/05_replication/04_credible_set_status.R")
}
r_all[, measured := assay_match_type != "not_measured_same_fluid"]
r_all[, testable := bf1_status %in% c("replicated", "tested, not replicated")]
r_all[, concordant := cs_concordant]
r_all[, passed := bf1_status == "replicated"]

funnel_cols <- function(d) d[, .(
  submitted = .N,
  measured_in_bf1 = sum(measured),
  testable = sum(testable),
  concordant = sum(concordant & testable, na.rm = TRUE),
  replicated = sum(passed & testable, na.rm = TRUE)
)]
funnel <- rbind(
  r_all[, funnel_cols(.SD), by = .(fluid_lab, cis_lab)],
  cbind(fluid_lab = "All", cis_lab = "", r_all[, funnel_cols(.SD)])
)
fwrite(funnel, file.path(OUT_DIR, "01a_replication_funnel.csv"))
message("Replication funnel (Table 1):")
print(funnel)

# The attrition is not neutral, and the manuscript says so, so the check that
# establishes it is run here rather than by hand. Untested sentinels are lower
# frequency in CSF and weaker in plasma -- both properties that make a signal harder
# to recover in a smaller cohort -- so the replication rate is a mild over-estimate.
r_all[, maf := pmin(bf2_eaf, 1 - bf2_eaf)]
r_all[, neglog10p := -log10(pmax(bf2_p, 1e-300))]
r_all[, drop_reason := fcase(
  testable == TRUE, "testable",
  assay_match_type == "not_measured_same_fluid", "protein not measured in BF1",
  default = "sentinel absent from BF1 sumstats"
)]
attrition <- r_all[is.finite(maf), .(
  n = .N, median_maf = round(median(maf), 3),
  median_neglog10p = round(median(neglog10p), 1)
), by = .(fluid_lab, drop_reason)][order(fluid_lab, drop_reason)]
attrition_test <- rbindlist(lapply(c("CSF", "Plasma"), function(f) {
  d <- r_all[fluid_lab == f & is.finite(maf)]
  data.table(
    fluid_lab = f,
    maf_wilcox_p = signif(wilcox.test(maf ~ testable, data = d)$p.value, 3),
    neglog10p_wilcox_p = signif(wilcox.test(neglog10p ~ testable, data = d)$p.value, 3)
  )
}))
fwrite(merge(attrition, attrition_test, by = "fluid_lab", all.x = TRUE),
       file.path(OUT_DIR, "01a_replication_attrition.csv"))
message("Replication attrition (untested sentinels vs testable):")
print(attrition); print(attrition_test)

r <- r_all[testable == TRUE]
fwrite(r[, .(replication_id, fluid_lab, cis_lab, assay_match_type,
             bf2_beta, bf1_beta, bf1_p, concordant, passed)],
       file.path(OUT_DIR, "01a_fig1_replication.csv"))
message("Saved: 01a_fig1_replication.csv (", nrow(r), " testable sentinels, ",
        sum(r$passed), " replicated at P < 0.05)")

# Deduped discovery counts. These, not the per-assay table above, are the numbers
# quoted in the manuscript: one protein contributes once per fluid.
dedup_summary <- function(d, signal_label) {
  # d already carries platform/fluid from classify_panel(); take only the gene key.
  x <- merge(
    d[paste(cohort, protein_assay) %in% kept_assays[, paste(cohort, protein_assay)]],
    kept_assays[, .(cohort, protein_assay, ensembl_gene_id)],
    by = c("cohort", "protein_assay")
  )
  by_platform <- x[, .(n_signals = .N, n_cis = sum(is_cis), n_trans = sum(!is_cis),
                       n_proteins = uniqueN(ensembl_gene_id)), by = .(fluid, platform)]
  by_fluid <- x[, .(platform = "All", n_signals = .N, n_cis = sum(is_cis),
                    n_trans = sum(!is_cis), n_proteins = uniqueN(ensembl_gene_id)),
                by = .(fluid)]
  out <- rbind(by_platform, by_fluid, use.names = TRUE)
  out[, signal_type := signal_label]
  out[order(fluid, platform)]
}

dedup <- rbind(dedup_summary(lead, "clumped lead"),
               dedup_summary(finemap, "95% credible set"))
assayed <- rbind(
  kept_assays[, .(n_proteins = uniqueN(ensembl_gene_id)), by = .(fluid, platform)],
  kept_assays[, .(platform = "All", n_proteins = uniqueN(ensembl_gene_id)), by = .(fluid)],
  use.names = TRUE
)
assayed[, `:=`(signal_type = "assayed", n_signals = NA_integer_,
               n_cis = NA_integer_, n_trans = NA_integer_)]

dedup_out <- rbind(assayed, dedup, use.names = TRUE, fill = TRUE)[
  order(signal_type, fluid, platform)]
dedup_out[, set := "deduped assays"]
# Deduped analyzed assays only; the per-assay "all assays" half had no reader and was
# dropped (code review 2026-10-30, R17). `set` stays, as readers filter on it.
fwrite(dedup_out, summary_csv)
message("Saved summary: ", summary_csv, " (proteins assayed in both fluids: ",
        length(intersect(kept_assays[fluid == "CSF", ensembl_gene_id],
                         kept_assays[fluid == "Plasma", ensembl_gene_id])), ")")
