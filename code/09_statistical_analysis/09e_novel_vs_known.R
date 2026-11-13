#!/usr/bin/env Rscript
# Novel versus previously reported CSF cis-pQTLs (docs/pipeline.md): are the novel
# signals weaker, rarer, more artefact-prone or less useful than the reported ones?
#
# Unit: the manuscript's distinct signals (lib/cs_signals.R), CSF cis, outside the MHC.
# Signal-level rows (genetic association characteristics) compare primary cis
# associations (09a assoc_class): novel = unmatched, A or B (the 09a cis candidates),
# known = matched. Secondary signals, unresolved, C and B detected by Western are out.
# Protein-level rows (protein characteristics) compare proteins by 07b's protein-level
# call: novel = no reported cis-pQTL (A), known = a reported one; C is out. B proteins
# have a reported cis-pQTL, so they are known here.
# MAF comes from 10a (the paired sample), so it covers proteins assayed in both fluids.
#
# Outputs: output/09e_novel_vs_known.csv, one row per metric, with the count (k) and
# denominator (n) for proportions or the median (mean for pav_mean) for continuous
# metrics, and a Fisher or Wilcoxon P; output/09e_platform_coverage.csv, every CSF protein
# with a BioFINDER cis-pQTL by which CSF reference reported one and which measured it
# (09f draws the first as panel a).
#
# Coverage: 07b's protein-level call per protein, split by ref_protein_cis_by: reported by
# both references (Cruchaga SomaScan, Western Olink), by one, by neither (tier A, the
# novel proteins) or measured by neither (tier C). The counts match the text's.
#
# Usage: Rscript code/09_statistical_analysis/09e_novel_vs_known.R

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/common.R")      # csq_class, COLOC_PP4
source("code/09_statistical_analysis/lib/cs_signals.R")  # distinct_signals

OUT_DIR <- "output"
S1   <- file.path(OUT_DIR, "S1_pqtl_finemapping.csv.gz")
CAND <- file.path(OUT_DIR, "09a_novel_candidates.csv")
NOV7 <- file.path(OUT_DIR, "07_novelty_credible_sets.csv.gz")
HET  <- file.path(OUT_DIR, "10a_fluid_heterogeneity.csv")
TIS  <- file.path(OUT_DIR, "04_protein_tissue_pair_evidence.csv.gz")
FEAT <- file.path(OUT_DIR, "02_uniprot_feature_data.csv.gz")
DRUG <- "info/bf_protein/drug_targets.tsv.gz"
VEP  <- file.path(GWAS_BASE, "bf2/annotation/protein/variant_annotation.tsv.gz")
MHC  <- "code/lib/mhc_hg19.range"
OUT  <- file.path(OUT_DIR, "09e_novel_vs_known.csv")
OUT_COV <- file.path(OUT_DIR, "09e_platform_coverage.csv")

for (f in c(S1, CAND, NOV7, HET, TIS, FEAT, DRUG, VEP, MHC)) if (!file.exists(f)) stop("missing input: ", f)

# ---- signals ------------------------------------------------------------------
s1 <- fread(S1)
mhc <- fread(MHC, header = FALSE, col.names = c("chr", "start", "end", "label"))
cis <- distinct_signals(s1[is_cs_lead == TRUE])
cis <- cis[fluid == "CSF" & is_cis == TRUE &
             !(chr == mhc$chr & pos >= mhc$start & pos <= mhc$end)]
n0 <- nrow(cis)
cis <- merge(cis, fread(CAND, select = c("cs_key", "assoc_class", "is_candidate")), by = "cs_key")
if (nrow(cis) != n0) stop("09a does not list every CSF cis signal")
cis[, pgrp := fcase(startsWith(novelty, "A"), "novel", startsWith(novelty, "known"), "known")]
if (cis[!is.na(pgrp), uniqueN(pgrp), by = protein_assay][V1 > 1, .N]) stop("a protein has both novel and known cis signals")
sig <- cis[, grp := fcase(assoc_class %in% c("A", "B"), "novel", assoc_class == "known", "known")][!is.na(grp)]
if (!setequal(sig[grp == "novel", cs_key], cis[is_candidate == TRUE, cs_key]))
  stop("novel primary associations differ from the 09a cis candidates")
if (anyDuplicated(sig$protein_assay)) stop("more than one primary association per protein")
# Novel means measured by the reference without a reported cis-pQTL (tier A), not
# missing from its panel (tier C): the text calls these platform differences.
ref_assayed <- fread(NOV7, select = c("cs_key", "ref_assayed"))
if (!all(ref_assayed[cs_key %chin% c(sig$cs_key, cis[!is.na(pgrp), cs_key]), ref_assayed]))
  stop("a novel or known signal's protein was not assayed by the reference")
message("primary associations: ", sig[grp == "novel", .N], " unmatched (",
        sig[assoc_class == "B", .N], " B), ", sig[grp == "known", .N], " matched; proteins: ",
        cis[pgrp == "novel", uniqueN(protein_assay)], " novel, ", cis[pgrp == "known", uniqueN(protein_assay)], " known")

# ---- tests ----------------------------------------------------------------------
GRP <- c("novel", "known")
prop_row <- function(metric, level, x, g) {
  t <- table(factor(g, GRP), factor(x, c(TRUE, FALSE)))
  data.table(metric, level, novel_k = t["novel", "TRUE"], novel_n = sum(t["novel", ]),
             known_k = t["known", "TRUE"], known_n = sum(t["known", ]),
             novel_value = t["novel", "TRUE"] / sum(t["novel", ]),
             known_value = t["known", "TRUE"] / sum(t["known", ]),
             p = fisher.test(t)$p.value, test = "Fisher")
}
cont_row <- function(metric, level, x, g, stat = median) {
  a <- x[g == "novel" & !is.na(x)]; b <- x[g == "known" & !is.na(x)]
  data.table(metric, level, novel_k = NA_integer_, novel_n = length(a),
             known_k = NA_integer_, known_n = length(b),
             novel_value = stat(a), known_value = stat(b),
             p = wilcox.test(a, b)$p.value, test = "Wilcoxon")
}

out <- list()
# signal level
out$beta <- cont_row("abs_beta", "signal", abs(sig$beta), sig$grp)
out$p    <- cont_row("neg_log10_p", "signal", -log10(sig$p), sig$grp)
out$cs   <- cont_row("cs_size", "signal", sig$cs_size, sig$grp)
het <- fread(HET, select = c("cs_key", "eaf"))
m <- merge(sig[, .(cs_key, grp)], het, by = "cs_key")
out$maf <- cont_row("maf_paired", "signal", pmin(m$eaf, 1 - m$eaf), m$grp)
te <- sig[!startsWith(bf1_status, "not")]
out$rep <- prop_row("bf1_replicated", "signal", te$bf1_status == "replicated", te$grp)

# protein-altering share of each signal's PIP (03b's definition)
vep <- unique(fread(VEP, select = c("SNP", "Consequence")), by = "SNP")
mem <- s1[cs_key %chin% sig$cs_key, .(cs_key, snp, pip)]
mem <- merge(mem, vep[, .(snp = SNP, csq = csq_class(Consequence))], by = "snp", all.x = TRUE)
miss <- mem[is.na(csq), sum(pip)] / mem[, sum(pip)]
if (miss > 0.01) stop(sprintf("%.1f%% of PIP mass has no VEP record", 100 * miss))
mem[, w := pip / sum(pip), by = cs_key]
pav <- merge(mem[, .(pav = sum(w[csq %in% "Protein-altering"])), by = cs_key],
             sig[, .(cs_key, grp)], by = "cs_key")
if (nrow(pav) != nrow(sig)) stop("protein-altering share missing for some signals")
out$pav <- cont_row("pav_mean", "signal", pav$pav, pav$grp, stat = mean)

# protein level
prot <- cis[!is.na(pgrp), .(grp = pgrp[1], ensg = ensembl_gene_id[1],
                            joint = any(n_diseases_joint > 0, na.rm = TRUE)), by = protein_assay]
feat <- merge(prot[, .(protein_assay, grp)],
              fread(FEAT)[fluid == "CSF", .(protein_assay, signal_peptide, has_glycosylation,
                                             is_membrane, is_secreted)],
              by = "protein_assay")
if (nrow(feat) < 0.95 * nrow(prot)) stop("UniProt features missing for over 5% of proteins")
for (v in c("signal_peptide", "has_glycosylation", "is_membrane", "is_secreted"))
  out[[v]] <- prop_row(v, "protein", feat[[v]] == 1, feat$grp)

# brain eQTL coloc of the protein's own gene, among proteins where one ran (09d's rule)
bc <- brain_coloc_sets(TIS)   # lib/common.R, shared with 09d
tested <- bc$tested
hit <- bc$hit
pb <- prot[protein_assay %chin% tested]
out$brain <- prop_row("brain_eqtl_coloc", "protein", pb$protein_assay %chin% hit, pb$grp)

out$joint <- prop_row("disease_mr_coloc", "protein", prot$joint, prot$grp)
drug <- fread(DRUG)
out$drug <- prop_row("drug_target", "protein", prot$ensg %chin% drug$ensembl_gene_id, prot$grp)
out$appr <- prop_row("approved_drug_target", "protein",
                     prot$ensg %chin% drug[stage_label == "Approved", ensembl_gene_id], prot$grp)

res <- rbindlist(out)
fwrite(res, OUT)

# ---- reference coverage: every CSF protein with a BioFINDER cis-pQTL ----------------
all_cis <- distinct_signals(s1[is_cs_lead == TRUE])[fluid == "CSF" & is_cis == TRUE]
all_cis[, tier := sub(":.*", "", novelty)]
if (!all(all_cis$tier %in% c("known", "A", "C"))) stop("unexpected novelty tier: ", paste(setdiff(all_cis$tier, c("known", "A", "C")), collapse = ", "))
n0 <- nrow(all_cis)
# na.strings: a protein no reference reported has an empty field, which is NA here
all_cis <- merge(all_cis, fread(NOV7, select = c("cs_key", "ref_assayed_by", "ref_protein_cis_by"), na.strings = c("", "NA")),
                 by = "cs_key", all.x = TRUE)
if (nrow(all_cis) != n0) stop("07 novelty merge changed the signal count")
prot <- unique(all_cis[, .(protein_assay, tier, measured_by = fcoalesce(ref_assayed_by, "neither"), ref_protein_cis_by)])
if (anyDuplicated(prot$protein_assay)) stop("a protein has cis signals with different reference calls")
if (prot[(tier == "known") != !is.na(ref_protein_cis_by), .N]) stop("tier and ref_protein_cis_by disagree")
COV_GROUPS <- c("Both studies", "Olink study only", "SomaScan study only", "BioFINDER only", "Not measured in either")
prot[, group := fcase(tier == "C", "Not measured in either",
                      tier == "A", "BioFINDER only",
                      ref_protein_cis_by == "Olink;SomaScan", "Both studies",
                      ref_protein_cis_by == "Olink", "Olink study only",
                      ref_protein_cis_by == "SomaScan", "SomaScan study only")]
if (anyNA(prot$group)) stop("unexpected ref_protein_cis_by: ", paste(unique(prot[is.na(group), ref_protein_cis_by]), collapse = ", "))
# measured_by: which reference measured the protein (Olink, SomaScan, both, neither), so the
# text can compare the platforms on the proteins both measured.
cov_out <- prot[, .N, by = .(group, measured_by)]
if (!all(cov_out$group %in% COV_GROUPS)) stop("unexpected coverage group")
cov_out <- cov_out[order(match(group, COV_GROUPS), measured_by)]
fwrite(cov_out, OUT_COV)
print(cov_out)
print(res[, .(metric, novel = sprintf("%.3g (%s)", novel_value, novel_n),
              known = sprintf("%.3g (%s)", known_value, known_n), p = signif(p, 2))])
message("wrote ", OUT, ", ", OUT_COV)
