#!/usr/bin/env Rscript
# Evidence table for the novel CSF pQTL signals (docs/pipeline.md):
# one row per 09a candidate, joining FinnGen MR/coloc (S2), brain eQTL coloc (04) and
# drug targets. (The Mac phenotype results 09b/09c were archived on 2026-10-31.)
# No script reads the S2 and drug columns; they stay on purpose, because the table is
# the readable evidence summary of the novel signals (review 20261030 R17).
#
# Grain of each block, which is not always the signal:
#   - S2 and drug targets are per protein: every signal of a protein gets the
#     same values.
#   - 04 brain coloc is cis + cognate only and runs on the protein's cis locus, so it
#     is per protein for cis signals and NA (not tested) for trans signals.
#   - NA also marks a protein that was never tested: no CSF cis instrument and no coloc
#     attempt for S2, no brain eQTL coloc run on its own gene for 04.
#
# Usage: Rscript code/09_statistical_analysis/09d_novel_evidence.R

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/common.R")
source("code/09_statistical_analysis/lib/eqtl_coloc_scope.R")  # COLOC_PP4

OUT_DIR <- "output"
A    <- file.path(OUT_DIR, "09a_novel_candidates.csv")
S2   <- file.path(OUT_DIR, "S2_protein_disease.csv.gz")
PAIR <- "output/05a_disease_pair_evidence.csv.gz"  # 05a cache
TIS  <- file.path(OUT_DIR, "04_protein_tissue_pair_evidence.csv.gz")
DRUG <- "info/bf_protein/drug_targets.tsv.gz"
OUT  <- file.path(OUT_DIR, "09d_novel_evidence.csv")

for (f in c(A, S2, PAIR, TIS, DRUG)) if (!file.exists(f)) stop("missing input: ", f)

# Already reported for the same protein at the same rsID (GWAS Catalog lookup,
# docs/issues.md); not to be called novel. SEMA3G (rs5848) and PVR (rs714948) left the
# list on 2026-10-30: Western 2026 CSF Olink, now a 07b reference, reports both.
CATALOG_KNOWN <- data.table(gene = "TREM2", snp = "rs7232")

# APOE_REGION (lib/common.R); positions in 09a are hg19.

# Join a per-key summary onto `x` and fail if the join changes the row count.
join_n <- function(x, y, by) {
  if (anyDuplicated(y, by = by)) stop("summary not unique on ", paste(by, collapse = "+"))
  n0 <- nrow(x)
  x <- merge(x, y, by = by, all.x = TRUE, sort = FALSE)
  if (nrow(x) != n0) stop("join on ", paste(by, collapse = "+"), " changed rows ", n0, " -> ", nrow(x))
  x
}
collapse <- function(v) if (length(v)) paste(sort(unique(v)), collapse = "; ") else ""

# ---- candidates ----
cand <- fread(A)[is_candidate == TRUE]
if (anyDuplicated(cand$cs_key)) stop("09a candidates not unique on cs_key")
message("candidates: ", nrow(cand), " signals on ", uniqueN(cand$protein_assay), " proteins")

hit <- merge(CATALOG_KNOWN, cand[, .(gene, snp)], by = c("gene", "snp"))
if (nrow(hit) != nrow(CATALOG_KNOWN)) stop("catalog-known signals not all among the candidates")
cand[, catalog_known := paste(gene, snp) %in% CATALOG_KNOWN[, paste(gene, snp)]]
cand[, apoe_region := in_apoe_region(chr, pos)]

ev <- cand[, .(cs_key, protein_assay, gene, ensembl_gene_id, platform, is_cis,
               assoc_class, snp, chr, pos, ea, nea, beta, se, p, catalog_known, apoe_region)]

# ---- S2: FinnGen MR / coloc, CSF, significant pairs only ----
# S2 holds only protein x disease pairs with MR or coloc evidence. A protein absent
# from it was linked to nothing if 05a tested it (MR run or coloc attempted), and was
# not tested (NA) otherwise.
s2 <- fread(S2)[fluid == "CSF"]
s2_lab <- function(disease, mr, co) paste0(disease, " (",
  fifelse(mr & co, "MR+coloc", fifelse(mr, "MR", "coloc")), ")")
s2_sum <- s2[, .(n_disease_mr    = uniqueN(disease[mr_sig %in% TRUE]),
                 n_disease_coloc = uniqueN(disease[coloc_sig %in% TRUE]),
                 n_disease_joint = uniqueN(disease[joint_evidence %in% TRUE]),
                 neuro_diseases  = collapse(s2_lab(disease, mr_sig %in% TRUE, coloc_sig %in% TRUE)[
                   disease_category %in% c("Neurological", "Psychiatric")])),
             by = protein_assay]
ev <- join_n(ev, s2_sum, "protein_assay")
pair <- fread(PAIR, select = c("protein_assay", "fluid", "mr_ok", "coloc_tested"))[fluid == "CSF"]
s2_tested <- pair[, any(mr_ok | coloc_tested), by = protein_assay][V1 == TRUE, protein_assay]
if (length(setdiff(s2$protein_assay, s2_tested))) stop("S2 proteins missing from the 05a tested set")
untested <- !ev$protein_assay %chin% s2_tested
for (v in c("n_disease_mr", "n_disease_coloc", "n_disease_joint")) set(ev, which(is.na(ev[[v]]) & !untested), v, 0L)
ev[is.na(neuro_diseases) & !untested, neuro_diseases := ""]

# ---- 04: brain eQTL coloc of the protein's own gene, cis only ----
# Tested and hit as in 09e (lib/common.R brain_coloc_sets); untested stays NA, not FALSE.
bc <- brain_coloc_sets(TIS)
tis_all <- bc$tis
brain_tested <- bc$tested
tis <- tis_all[coloc_retained == TRUE]
if (tis[is_cis != TRUE | cognate != TRUE, .N]) stop("04 retained rows outside the cis+cognate scope")
tis_sum <- tis[, .(brain_pp4_max      = max(PP4),
                   brain_coloc        = any(PP4 > COLOC_PP4),
                   brain_coloc_tissue = collapse(tissue[PP4 > COLOC_PP4])),
               by = protein_assay]
ev <- join_n(ev, tis_sum, "protein_assay")
ev[is_cis == TRUE & is.na(brain_coloc) & protein_assay %chin% brain_tested,
   `:=`(brain_coloc = FALSE, brain_coloc_tissue = "")]
ev[is_cis == FALSE, c("brain_pp4_max", "brain_coloc", "brain_coloc_tissue") :=
     .(NA_real_, NA, NA_character_)]

# ---- drug targets (Open Targets), per gene ----
drug <- fread(DRUG)
drug_sum <- drug[, .(n_drugs = uniqueN(drug_id),
                     drug_max_stage = if (all(is.na(max_stage))) NA_real_ else as.numeric(max(max_stage, na.rm = TRUE))),
                 by = ensembl_gene_id]
stage_lab <- unique(drug[!is.na(max_stage), .(drug_max_stage = as.numeric(max_stage), drug_stage = stage_label)])
drug_sum <- join_n(drug_sum, stage_lab, "drug_max_stage")
ev <- join_n(ev, drug_sum[, .(ensembl_gene_id, n_drugs, drug_stage)], "ensembl_gene_id")
ev[is.na(n_drugs), `:=`(n_drugs = 0L, drug_stage = "")]

setorder(ev, gene, chr, pos)
fwrite(ev, OUT)

# ---- summary for the log ----
pr <- unique(ev, by = "protein_assay")
message("\nper protein (", nrow(pr), "):")
message("  S2 not tested (no CSF cis instrument or coloc): ", pr[is.na(n_disease_mr), .N])
message("  S2 MR or coloc, any disease:    ", pr[n_disease_mr + n_disease_coloc > 0, .N],
        " (neuro/psych: ", pr[neuro_diseases != "", .N], ")")
pc <- unique(ev[is_cis == TRUE], by = "protein_assay")
message("  brain cognate coloc (cis):      ", pc[brain_coloc %in% TRUE, .N], " of ", pc[!is.na(brain_coloc), .N],
        " tested; ", pc[is.na(brain_coloc), .N], " of ", nrow(pc), " cis proteins not tested")
message("  >= 1 drug / approved drug:      ", pr[n_drugs > 0, .N], " / ",
        pr[drug_stage == "Approved", .N])
message("per signal (", nrow(ev), "):")
message("  catalog-known: ", ev[catalog_known == TRUE, .N], "; APOE region: ", ev[apoe_region == TRUE, .N])
message("wrote ", OUT)
