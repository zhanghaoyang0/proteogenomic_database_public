# Shared loading + scoping for the protein x eQTL coloc analyses.
#
# 04_protein_tissue.R (tissue-of-origin) and 04b_protein_celltype.R
# (cell-type-of-origin) must apply exactly the same assay set and the same
# cis+cognate filter, or their counts are not comparable and any sentence that puts
# them side by side is wrong. That is the reason this is shared code rather than
# copied: it is a correctness constraint, not tidiness.
#
# Paths are relative to the repo root, so callers must run from there.
suppressPackageStartupMessages(library(data.table))

# COLOC_PP4 and coloc_reported() come from lib/common.R, which callers source first.
# A coloc row's lead_snp is the protein locus it ran on (the clump lead); the locus
# is cis when any of its credible sets is (finemap/05), the same rule coloc used to
# pick it. The clump lead's own position is not the test.
SIGNAL_FILE  <- "info/bf_protein/signals.tsv.gz"
COGNATE_FILE <- "info/bf_protein/ensg/assay_ensg.csv"
ANNOT_FILE   <- "info/bf_protein/assay_annot.csv"
KEPT_ASSAYS  <- "output/00_kept_assays.csv"

# SingleBrain cell-type codes -> readable labels. Shared by 04 (its cell-type
# block) and 04b (the cell-type-of-origin figure) so both spell them the same way.
CELLTYPE_LAB <- c(
  Ast = "Astrocyte", End = "Endothelial", Ext = "Excitatory neuron",
  IN = "Inhibitory neuron", MG = "Microglia", OD = "Oligodendrocyte",
  OPC = "OPC"
)

COLOC_RESULT <- function(chain)   # COLOC_ROOT from lib/common.R, which callers source first
  file.path(COLOC_ROOT, sprintf("protein_x_%s/result/coloc_protein_x_%s.tsv", chain, chain))

# BF2 assay annotation: fluid, platform, SomaLogic dropped (plasma-only and
# trans-inflated, so it would break the CSF-vs-plasma comparison), restricted to the
# one highest-N assay per protein and fluid that the rest of the report uses.
load_assay_annot <- function() {
  annot <- fread(ANNOT_FILE,
                 select = c("cohort", "protein_assay", "gene", "tag"))[cohort == "bf2"]
  annot <- unique(annot, by = "protein_assay")
  annot[, exposure_id := paste0("bf2__", protein_assay)]
  annot[, fluid := assay_fluid(tag)]
  annot[, platform := assay_platform(tag)]
  annot <- annot[platform != "SomaLogic"]
  annot <- annot[!is.na(fluid)]
  if (!file.exists(KEPT_ASSAYS))
    stop("Missing ", KEPT_ASSAYS, " -- run the upstream report scripts first.")
  annot[protein_assay %in% fread(KEPT_ASSAYS)$protein_assay]
}

read_coloc <- function(f, res) {
  if (!file.exists(f)) stop("Missing coloc result: ", f)
  d <- fread(f, select = c(
    "exposure_id", "status", "PP4", "tissue", "ENSG", "GENE",
    "cs_id", "lead_snp", "o_lead", "chr", "pos", "coloc_method", "primary"
  ))
  d[, resource := res]
  d
}

# The manuscript estimand: a CIS pQTL lead paired with the protein's COGNATE gene.
# MetaBrain was run broadly for reuse, so it needs the filter; the gtex and sb chains
# are cognate-scoped already and pass it unchanged. Status-only rows survive because
# they record an attempted test, which is what the denominators are built from.
apply_cis_cognate_scope <- function(cb, annot, verbose = TRUE) {
  cb <- merge(cb, annot, by = "exposure_id", all = FALSE)

  lead_scope <- fread(SIGNAL_FILE,
    select = c("trait", "cohort", "signal_id", "is_cis"))[cohort == "bf2"]
  lead_scope <- lead_scope[, .(exposure_id = trait, lead_snp = signal_id, is_cis)]

  cognate <- fread(COGNATE_FILE,
    select = c("cohort", "protein_assay", "ensembl_gene_id"))[cohort == "bf2"]
  cognate[, `:=`(exposure_id = paste0(cohort, "__", protein_assay),
                 ensg_key = sub("\\..*$", "", ensembl_gene_id))]
  cognate <- unique(cognate[, .(exposure_id, ensg_key)])

  cb[, ensg_key := sub("\\..*$", "", ENSG)]
  cb <- merge(cb, lead_scope, by = c("exposure_id", "lead_snp"), all.x = TRUE)
  cb <- merge(cb, cognate[, .(exposure_id, ensg_key, cognate = TRUE)],
              by = c("exposure_id", "ensg_key"), all.x = TRUE)
  cb[is.na(cognate), cognate := FALSE]

  if (verbose) {
    n_out <- cb[status == "ok" & (is.na(is_cis) | is_cis != TRUE | cognate != TRUE), .N]
    if (n_out > 0L)
      message("Dropping retained eQTL coloc rows outside cis+cognate scope: ", n_out)
  }
  cb <- cb[status != "ok" | (is_cis == TRUE & cognate == TRUE)]
  # Reported rows (coloc_reported, lib/common.R); other ok rows count as tested, never
  # as colocalizing.
  cb[, coloc_retained := coloc_reported(status, primary, PP4)]
  cb[]
}
