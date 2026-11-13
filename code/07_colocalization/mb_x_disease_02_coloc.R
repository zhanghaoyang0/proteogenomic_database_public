#!/usr/bin/env Rscript
# Per-chunk coloc runner: MetaBrain brain eQTL x disease GWAS, coloc.susie only —
# the replacement for the legacy coloc.abf mb_x_disease. Exposure loci are the
# gene's credible sets, ALL of them (cis and trans); the body is in
# coloc_x_disease.R, shared with protein_x_disease_02_coloc.R.

SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_coloc.R"))
source(file.path(SCRIPT_DIR, "coloc_x_disease.R"))

ereg_cache <- list()   # MetaBrain loci, cached per tissue_ENSG

run_coloc_x_disease(list(
  coloc_dir = file.path(COLOC_ROOT, "mb_x_disease"),
  pair_ok   = function(pair) {
    pair <- pair[!is.na(tissue) & !is.na(ENSG) & !is.na(outcome_id)]
    pair[, exposure_id := paste0(tissue, "_", ENSG)]
    pair
  },
  exposure_loci = function(p) {
    key <- paste0(p$tissue, "_", p$ENSG)
    if (is.null(ereg_cache[[key]]))
      ereg_cache[[key]] <<- finemap_loci(metabrain_finemap_path(p$tissue, p$ENSG), cis_only = FALSE)
    ereg_cache[[key]]
  },
  exposure_rds = function(p, lead) metabrain_susie_path(p$tissue, p$ENSG, lead),
  id_cols      = function(p) list(exposure_id = p$exposure_id, tissue = p$tissue,
                                  ENSG = p$ENSG, GENE = p$GENE),
  progress_every  = 200L
))
