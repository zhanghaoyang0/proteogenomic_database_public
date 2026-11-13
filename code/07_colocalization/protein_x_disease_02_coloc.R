#!/usr/bin/env Rscript
# Per-chunk coloc runner: BF protein cis-pQTL x disease GWAS, coloc.susie only.
# Exposure loci are the protein's CIS credible sets; the body is in
# coloc_x_disease.R, shared with mb_x_disease_02_coloc.R.

SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_coloc.R"))
source(file.path(SCRIPT_DIR, "coloc_x_disease.R"))

ereg_cache <- list()   # protein cis loci (clump_lead_snp, chr, pos)

run_coloc_x_disease(list(
  coloc_dir = file.path(COLOC_ROOT, "protein_x_disease"),
  pair_ok   = function(pair) pair[!is.na(exposure_id) & !is.na(outcome_id)],
  exposure_loci = function(p) {
    id <- p$exposure_id
    if (is.null(ereg_cache[[id]]))
      ereg_cache[[id]] <<- finemap_loci(finemap_path(p$exposure_source, id), cis_only = TRUE)
    ereg_cache[[id]]
  },
  exposure_rds = function(p, lead) protein_susie_path(p$exposure_source, p$exposure_id, lead),
  id_cols      = function(p) list(exposure_id = p$exposure_id),
  progress_every  = 50L
))
