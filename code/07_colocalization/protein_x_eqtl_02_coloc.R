#!/usr/bin/env Rscript
# Per-chunk coloc runner for the three COGNATE-scope eQTL chains, coloc.susie only:
#   gtex  GTEx multi-tissue eQTL
#   sb    SingleBrain single-nucleus eQTL ("tissue" is a brain cell type)
#   mb    MetaBrain bulk brain eQTL (had its own copy of this file until 2026-10-06)
# A protein's CIS credible sets are colocalized against the eQTL of ITS OWN gene
# only (cognate ENSG from assay_ensg.csv), across all tissues / cell types. The
# body is in coloc_protein_x_qtl.R.
# Usage: called with the chain as argv[1]; see submit_coloc_chain.sh

SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_coloc.R"))
source(file.path(SCRIPT_DIR, "functions_db.R"))
source(file.path(SCRIPT_DIR, "coloc_protein_x_qtl.R"))

CHAIN    <- commandArgs(trailingOnly = TRUE)[1L]
RESOURCE <- c(gtex = "gtex", sb = "singlebrain", mb = "metabrain")[CHAIN]
if (is.na(RESOURCE)) stop("usage: protein_x_eqtl_02_coloc.R <gtex|sb|mb>  (got: ",
                          if (is.na(CHAIN)) "nothing" else CHAIN, ")")

run_coloc_protein_x_qtl(list(
  coloc_dir = file.path(COLOC_ROOT, sprintf("protein_x_%s", CHAIN)),
  qtl_info  = sprintf("./info/%s/trait_lowest_p.csv", RESOURCE),
  qtl_label = RESOURCE,
  qtl_finemap = function(tissue, ensg) eqtl_finemap_path(RESOURCE, tissue, ensg),
  qtl_rds     = function(tissue, ensg, lead) eqtl_susie_path(RESOURCE, tissue, ensg, lead),
  # Status string stays per-chain so it keeps matching the rows already in each
  # ledger: "no_gtex_cs", "no_sb_cs", "no_mb_cs".
  no_cs_status = sprintf("no_%s_cs", CHAIN)
))
