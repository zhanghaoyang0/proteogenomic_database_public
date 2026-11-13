#!/usr/bin/env Rscript
# Collect coloc protein x disease chunk results into the SQLite ledger, then export the
# flat result tsv that downstream already reads.
#
# The work is in collect_chunked.R, shared with the other chunk-based collector.
SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_db.R"))
source(file.path(SCRIPT_DIR, "collect_chunked.R"))
source(file.path(SCRIPT_DIR, "functions_coloc.R"))   # add_component_flags

t0 <- proc.time()[[3L]]
collect_chunked(
  coloc_dir = file.path(COLOC_ROOT, "protein_x_disease"),
  db_table  = "coloc_protein_x_disease",
  key_cols  = c("exposure_id", "outcome_id"),
  annotate  = function(all, pair) {
    osrc <- unique(pair[, .(outcome_id, outcome_source)])
    if (anyDuplicated(osrc$outcome_id)) stop("an outcome_id with two outcome_source values")
    o_src <- osrc$outcome_source[match(all$outcome_id, osrc$outcome_id)]
    e_src <- protein_src(all$exposure_id)
    add_component_flags(all,
      e_rds = protein_susie_path(e_src, all$exposure_id, lead_or_na(all$lead_snp)),
      e_fm  = finemap_path(e_src, all$exposure_id),
      o_rds = disease_susie_path(o_src, all$outcome_id, lead_or_na(all$o_lead)),
      o_fm  = disease_finemap_path(o_src, all$outcome_id),
      e_cis_only = TRUE)
  }
)
