#!/usr/bin/env Rscript

# Build the protein list for a protein x molecular-QTL coloc chain.
# Source: all BF proteins with >=1 finemap CS row (independent of MR sig).
# Coloc results are reusable downstream (intersect with MR later).
# Output: protein_list.csv (full), protein_list_todo.csv (not yet in the coloc.db ledger).

SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_db.R"))
source(file.path(SCRIPT_DIR, "functions_coloc.R"))   # GWAS_BASE, finemap_path

# Which molecular-QTL chain to build for. The three differ only in where their
# results live -- the protein selection, the ledger semantics and the todo logic are
# identical -- so they share one script rather than three 89%-identical copies.
#   mb   MetaBrain bulk brain eQTL
#   gtex GTEx multi-tissue eQTL
#   sb   SingleBrain single-nucleus (cell-type) brain eQTL
CHAIN <- commandArgs(trailingOnly = TRUE)[1L]
CHAIN_LABEL <- c(mb = "MetaBrain", gtex = "GTEx eQTL",
                 sb = "SingleBrain single-nucleus eQTL")[CHAIN]
if (is.na(CHAIN_LABEL)) stop("usage: <script> <mb|gtex|sb>  (got: ",
                             if (is.na(CHAIN)) "nothing" else CHAIN, ")")

INFO_DIR    <- path.expand("./info")
COLOC_DIR   <- file.path(COLOC_ROOT, sprintf("protein_x_%s", CHAIN))
OUT_INFO    <- file.path(COLOC_DIR, "info")
DB_FILE     <- file.path(COLOC_DIR, "coloc.db")
DB_TABLE    <- sprintf("coloc_protein_x_%s", CHAIN)
BF1_FINEMAP <- file.path(GWAS_BASE, "bf1/finemapping/protein")
BF2_FINEMAP <- file.path(GWAS_BASE, "bf2/finemapping/protein")
dir.create(OUT_INFO, showWarnings = FALSE, recursive = TRUE)

t0 <- proc.time()[[3L]]

finemap_has_cs <- function(dir) {
  files <- list.files(dir, pattern = "\\.finemapped\\.tsv\\.gz$", full.names = TRUE)
  ids   <- sub("\\.finemapped\\.tsv\\.gz$", "", basename(files))
  keep  <- vapply(files, function(f) length(readLines(gzfile(f), n = 2)) >= 2, logical(1L))
  ids[keep]
}

# BF2 is the discovery/exposure cohort for coloc, matching protein_x_disease and
# the MR arms (docs/project_status.md). BF1 finemapping exists and this script used to
# pick it up -- the selection was never cohort-filtered -- which put 1,060 BF1
# exposures into this table. BF1 CSF pQTL have per-variant N of 525-792, an order
# of magnitude below BF2, so mixing them into one result table forces every
# downstream count to be reported per cohort. Set COLOC_INCLUDE_BF1=1 to restore.
INCLUDE_BF1 <- identical(Sys.getenv("COLOC_INCLUDE_BF1"), "1")
bf1 <- if (INCLUDE_BF1) finemap_has_cs(BF1_FINEMAP) else character(0)
bf2 <- finemap_has_cs(BF2_FINEMAP)
if (!INCLUDE_BF1) msg("BF2-exposure only (COLOC_INCLUDE_BF1=1 to include BF1)")
msg("bf1 finemap with CS: ", length(bf1), " | bf2: ", length(bf2))

ti <- fread(file.path(INFO_DIR, "trait_info.csv"))
is_bf_prot <- ti$source %in% c("BioFINDER1", "BioFINDER2") &
              ti$category %in% c("Protein CSF", "Protein plasma")
prot <- ti[is_bf_prot & trait_id %in% c(bf1, bf2),
           .(exposure_id = trait_id, exposure_source = source)]
# trait_info can carry more than one row per trait_id; without this the same protein
# lands in several chunks and gets colocalised twice.
prot <- unique(prot, by = "exposure_id")
setorder(prot, exposure_id)
msg("chain: protein x ", CHAIN_LABEL, " -> ", COLOC_DIR)
msg("protein list: ", nrow(prot), " unique proteins")

out_full <- file.path(OUT_INFO, "protein_list.csv")
fwrite(prot, out_full)
msg("wrote -> ", out_full)

# Done = present in the ledger. (Legacy data/ files no longer count as done.)
done <- db_done_keys(DB_FILE, DB_TABLE, "exposure_id")
msg("ledger done: ", length(done))
# Staged but not yet collected counts as done too, as un-imported chunks do in
# write_coloc_todo: a resubmit after one failed job would otherwise redo the chain.
staged <- list.files(file.path(COLOC_DIR, "staging"), pattern = "^job_.*\\.tsv$", full.names = TRUE)
staged_rows <- if (length(staged))
  rbindlist(lapply(staged, fread, select = c("exposure_id", "computed_at"))) else NULL
# Stale = in the ledger but computed before the protein's fine-map, or the fine-map of
# one of its cognate genes in any tissue of this resource, changed.
RESOURCE <- c(mb = "metabrain", gtex = "gtex", sb = "singlebrain")[[CHAIN]]
qtl_files <- list.files(file.path(GWAS_BASE, RESOURCE, "finemapping"),
                        pattern = "\\.finemapped\\.tsv\\.gz$", full.names = TRUE)
qtl_m <- data.table(ENSG = sub("^.*_(ENSG[0-9]+).*$", "\\1", basename(qtl_files)),
                    m = as.numeric(file.mtime(qtl_files)))[, .(m_qtl = max(m)), by = ENSG]
ae <- fread(file.path(INFO_DIR, "bf_protein/ensg/assay_ensg.csv"))
ae[, exposure_id := paste0(cohort, "__", protein_assay)]
ae <- ae[exposure_id %in% prot$exposure_id & !is.na(ensembl_gene_id),
         .(ENSG = sub("\\..*$", "", trimws(unlist(strsplit(ensembl_gene_id, ","))))), by = exposure_id]
inputs <- merge(prot, merge(ae, qtl_m, by = "ENSG")[, .(m_qtl = max(m_qtl)), by = exposure_id],
                by = "exposure_id", all.x = TRUE)
inputs[, input_mtime := pmax(as.numeric(file.mtime(finemap_path(exposure_source, exposure_id))),
                             m_qtl, na.rm = TRUE)]
stale <- db_stale_keys(DB_FILE, DB_TABLE, "exposure_id", inputs[, .(exposure_id, input_mtime)])
if (!is.null(staged_rows)) {
  # Staged rows count as done only if not older than their fine-maps (as the ledger).
  s <- merge(staged_rows[, .(computed_at = min(computed_at)), by = exposure_id],
             inputs[, .(exposure_id, input_mtime)], by = "exposure_id", all.x = TRUE)
  ok <- s[is.na(input_mtime) | (!is.na(computed_at) & computed_at >= input_mtime), exposure_id]
  done <- unique(c(setdiff(done, stale), ok)); stale <- character(0)
  msg("plus ", length(staged), " uncollected staging file(s): ", length(ok), " current exposures")
}
msg("stale (a fine-map is newer than its ledger rows): ", length(stale))
todo <- prot[!exposure_id %in% setdiff(done, stale)]
out_todo <- file.path(OUT_INFO, "protein_list_todo.csv")
fwrite(todo, out_todo)
msg("todo: ", nrow(todo), "/", nrow(prot), " -> ", out_todo)
