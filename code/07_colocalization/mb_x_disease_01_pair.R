#!/usr/bin/env Rscript

# Generate coloc pairs (MetaBrain eQTL x disease) — the coloc.susie replacement for
# the legacy coloc.abf mb_x_disease. Same regional pre-filter style as
# protein_x_disease_01_pair.R, with the exposure side = MetaBrain (tissue, ENSG)
# cis windows instead of protein cis windows.
#   Keep (mb gene, disease) pairs where a disease clump lead SNP (p<5e-8) falls in
#   the MetaBrain gene cis window (gene +/- 1 Mb, from the manifest start/end).
# FinnGen and Other outcomes, as in protein_x_disease (Other has been fine-mapped since
# 2026-10-06; it was left out here until 2026-11-01 on the older reason). Output:
#   info/mb_x_disease/coloc_pair.csv  and  coloc_pair_todo.csv

SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_db.R"))
source(file.path(SCRIPT_DIR, "functions_coloc.R"))   # eqtl_finemap_path, disease_finemap_path
source(path.expand("code/lib/gwas_functions.R"))   # P_EGENE_DEFAULT

MB_INFO   <- path.expand("./info/metabrain/trait_lowest_p.csv")
COLOC_DIR <- file.path(COLOC_ROOT, "mb_x_disease")
INFO_DIR  <- file.path(COLOC_DIR, "info")
TMP_DIR   <- file.path(COLOC_DIR, "tmp_chunks")
DB_FILE   <- file.path(COLOC_DIR, "coloc.db")
DB_TABLE  <- "coloc_mb_x_disease"
KEY_COLS  <- c("tissue", "ENSG", "outcome_id")   # the unit of work for this chain
CLUMP_DIRS <- list(
  "FinnGen meta-analysis" = file.path(GWAS_BASE, "finngen/clump", CLUMP_SUBDIR),
  "Other"                 = file.path(GWAS_BASE, "other/clump", CLUMP_SUBDIR)
)
assert_clump_dirs(CLUMP_DIRS)
dir.create(INFO_DIR, showWarnings = FALSE, recursive = TRUE)

t0 <- proc.time()[[3L]]

# ---- 1. MetaBrain (tissue, ENSG) cis windows ----
mb <- fread(MB_INFO)
setnames(mb, c("REGION", "ENSG_id"), c("tissue", "ENSG"))
mb <- mb[!is.na(start) & !is.na(end) & !is.na(CHR)]
mb[, CHR := as.integer(CHR)]
# The same (tissue, gene) set the cognate fine-mapping manifest uses
# (finemap/prep_eqtl_cis_loci.R): eGenes at P < P_EGENE_DEFAULT, one row per tissue +
# gene, keeping the row with the smallest lowest_p (the true cis window; liftover
# artifact rows on other chromosomes can be wider). Without the eGene filter this
# picked up 2,349 tissue-genes whose only SuSiE fits were from the June overlap-scope
# run (no per-SNP N filter, max_iter 100).
mb <- mb[lowest_p < P_EGENE_DEFAULT]
setorder(mb, tissue, GENE, lowest_p)
mb <- unique(mb, by = c("tissue", "GENE"))
mb_win <- unique(mb[, .(tissue, ENSG, GENE, e_chr = CHR, cis_start = start, cis_end = end)])
msg("MetaBrain cis windows: ", nrow(mb_win), " (tissue,ENSG)")
# Genes that encode a BF2 protein only: the protein -> own-gene brain expression ->
# disease triangle. The scope used to be "whatever MetaBrain fit is on disk", i.e.
# every gene near any protein CS (cis or trans) from the old OVERLAP fine-mapping.
# Cognate genes of the proteins protein_x_mb tests, exactly as the fine-mapping
# manifest (finemap/prep_eqtl_cis_loci.R) defines them, so every pair here has a fit
# made by the current code. All BF2 assays' genes used to be taken, and 2,094 of
# those tissue-genes only had June fits.
PROT_LIST <- file.path(COLOC_ROOT, "protein_x_mb/info/protein_list.csv")
if (!file.exists(PROT_LIST)) stop("missing ", PROT_LIST, " -- run protein_x_eqtl_01_pair.R mb first")
prot <- fread(PROT_LIST)
ae <- fread(path.expand("./info/bf_protein/ensg/assay_ensg.csv"))
ae[, exposure_id := paste0(cohort, "__", protein_assay)]
cog <- unique(sub("\\..*$", "", trimws(unlist(strsplit(ae[exposure_id %in% prot$exposure_id & !is.na(ensembl_gene_id), ensembl_gene_id], ",")))))
mb_win <- mb_win[sub("\\..*$", "", ENSG) %in% cog]
msg("... restricted to genes encoding a BF2 protein: ", nrow(mb_win), " (tissue,ENSG)")

# ---- 2. Disease clump lead SNPs ----
disease_snps <- read_disease_leads(CLUMP_DIRS)
msg("Disease clump lead SNPs: ", nrow(disease_snps), " across ", uniqueN(disease_snps$outcome_id), " outcomes")

# ---- 3. Regional overlap: disease SNP in MetaBrain gene cis window ----
hits <- merge(mb_win, disease_snps, by.x = "e_chr", by.y = "CHR", allow.cartesian = TRUE)
hits <- hits[BP >= cis_start & BP <= cis_end]
# o_leads: the disease clump leads that put the pair here (see protein_x_disease_01_pair.R)
stopifnot(!anyNA(hits$SNP))
pair <- hits[, .(o_leads = paste(sort(unique(SNP)), collapse = ";")),
             by = .(tissue, ENSG, GENE, outcome_id, outcome_source)]
pair[, exposure_source := "metabrain"]
msg("Pairs passing regional filter: ", nrow(pair))

# Keep only mb (tissue,ENSG) whose current fine-mapping table has at least one
# credible-set row.  Old SuSiE RDS files can survive a rerun that produces no CS,
# so their mere presence is not evidence that the current locus is usable.
MB_FINEMAP_DIR <- file.path(GWAS_BASE, "metabrain/finemapping")
mb_fm_files <- list.files(MB_FINEMAP_DIR, "\\.finemapped\\.tsv\\.gz$", full.names = TRUE)
# A header-only table (no CS) has no row; an unreadable one stops the run rather than
# silently dropping the gene from every pair.
has_finemap_row <- vapply(mb_fm_files, function(f) {
  nrow(tryCatch(fread(f, nrows = 1L, showProgress = FALSE),
                error = function(e) stop("cannot read ", f, ": ", conditionMessage(e)))) > 0L
}, logical(1L))
mb_te_with_fit <- sub("\\.finemapped\\.tsv\\.gz$", "", basename(mb_fm_files[has_finemap_row]))
pair[, te := paste0(tissue, "_", ENSG)]
n_pre <- nrow(pair)
pair <- pair[te %in% mb_te_with_fit]
pair[, te := NULL]
msg("kept ", nrow(pair), "/", n_pre, " pairs whose current mb fine-map has rows")

msg("pairs by outcome source: ", paste(pair[, .N, by = outcome_source][, paste(outcome_source, N)], collapse = ", "))

setorder(pair, outcome_id, tissue, ENSG)
out_pair <- file.path(INFO_DIR, "coloc_pair.csv")
fwrite(pair, out_pair)
msg("Final coloc pairs: ", nrow(pair), " -> ", out_pair)

# ---- todo = pairs not yet in the ledger, or stale ----
pair[, input_mtime := pmax(as.numeric(file.mtime(eqtl_finemap_path("metabrain", tissue, ENSG))),
                           as.numeric(file.mtime(disease_finemap_path(outcome_source, outcome_id))),
                           na.rm = TRUE)]
write_coloc_todo(pair, DB_FILE, DB_TABLE, KEY_COLS, TMP_DIR, INFO_DIR)
