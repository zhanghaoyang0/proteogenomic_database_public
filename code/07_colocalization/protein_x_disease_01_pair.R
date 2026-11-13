#!/usr/bin/env Rscript

# Generate coloc pairs (protein x disease).
# Strategy: regional pre-filter (Cell-paper style, independent of MR results).
#   Keep (protein, disease) pairs where the disease has a clump lead SNP (p<5e-8)
#   within the protein's cis window (gene +/- 500kb).
# This lets coloc run in parallel with MR rather than depending on MR results.
# Output: coloc_pair.csv (full filtered pair list), coloc_pair_todo.csv (not in tmp_chunks yet).

SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_db.R"))
source(file.path(SCRIPT_DIR, "functions_coloc.R"))   # finemap_path, disease_finemap_path
source(file.path(SCRIPT_DIR, "..", "lib", "mhc.R"))

TRAIT_INFO <- file.path(REPO_ROOT, "info", "trait_info.csv")
COLOC_DIR <- file.path(COLOC_ROOT, "protein_x_disease")
INFO_DIR  <- file.path(COLOC_DIR, "info")
TMP_DIR   <- file.path(COLOC_DIR, "tmp_chunks")
DB_FILE   <- file.path(COLOC_DIR, "coloc.db")
DB_TABLE  <- "coloc_protein_x_disease"
KEY_COLS  <- c("exposure_id", "outcome_id")

PQTL_FILE <- file.path(REPO_ROOT, "info", "bf_protein", "clump_lead_pqtl.csv.gz")
CLUMP_DIRS <- list(
  "FinnGen meta-analysis" = file.path(GWAS_BASE, "finngen/clump", CLUMP_SUBDIR),
  "Other"                 = file.path(GWAS_BASE, "other/clump", CLUMP_SUBDIR)
)
assert_clump_dirs(CLUMP_DIRS)

dir.create(INFO_DIR, showWarnings = FALSE, recursive = TRUE)

t0 <- proc.time()[[3L]]

# ---- 1. Protein cis windows ----
msg("Loading protein cis windows from pqtl file...")
pqtl <- fread(PQTL_FILE, select = c("cohort", "protein_assay", "CHR", "BP", "is_cis"))
# cis windows, one per subunit; windows of a gene in the MHC are dropped (../lib/mhc.R).
source(file.path(REPO_ROOT, "code/lib/cis_windows.R"))
win <- read_cis_windows(REPO_ROOT)
win <- win[!mhc_cis_window(chr, cis_start, cis_end)]
# Exposures: bf2 proteins with a cis lead outside the MHC, in one of those windows --
# the same set MR uses, stated here rather than taken from MR's mr_pair.csv.gz, so
# coloc does not depend on MR having run. A pair can come from any of their windows.
pqtl <- pqtl[cohort == "bf2" & is_cis == TRUE & !mhc_variant(CHR, BP)]
pqtl <- pqtl[in_cis_window(protein_assay, CHR, BP, win)]
cis_win <- unique(win[protein_assay %chin% pqtl$protein_assay,
                      .(exposure_id = paste0("bf2__", protein_assay), p_chr = chr, cis_start, cis_end)])
msg("Cis windows: ", nrow(cis_win), " for ", uniqueN(cis_win$exposure_id), " unique proteins")

# ---- 2. Disease clump lead SNPs ----
msg("Loading disease clump lead SNPs...")
disease_snps <- read_disease_leads(CLUMP_DIRS)
msg("Disease clump lead SNPs: ", nrow(disease_snps), " across ", uniqueN(disease_snps$outcome_id), " outcomes")

# ---- 3. Regional overlap: disease SNP in protein cis window ----
msg("Computing regional overlap (disease SNP in protein cis window)...")
# Join on chromosome; then filter BP within window
cis_win[, p_chr := as.integer(p_chr)]

# Cross-join within each chromosome, then filter
hits <- merge(cis_win, disease_snps, by.x = "p_chr", by.y = "CHR", allow.cartesian = TRUE)
hits <- hits[BP >= cis_start & BP <= cis_end]
# o_leads: the disease clump leads that put the pair here; the runner uses them to
# tell "no disease CS" from "a disease CS, but not near the pQTL lead".
stopifnot(!anyNA(hits$SNP))
pair_keys <- hits[, .(o_leads = paste(sort(unique(SNP)), collapse = ";")),
                  by = .(exposure_id, outcome_id, outcome_source)]
msg("Pairs passing regional filter: ", nrow(pair_keys))

# ---- 4. Attach exposure_source from trait_info.csv ----
ti <- fread(TRAIT_INFO, select = c("trait_id", "source"))   # trait_id is unique there
n_pre <- nrow(pair_keys)
pair <- merge(pair_keys, ti[, .(exposure_id = trait_id, exposure_source = source)], by = "exposure_id")
if (nrow(pair) != n_pre) stop(n_pre - nrow(pair), " exposures are not in ", TRAIT_INFO)
pair <- pair[, .(exposure_id, outcome_id, exposure_source, outcome_source, o_leads)]

# `Other` outcomes are fine-mapped by `submit_01_finemapping_external.sh other`
# (since 2026-10-06, for Bellenguez AD / Jansen CSF); run that before this script.
msg("pairs by outcome source: ",
    paste(pair[, .N, by = outcome_source][, paste0(outcome_source, "=", N)], collapse = ", "))

setorder(pair, outcome_id, exposure_id)
msg("Final coloc pairs: ", nrow(pair))

out_pair <- file.path(INFO_DIR, "coloc_pair.csv")
fwrite(pair, out_pair)
msg("wrote -> ", out_pair)

# ---- todo = pairs not yet in the ledger, or stale ----
pair[, input_mtime := pmax(as.numeric(file.mtime(finemap_path(exposure_source, exposure_id))),
                           as.numeric(file.mtime(disease_finemap_path(outcome_source, outcome_id))),
                           na.rm = TRUE)]
write_coloc_todo(pair, DB_FILE, DB_TABLE, KEY_COLS, TMP_DIR, INFO_DIR)
