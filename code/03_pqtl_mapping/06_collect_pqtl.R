#!/usr/bin/env Rscript
# RUN LOCATION: CLUSTER — reads clump output under ~/db/shared_data/gwas/{bf1,bf2}/clump/.
# Depends on: info/bf_protein/ensg/assay_ensg.csv (built locally, rsync'd to cluster).
# Output: info/bf_protein/clump_lead_pqtl.csv.gz (gzipped; ~470MB raw -> ~60MB)
#
# Collect genome-wide clump lead SNPs for bf2 protein; label cis/trans.
# is_cis = lead inside any of the assay's windows (info/bf_protein/ensg/assay_cis_windows.csv,
# one per subunit). gene_*/cis_start/cis_end are the named gene's, kept for display only.

library(data.table)

source(path.expand("code/lib/common.R"))   # REPO_ROOT, GWAS_BASE, CLUMP_SUBDIR
CLUMP_DIR <- file.path("clump/protein", CLUMP_SUBDIR)
INFO_FILE <- file.path(REPO_ROOT, "info/bf_protein/ensg/assay_ensg.csv")
OUT_FILE  <- file.path(REPO_ROOT, "info/bf_protein/clump_lead_pqtl.csv.gz")
stopifnot(grepl("\\.csv\\.gz$", OUT_FILE))

ensg <- fread(INFO_FILE)[, .(cohort, protein_assay, gene, ensembl_gene_id,
                              gene_chr, gene_start, gene_end, cis_start, cis_end)]

INPUT_DIR <- "clump/protein/input_5e-8Pthres"
NO_SIG_RE <- "^No significant SNPs in input|No significant --clump results"

# BF1 is a replication lookup and is never clumped for the resource (docs/project_status.md);
# a bf1/clump/protein/ dir from a side analysis (lina_rs5848) is not read.
if (dir.exists(file.path(GWAS_BASE, "bf1", CLUMP_DIR)))
  message("bf1: clump dir exists but is skipped (bf2 only)")
rows <- list()
for (cohort in "bf2") {
  clump_path <- file.path(GWAS_BASE, cohort, CLUMP_DIR)
  if (!dir.exists(clump_path)) { message(cohort, ": no clump dir, skipped (", clump_path, ")"); next }
  # Every clump input must have a .clumped or a "nothing to clump" log.
  stems <- sub("\\.tsv$", "", list.files(file.path(GWAS_BASE, cohort, INPUT_DIR), pattern = "\\.tsv$"))
  done <- file.exists(file.path(clump_path, paste0(stems, ".clumped"))) |
    vapply(file.path(clump_path, paste0(stems, ".log")), function(l)
      file.exists(l) && any(grepl(NO_SIG_RE, readLines(l, warn = FALSE))), logical(1))
  if (!all(done)) stop(cohort, ": ", sum(!done), " of ", length(stems), " clump inputs have no result, e.g. ",
                       paste(head(stems[!done], 3), collapse = ", "))
  files <- list.files(clump_path, pattern = "\\.clumped$", full.names = TRUE)
  message(cohort, ": ", length(files), " clumped files, ", length(stems), " inputs")
  for (f in files) {
    protein_assay <- sub(paste0("^", cohort, "__"), "", sub("\\.clumped$", "", basename(f)))
    dt <- setDT(read.table(f, header = TRUE))
    if (nrow(dt) == 0) stop("empty .clumped: ", f)
    dt[, `:=`(cohort = cohort, protein_assay = protein_assay)]
    rows[[length(rows) + 1]] <- dt[, .(cohort, protein_assay, SNP, CHR, BP, P)]
  }
}

pqtl <- rbindlist(rows)
message("Total lead SNPs: ", nrow(pqtl))

# Join gene coords
stopifnot(!anyDuplicated(ensg, by = c("cohort", "protein_assay")))
n_before <- nrow(pqtl)
pqtl <- merge(pqtl, ensg, by = c("cohort", "protein_assay"), all.x = TRUE)
stopifnot(nrow(pqtl) == n_before)
n_no_match <- sum(is.na(pqtl$gene_chr))
if (n_no_match > 0) message("WARNING: ", n_no_match, " lead SNPs with no gene coords")

# Label cis: lead inside any of the assay's windows
source(file.path(REPO_ROOT, "code/lib/cis_windows.R"))
pqtl[, is_cis := in_cis_window(protein_assay, CHR, BP, read_cis_windows(REPO_ROOT))]

# MHC loci are counted but never fine-mapped or used for MR (../lib/mhc.R).
source(path.expand("code/lib/mhc.R"))
pqtl[, in_mhc := mhc_variant(CHR, BP)]

setorder(pqtl, cohort, protein_assay, P)
message("Cis: ", sum(pqtl$is_cis), "  Trans: ", sum(!pqtl$is_cis), "  MHC: ", sum(pqtl$in_mhc))

fwrite(pqtl, OUT_FILE, compress = "gzip")
message("Written: ", OUT_FILE)
