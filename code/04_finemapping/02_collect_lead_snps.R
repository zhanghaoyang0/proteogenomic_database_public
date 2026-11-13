#!/usr/bin/env Rscript
# RUN LOCATION: CLUSTER — reads finemap output under ~/db/shared_data/gwas/bf2/finemapping/<subdir> (BF1 is not fine-mapped).
# Depends on: info/bf_protein/ensg/assay_ensg.csv (cohort, protein_assay -> ensg, gene).
# Output: info/bf_protein/finemap_lead_snps.csv — one row per credible set
#         using PIP-max SNP as the "lead" (robust to truncated lead_snp column upstream).
# Usage: Rscript 02_collect_lead_snps.R [finemap_subdir]
#   finemap_subdir defaults to "protein"; pass "protein" explicitly or omit.

suppressPackageStartupMessages(library(data.table))

args          <- commandArgs(trailingOnly = TRUE)
FINEMAP_SUBDIR <- if (length(args) >= 1) args[1] else "protein"

source(path.expand("code/lib/common.R"))   # GWAS_BASE, msg
INFO_DIR  <- path.expand("./info")
ENSG      <- file.path(INFO_DIR, "bf_protein", "ensg", "assay_ensg.csv")
OUT       <- file.path(INFO_DIR, "bf_protein", "finemap_lead_snps.csv")

t0  <- proc.time()[[3L]]

ensg <- fread(ENSG, select = c("cohort", "protein_assay", "gene", "ensembl_gene_id"))
setnames(ensg, "ensembl_gene_id", "ensg")
msg("ensg map: ", nrow(ensg), " rows")

collect_cohort <- function(coh) {
  fm_dir <- file.path(GWAS_BASE, coh, "finemapping", FINEMAP_SUBDIR)
  files  <- list.files(fm_dir, pattern = "\\.finemapped\\.tsv\\.gz$", full.names = TRUE)
  msg(coh, ": ", length(files), " finemap files")
  if (!length(files)) return(NULL)
  fm <- rbindlist(lapply(files, fread), fill = TRUE)
  if (!nrow(fm)) return(NULL)
  # PIP-max SNP per (trait, cs_id) as the lead; also keep cs_size.
  setorder(fm, trait, cs_id, -pip)
  cs_size <- fm[, .(cs_size = .N), by = .(trait, cs_id)]
  lead    <- unique(fm, by = c("trait", "cs_id"))
  lead    <- merge(lead, cs_size, by = c("trait", "cs_id"))
  prefix  <- paste0("^", coh, "__")
  lead[, cohort := coh]
  lead[, protein_assay := sub(prefix, "", trait)]
  lead[, .(cohort, protein_assay, trait, cs_id, cs_size, cs_purity,
           chr, pos, snp, ea, nea, beta, se, p, max_pip = pip, is_cis)]
}

out <- rbindlist(lapply("bf2", collect_cohort), fill = TRUE)
msg("combined leads: ", nrow(out), " rows")

# is_cis comes from 01_finemapping.R (top-PIP variant in the assay's cis window), the
# one definition; it is not recomputed here.
if (anyNA(out$is_cis)) stop(sum(is.na(out$is_cis)), " credible sets have no is_cis")
n0  <- nrow(out)
out <- merge(out, ensg, by = c("cohort", "protein_assay"), all.x = TRUE)
if (nrow(out) != n0) stop("assay_ensg.csv is not unique per (cohort, protein_assay)")
out[, hc_p5e8 := !is.na(p) & p < P_GWS]

setcolorder(out, c("cohort", "protein_assay", "ensg", "gene", "trait",
                   "chr", "pos", "snp", "ea", "nea", "beta", "se", "p",
                   "cs_id", "cs_size", "cs_purity", "max_pip",
                   "is_cis", "hc_p5e8"))
setorder(out, cohort, protein_assay, cs_id)

# Round numeric columns to 3 significant figures.
for (col in c("beta", "se", "p", "cs_purity", "max_pip")) {
  if (is.numeric(out[[col]])) out[, (col) := signif(get(col), 3)]
}

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
fwrite(out, OUT)
msg("wrote -> ", OUT)

# ---- breakdown ----
msg("")
msg("--- breakdown ---")
cat("rows by cohort:\n");           print(out[, .N, by = cohort])
cat("\nrows by is_cis:\n");          print(out[, .N, by = is_cis])
cat("\nrows by hc_p5e8:\n");         print(out[, .N, by = hc_p5e8])
cat("\ncis × hc_p5e8 cross:\n");     print(out[, .N, by = .(is_cis, hc_p5e8)][order(-N)])
cat("\nNA ensg (unmapped):\n");      print(out[is.na(ensg), .N, by = cohort])
cat("\ncs_size quantiles:\n");       print(quantile(out$cs_size, c(0,0.25,0.5,0.75,0.95,1)))
cat("\nmax_pip quantiles:\n");       print(round(quantile(out$max_pip, c(0,0.25,0.5,0.75,0.95,1), na.rm=TRUE), 3))

msg("done in ", round(proc.time()[[3L]] - t0, 1), "s")
