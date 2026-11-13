#!/usr/bin/env Rscript
# Merge the per-cell-type outputs of singlebrain_clean.R into the two repo-tracked
# files, matching how metabrain/gtex are laid out:
#
#   info/singlebrain/trait_lowest_p.csv   <- what prep_*_cis_loci.R reads (hg19)
#   info/singlebrain/egene_audit.csv.gz   <- every tested gene, for absence claims
#
# The per-cell-type directories that feed these stay gitignored (same as
# info/gtex/trait_lowest_p/); only the merged files are committed.
# Run after the clean array finishes: Rscript code/06_external_data/singlebrain_collect_info.R
suppressPackageStartupMessages(library(data.table))

project_root <- "."
info <- file.path(project_root, "info/singlebrain")
cts  <- fread(file.path(info, "datasets.csv"))$cell_type

read_all <- function(sub, pat) {
  f <- list.files(file.path(info, sub), pattern = pat, full.names = TRUE)
  if (length(f) == 0) stop("No files matching ", pat, " under ", file.path(info, sub))
  rbindlist(lapply(f, fread), use.names = TRUE, fill = TRUE)
}

# 1. Trait manifest (kept genes, hg19 positions)
tl <- read_all("trait_lowest_p", "\\.csv$")
missing <- setdiff(cts, unique(tl$REGION))
if (length(missing)) stop("trait_lowest_p missing cell type(s): ", paste(missing, collapse = ", "))
setorder(tl, REGION, GENE, CHR)
fwrite(tl, file.path(info, "trait_lowest_p.csv"))

# 2. Audit of every tested gene (from top_assoc, so it covers dropped genes too)
au <- rbindlist(lapply(cts, function(c) {
  f <- file.path(info, "egene_audit", paste0(c, "_all_tested.csv"))
  if (!file.exists(f)) stop("Missing audit file: ", f)
  d <- fread(f, select = c("feature", "Fixed_P", "Random_P", "qval"))
  d[, `:=`(cell_type = c, GENE = sub("\\..*$", "", feature))]
  d
}), use.names = TRUE)
setcolorder(au, c("cell_type", "GENE", "feature"))
fwrite(au, file.path(info, "egene_audit.csv.gz"))

message(sprintf("trait_lowest_p.csv: %d rows, %d cell types | egene_audit.csv.gz: %d rows",
                nrow(tl), uniqueN(tl$REGION), nrow(au)))
