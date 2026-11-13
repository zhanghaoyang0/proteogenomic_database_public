#!/usr/bin/env Rscript
# Build info/singlebrain/datasets.csv from the paper's own Supplementary Table 4.
#
# Why this exists: the Zenodo sumstats carry no per-SNP N, and the record's
# description only gives the overall "983 donors". That number is wrong for every
# cell type -- donors contributing <10 nuclei of a type are dropped for that type,
# so the real n ranges 677 (endothelial) to 945 (astrocyte/oligodendrocyte). N goes
# into SuSiE's effective sample size, so the per-cell-type value matters.
#
# Source: Nat Genet 2026, doi:10.1038/s41588-026-02541-x, Supplementary Table 4
# ("Table4_eGenes": Cohort, Cell type, Samples (n), eGenes (n), Test genes (n)).
# The ESM xlsx is openly downloadable, no auth needed.
#
# Usage: Rscript code/06_external_data/singlebrain_sample_size.R
# Re-run to refresh; it overwrites datasets.csv and prints what changed.
suppressPackageStartupMessages({ library(data.table); library(readxl) })

project_root <- "."
out_csv  <- file.path(project_root, "info/singlebrain/datasets.csv")
esm_dir  <- path.expand("~/db/shared_data/gwas/singlebrain/raw")
esm_file <- file.path(esm_dir, "41588_2026_2541_MOESM4_ESM.xlsx")
esm_url  <- paste0("https://static-content.springer.com/esm/",
                   "art%3A10.1038%2Fs41588-026-02541-x/MediaObjects/",
                   "41588_2026_2541_MOESM4_ESM.xlsx")

dir.create(esm_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(esm_file)) {
  cat("Downloading Supplementary Tables ...\n")
  if (download.file(esm_url, esm_file, mode = "wb", quiet = TRUE) != 0)
    stop("Failed to download ", esm_url)
}
if (file.size(esm_file) < 1e6) stop("ESM file looks truncated: ", esm_file)

# Row 1 is a title line, row 2 the real header -> skip = 1.
tab <- as.data.table(read_excel(esm_file, sheet = "Table4_eGenes", skip = 1))
setnames(tab, make.names(names(tab)))
tab <- tab[Cohort == "SingleBrain"]
if (nrow(tab) == 0) stop("No SingleBrain rows in Table4_eGenes -- sheet layout changed?")

# Paper's cell-type names -> the short codes used in the Zenodo file names.
# Matched exactly (after trimming) so a renamed row fails loudly instead of
# silently landing on the wrong cell type -- note "Oligodendrocyte" is a prefix of
# "Oligodendrocyte progenitor cell", which a grepl-based match would confuse.
map <- data.table(
  paper_name = c("Astrocyte", "Endothelial cell", "Excitatory neuron",
                 "Inhibitory neuron", "Microglia", "Oligodendrocyte",
                 "Oligodendrocyte progenitor cell"),
  cell_type  = c("Ast", "End", "Ext", "IN", "MG", "OD", "OPC"),
  label      = c("Astrocytes", "Endothelial cells", "Excitatory neurons",
                 "Inhibitory neurons", "Microglia", "Oligodendrocytes",
                 "Oligodendrocyte progenitor cells"))

tab[, paper_name := trimws(Cell.type)]
d <- merge(map, tab, by = "paper_name", all.x = TRUE)
if (anyNA(d$Samples..n.))
  stop("Unmatched cell type(s): ", paste(d[is.na(Samples..n.), paper_name], collapse = ", "),
       "\nTable4 has: ", paste(sort(unique(tab$paper_name)), collapse = " | "))

out <- d[, .(cell_type, label,
             n_donors    = as.integer(Samples..n.),
             egenes_paper = as.integer(eGenes..n.),
             genes_tested = as.integer(Test.genes..n.))]
setorder(out, cell_type)

if (file.exists(out_csv)) {
  old <- fread(out_csv)
  if ("n_donors" %in% names(old)) {
    cmp <- merge(old[, .(cell_type, n_old = n_donors)], out[, .(cell_type, n_new = n_donors)],
                 by = "cell_type")
    chg <- cmp[n_old != n_new]
    if (nrow(chg)) { cat("n_donors changed:\n"); print(chg) } else cat("n_donors unchanged\n")
  }
}

fwrite(out, out_csv)
print(out)
message("Wrote ", out_csv, " from Supplementary Table 4 (", basename(esm_file), ")")
