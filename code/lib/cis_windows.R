# cis windows of protein assays, one row per window (gene +/- 500 kb, one per subunit of a
# complex, overlapping ones merged), written by code/02_protein_qc/03_annot_ensg.R.
# A variant is cis for an assay if it lies in any of that assay's windows. Use these
# instead of the single cis_start/cis_end of assay_ensg.csv, which cover the named gene only.
# Needs data.table loaded by the caller.

# Half-width of the window around each gene: 03_annot_ensg.R builds the windows with it,
# code/lib/mhc.R recovers the gene from a window with it.
CIS_WINDOW_BP <- 5e5

read_cis_windows <- function(repo_root) {
  f <- file.path(repo_root, "info/bf_protein/ensg/assay_cis_windows.csv")
  if (!file.exists(f)) stop("cis windows not found: ", f, " (run code/02_protein_qc/03_annot_ensg.R)")
  w <- fread(f, colClasses = c(chr = "character"))
  # windows are a property of the assay, so bf1 and bf2 rows of one assay must agree
  w <- unique(w[, .(protein_assay, chr, cis_start, cis_end, genes)])
  if (anyDuplicated(w[, .(protein_assay, chr, cis_start)]))
    stop("assay_cis_windows.csv: one assay has different windows in different cohorts")
  w
}

# TRUE where (assay, chr, bp) lies in one of the assay's windows; FALSE otherwise,
# including assays without coordinates.
in_cis_window <- function(assay, chr, bp, windows) {
  q <- data.table(i = seq_along(assay), protein_assay = as.character(assay),
                  chr = sub("^chr", "", as.character(chr)), bp = as.numeric(bp))
  hit <- q[windows, on = .(protein_assay, chr, bp >= cis_start, bp <= cis_end), nomatch = NULL, x.i]
  seq_along(assay) %in% hit
}
