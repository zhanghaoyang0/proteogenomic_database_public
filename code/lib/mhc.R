# The MHC (hg19), one definition for all of post-GWAS. The numbers live in
# mhc_hg19.range so the bash clump step reads the same ones.
#
# Rule: an MHC signal is discovered and counted (one locus per trait, clump keeps
# the minimum-P lead) but is not fine-mapped, and is not used for MR. Its
# long-range LD cannot be separated with a 1000G reference (HLA-A plasma clumped
# into 19 "independent" leads). Without a credible set it also gets no coloc and
# no BF1 replication.
#
#   source(".../code/lib/mhc.R")
#   mhc_variant(chr, pos)        TRUE when the variant is inside the MHC
#   mhc_gene(chr, start, end)    TRUE when the gene body overlaps the MHC
#   mhc_cis_window(chr, cis_start, cis_end)  TRUE when the gene of a cis window
#                                (gene +/- CIS_WINDOW_BP, assay_cis_windows.csv) overlaps it;
#                                needs code/lib/cis_windows.R sourced first
# Missing coordinates count as outside (FALSE), never NA.

MHC <- local({
  r <- read.table(path.expand("code/lib/mhc_hg19.range"),
                  col.names = c("chr", "start", "end", "label"))
  stopifnot(nrow(r) == 1L)
  list(chr = as.character(r$chr), start = r$start, end = r$end)
})

mhc_variant <- function(chr, pos)
  (as.character(chr) == MHC$chr & pos >= MHC$start & pos <= MHC$end) %in% TRUE

mhc_gene <- function(chr, start, end)
  (as.character(chr) == MHC$chr & end >= MHC$start & start <= MHC$end) %in% TRUE

# Per window, not per assay: a complex with one subunit in the MHC keeps the windows of
# its other subunits (CSNK2A1/CSNK2B). For a one-window assay this is mhc_gene().
mhc_cis_window <- function(chr, cis_start, cis_end, half = CIS_WINDOW_BP)
  mhc_gene(chr, cis_start + half, cis_end - half)
