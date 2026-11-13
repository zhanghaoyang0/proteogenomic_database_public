#!/usr/bin/env Rscript
# Filter disease clump leads (FinnGen, or `other` with COHORT=other) to those inside ANY
# bf2 protein cis window (gene +-500kb, one per subunit) widened by FM_WINDOW_BP. coloc
# accepts a disease credible set whose top-PIP lies within 500 kb of a cis pQTL lead, and
# a cis lead can sit at the edge of its window, so every disease locus coloc could use is
# fine-mapped (r013: PAM x T2D_WIDE existed only because a bogus wide window staged it).
# Writes per-trait filtered .clumped to a staging dir + a manifest of traits with >=1
# cis-overlap lead. This is the only FinnGen fine-mapping coloc.susie needs (matches
# Eldjarn et al., who fine-map FinnGen within the coloc regions, not genome-wide).
# Usage: prep_finngen_cis_loci.R <manifest_out> [finngen|other]

suppressPackageStartupMessages(library(data.table))
source(path.expand("code/lib/common.R"))   # write_if_changed, GWAS_BASE, CLUMP_SUBDIR

INFO_DIR      <- path.expand("./info")
args     <- commandArgs(trailingOnly = TRUE)
MANIFEST <- args[1]
COHORT   <- if (length(args) >= 2) args[2] else "finngen"
if (is.na(MANIFEST)) stop("usage: prep_finngen_cis_loci.R <manifest_out> [finngen|other]")
if (!COHORT %in% c("finngen", "other")) stop("cohort must be finngen or other, got: ", COHORT)

FINNGEN_CLUMP <- file.path(GWAS_BASE, COHORT, "clump", CLUMP_SUBDIR)
STAGE_DIR     <- file.path(GWAS_BASE, COHORT, "finemapping/cis_clumped")

dir.create(STAGE_DIR, recursive = TRUE, showWarnings = FALSE)

# bf2 protein cis windows (gene +-500kb), one per subunit (assay_cis_windows.csv), +- FM_WINDOW_BP.
source(path.expand("code/lib/cis_windows.R"))
bf2  <- fread(file.path(INFO_DIR, "bf_protein/ensg/assay_ensg.csv"))[cohort == "bf2", protein_assay]
win  <- unique(read_cis_windows(REPO_ROOT)[protein_assay %chin% bf2,
                                           .(chr = as.character(chr), start = pmax(0L, cis_start - FM_WINDOW_BP),
                                             end = cis_end + FM_WINDOW_BP)])
win  <- win[!is.na(chr) & !is.na(start) & !is.na(end)]
setkey(win, chr, start, end)

files <- list.files(FINNGEN_CLUMP, pattern = "\\.clumped$", full.names = TRUE)
if (!length(files)) stop("no .clumped files in ", FINNGEN_CLUMP)
man <- character(0)
for (f in files) {
  if (file.size(f) == 0) next
  # Stop, not skip: an unreadable file would leave the trait unstaged without a word.
  cl <- tryCatch(fread(f), error = function(e) stop("cannot read ", f, ": ", conditionMessage(e)))
  if (!nrow(cl)) next
  q <- cl[, .(chr = as.character(CHR), start = BP, end = BP, rid = .I)]
  setkey(q, chr, start, end)
  ov <- foverlaps(q, win, type = "within", nomatch = 0L)
  if (!nrow(ov)) next
  keep <- cl[sort(unique(ov$rid))]
  out  <- file.path(STAGE_DIR, basename(f))
  write_if_changed(keep, out)
  man <- c(man, out)
}
writeLines(man, MANIFEST)
cat(COHORT, "traits with >=1 cis-overlap lead:", length(man),
    "-> ", MANIFEST, "\n")
