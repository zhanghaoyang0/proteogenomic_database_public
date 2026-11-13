# Shared R helpers and constants for code. Source with
#   source(path.expand("code/lib/common.R"))
# Dependency-free on purpose: the app and projects/ read the coloc constants and
# coloc_reported() from here (sys.source into an env), so nothing here may load a package.

REPO_ROOT <- path.expand(".")
GWAS_BASE <- path.expand("~/db/shared_data/gwas")
# The clump set that annotation, fine-mapping, coloc and MR use (p 5e-8, 1000 kb,
# r2 0.01). The PRS sets (250 kb / 0.1) are separate and named in prs/.
CLUMP_SUBDIR <- "output_5e-8Pthres_1000kb_0.01r2"
COLOC_ROOT <- path.expand("~/db/shared_data/coloc")
MR_ROOT    <- path.expand("~/db/shared_data/mr")
# The LD reference (1000 Genomes phase 3 EUR, hg19; PLINK bfile prefix) for clump,
# fine-mapping, dedup and the locus plots.
LD_REF     <- path.expand("~/db/shared_data/ref_geno/1kg_v3/EUR")

P_GWS        <- 5e-8      # genome-wide significance: clump leads, kept credible sets
FM_WINDOW_BP <- 500000L   # fine-mapping half-window around a clump lead (finemap 01, 05)
PP4_STRONG   <- 0.8       # strong colocalisation: plots, app, manuscript
CIS_HALF_BP  <- 5e5L      # coloc: loci compared when their CS spans lie this close (coloc/loci_near)

# The coloc rows that are reported -- retained (PP4 >= PP4_KEEP), both compared SuSiE
# components kept GWS credible sets and both fits converged (coloc/README.md). The app
# and the manuscript call this; restating it is how the copies drifted before.
coloc_reported <- function(status, primary, PP4)
  as.character(status) == "ok" & primary %in% TRUE & is.finite(as.numeric(PP4))

# Elapsed-time logger. `t0` is the caller's start time (set per script before use),
# resolved lexically at call time.
msg <- function(...) message(sprintf("[%.1fs] ", proc.time()[[3L]] - t0), ...)

# Clump output dirs are named after the clump parameters, so a parameter change -- or a
# stale hardcoded subdir name -- silently points at a directory that does not exist.
# list.files()/file.exists() then report "nothing there", which downstream reads as
# "this trait had no significant SNPs" rather than as a misconfiguration. Fail loudly.
assert_clump_dirs <- function(...) {
  d <- unlist(list(...), use.names = FALSE)
  miss <- d[!dir.exists(d)]
  if (length(miss)) stop("clump dir not found: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

# Write a table only when its content changed, so the file's mtime moves only when the
# content does. Fine-mapping skips a trait whose output is newer than its staged
# .clumped (lib/common.sh fm_is_current); a prep step that rewrote every staged file
# on each run would make every trait look stale.
write_if_changed <- function(dt, out, sep = "\t") {
  tmp <- paste0(out, ".tmp", Sys.getpid())
  data.table::fwrite(dt, tmp, sep = sep)
  if (file.exists(out) && identical(unname(tools::md5sum(tmp)), unname(tools::md5sum(out)))) {
    unlink(tmp)
  } else if (!file.rename(tmp, out)) {
    stop("could not write ", out)
  }
  invisible(out)
}
