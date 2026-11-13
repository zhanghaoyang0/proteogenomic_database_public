#!/usr/bin/env Rscript
# Stage cis-eQTL loci for fine-mapping, restricted to COGNATE genes of tested proteins
# (cis pQTL-eQTL coloc: a protein's cis-pQTL vs the eQTL of its OWN gene). For each
# (tissue, ENSG) where ENSG encodes a tested protein, write a 1-lead .clumped
# (lead = top cis-eQTL SNP = min-P row of the clean file) so 01_finemapping.R can
# fine-map it over lead +-500kb -> susie_rds for coloc.susie (aligned by rsid).
#
# Serves all three eQTL resources, which share the cognate scope (since 2026-10-01
# MetaBrain too; before, it was staged by OVERLAP with any protein CS, cis or trans,
# in prep_metabrain_cis_loci.R, so its coloc asked a different question):
#   gtex       GTEx multi-tissue eQTL
#   sb         SingleBrain single-nucleus eQTL ("tissue" is a cell type)
#   metabrain  MetaBrain bulk brain eQTL
#
# Usage: prep_eqtl_cis_loci.R <gtex|sb|metabrain> <manifest_out>

suppressPackageStartupMessages({ library(data.table); library(parallel) })

args     <- commandArgs(trailingOnly = TRUE)
CHAIN    <- args[1]
MANIFEST <- args[2]
RESOURCE <- c(gtex = "gtex", sb = "singlebrain", metabrain = "metabrain")[CHAIN]   # dir under gwas/ and info/
COLOC    <- c(gtex = "gtex", sb = "sb", metabrain = "mb")[CHAIN]                   # coloc/protein_x_<COLOC>
if (is.na(RESOURCE) || is.na(MANIFEST))
  stop("usage: prep_eqtl_cis_loci.R <gtex|sb|metabrain> <manifest_out>")

source(path.expand("code/lib/common.R"))   # GWAS_BASE, write_if_changed
PROT_LIST <- file.path(COLOC_ROOT, sprintf("protein_x_%s/info/protein_list.csv", COLOC))
ASSAY_ENSG <- path.expand("./info/bf_protein/ensg/assay_ensg.csv")
QTL_MANIFEST <- path.expand(sprintf("./info/%s/trait_lowest_p.csv", RESOURCE))
QTL_CLEAN <- file.path(GWAS_BASE, RESOURCE, "clean")
STAGE_DIR <- file.path(GWAS_BASE, RESOURCE, "finemapping/cis_clumped")
CORES     <- as.integer(Sys.getenv("PREP_CORES", "8"))
for (f in c(PROT_LIST, ASSAY_ENSG, QTL_MANIFEST)) if (!file.exists(f)) stop("missing input: ", f)
source(path.expand("code/lib/gwas_functions.R"))   # P_EGENE_DEFAULT
dir.create(STAGE_DIR, recursive = TRUE, showWarnings = FALSE)

# 1. Cognate ENSGs of tested proteins. exposure_id = <cohort>__<protein_assay>.
prot <- fread(PROT_LIST)
ae   <- fread(ASSAY_ENSG)
ae[, exposure_id := paste0(cohort, "__", protein_assay)]
# ensembl_gene_id can list several genes ("ENSG1,ENSG2", 128 bf2 assays): split them,
# or the whole string matches no eQTL gene and the protein silently has no cognate.
cog  <- unique(sub("\\..*$", "", trimws(unlist(strsplit(ae[exposure_id %in% prot$exposure_id, ensembl_gene_id], ",")))))
cog  <- cog[!is.na(cog) & nzchar(cog)]
cat("tested proteins:", uniqueN(prot$exposure_id), "| cognate ENSGs:", length(cog), "\n")

# 2. GTEx (tissue, ENSG) loci for cognate genes only.
gx <- fread(QTL_MANIFEST)
setnames(gx, c("REGION", "ENSG_id"), c("tissue", "ENSG"))
# One row per (tissue, ENSG). trait_lowest_p carries a row per (REGION, GENE, CHR),
# and hg38->hg19 liftover strands a few variants of a gene on another chromosome, so a
# gene can appear 2-3 times -- the real cis window plus single-SNP artifacts (measured:
# 48 genes in singlebrain, 200 in metabrain, 13,932 in gtex). Undeduplicated, the same
# trait is staged and fine-mapped more than once: duplicate manifest lines, two array
# tasks racing on one output, and the loser erroring. Since c199ec0 makes any per-trait
# error fail the whole array task, that is what turned a harmless duplicate into
# 20 coloc chunks stuck on DependencyNeverSatisfied. Keep the row with the smallest
# lowest_p: the true cis window holds the gene's strongest eQTL. The widest span is
# NOT safe -- an artifact row can span tens of Mb (LGALS9 adipose: chr1 148-222 Mb,
# P 0.06, vs its chr17 cis window, P 2e-10), and keeping it dropped 1,258 GTEx
# eGenes from the manifest (2026-10-05).
cand <- gx[GENE %in% cog]
setorder(cand, tissue, GENE, lowest_p)
cand <- unique(cand, by = c("tissue", "GENE"))
# The shared eGene threshold (5e-8 since 2026-10-04). Clean trees built under the old
# 1e-5 still hold weaker genes; their credible sets would fail the 5e-8 gate in
# 01_finemapping.R anyway, so drop them before fine-mapping.
n_cand <- nrow(cand)
cand <- cand[lowest_p < P_EGENE_DEFAULT]
cat("eGene filter P <", P_EGENE_DEFAULT, ":", nrow(cand), "of", n_cand, "loci kept\n")
cat("candidate", RESOURCE, "(tissue,ENSG) cognate loci:", nrow(cand), "/", nrow(gx),
    "(deduplicated to one row per tissue+gene)\n")

# 3. Per candidate: lead = min-P SNP from its clean file; write 1-lead .clumped.
qtl_path <- function(t, e) file.path(QTL_CLEAN, t, paste0(e, ".tsv.gz"))  # clean/<tissue>/<ENSG>.tsv.gz
# Each candidate returns its staged path, or "" when it has no clean data. Anything else
# -- a try-error from a failed fork, NULL from a killed one -- stops the run: unlisted,
# an error became a manifest line and a killed fork silently dropped its genes.
res <- mclapply(seq_len(nrow(cand)), function(j) {
  t <- cand$tissue[j]; e <- cand$ENSG[j]
  cf <- qtl_path(t, e); if (!file.exists(cf)) return("")
  d  <- fread(cf, select = c("CHR", "POS", "SNP", "P"))
  if (!nrow(d)) return("")
  lead <- d[which.min(P)]
  out  <- file.path(STAGE_DIR, paste0(t, "_", e, ".clumped"))
  write_if_changed(data.table(CHR = lead$CHR, SNP = lead$SNP, BP = lead$POS, P = lead$P), out)
  out
}, mc.cores = CORES)
bad <- which(!vapply(res, function(x) is.character(x) && length(x) == 1L, logical(1)))
if (length(bad)) stop(length(bad), " candidates failed, e.g. ", cand$tissue[bad[1]], " ", cand$ENSG[bad[1]],
                      ": ", if (inherits(res[[bad[1]]], "try-error")) as.character(res[[bad[1]]]) else "no result")
man <- unlist(res); man <- man[nzchar(man)]
cat("candidates without a clean file:", sum(!file.exists(qtl_path(cand$tissue, cand$ENSG))), "/", nrow(cand), "\n")
if (nrow(cand) > 0 && !length(man)) stop("0 of ", nrow(cand), " candidates staged; check ", QTL_CLEAN)
writeLines(man, MANIFEST)
cat("staged", length(man), RESOURCE, "cognate loci ->", MANIFEST, "\n")
