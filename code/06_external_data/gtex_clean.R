#!/usr/bin/env Rscript
# Clean GTEx V8 cis-eQTL full sumstats (eQTL Catalogue reprocessing) into the
# standard schema, one tissue (dataset) per run.
# Output: {sample_group}_{ENSG}.tsv.gz, cols CHR POS SNP EA NEA EA_FREQ BETA SE P N (hg19).
# Dataset chosen by CLI arg (sample_group or QTD id) or SLURM_ARRAY_TASK_ID
# (1-based row of info/gtex/datasets.csv).
library(data.table)
library(GenomicRanges)
library(rtracklayer)

project_root <- "."
setwd(project_root)
source(file.path(project_root, "code/lib/gwas_functions.R"))

# Directories
path_gtex  <- path.expand("~/db/shared_data/gwas/gtex")
path_raw   <- file.path(path_gtex, "raw")
path_clean <- file.path(path_gtex, "clean")
dir.create(path_clean, recursive = TRUE, showWarnings = FALSE)
path_info_repo <- file.path(project_root, "info", "gtex")
path_summary   <- file.path(path_info_repo, "trait_lowest_p")   # per-tissue, array-safe
dir.create(path_summary, recursive = TRUE, showWarnings = FALSE)
path_audit     <- file.path(path_info_repo, "egene_audit")      # ALL tested genes
dir.create(path_audit, recursive = TRUE, showWarnings = FALSE)
P_EGENE <- P_EGENE_DEFAULT   # shared eGene rule, see code/lib/gwas_functions.R

# Dataset manifest (49 GTEx ge datasets)
manifest <- fread(file.path(path_info_repo, "datasets.csv"))

# --- select dataset ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
slurm_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (length(args) >= 1 && nzchar(args[1])) {
  row <- manifest[sample_group == args[1] | dataset_id == args[1]]
} else if (nzchar(slurm_id)) {
  row <- manifest[as.integer(slurm_id)]
} else {
  stop("Provide a sample_group or QTD id as arg, or run under SLURM_ARRAY_TASK_ID")
}
if (nrow(row) != 1) stop("Could not resolve a single dataset from args/array id")
tissue <- row$sample_group
qtd    <- row$dataset_id
cat("Dataset:", qtd, "| tissue:", tissue, "| sample_size:", row$sample_size, "\n")

# reference afreq (hg19) for harmonization + liftover chain hg38 -> hg19
freq <- fread(path.expand("~/db/shared_data/ref_geno/1kg_v3/EUR.afreq"))
chain_hg38_to_hg19 <- rtracklayer::import(
  path.expand("~/db/shared_data/ref_geno/liftover/hg38ToHg19.over.chain"), format = "chain")

# --- read full cis sumstats for this dataset ---------------------------------
file_in <- file.path(path_raw, paste0(qtd, ".all.tsv.gz"))
if (!file.exists(file_in)) stop("Missing sumstats file: ", file_in)
df <- fread(file_in, select = c("gene_id", "chromosome", "position", "ref", "alt",
                                "beta", "se", "pvalue", "ac", "an", "rsid"))

# --- map to standard schema --------------------------------------------------
# ref->NEA, alt->EA, beta is per-alt; EA_FREQ = alt-allele freq = ac/an; N = an/2
setnames(df, c("chromosome", "position", "ref", "alt", "beta", "se", "pvalue", "rsid"),
             c("CHR",        "POS",      "NEA", "EA",  "BETA", "SE", "P",      "SNP"))
df[, `:=`(
  CHR     = suppressWarnings(as.integer(CHR)),      # X/Y -> NA, dropped below
  POS     = as.integer(POS),
  GENE    = sub("\\..*$", "", gene_id),             # ENSG (already unversioned)
  EA_FREQ = fifelse(an > 0, ac / an, NA_real_),
  N       = as.integer(round(an / 2))
)]
df[, ENSG_id := GENE]
df <- df[!is.na(CHR) & !is.na(POS) & !is.na(EA) & !is.na(NEA) & !is.na(EA_FREQ)]
if (nrow(df) == 0) stop("No parseable autosomal rows for ", tissue)

# --- audit: best P for EVERY tested gene, written BEFORE the eGene filter -----
# The trait manifest below only covers kept genes, so without this there is no
# record of what the filter dropped, and "gene X has no eQTL in this tissue"
# cannot be distinguished from "gene X was dropped at P = 3e-5". metabrain and
# singlebrain keep the same record; check it before writing any absence claim.
df[, REGION := tissue]
gene_minp <- df[, .(lowest_p = min(P, na.rm = TRUE)), by = .(REGION, GENE)]
gene_minp[, kept := lowest_p < P_EGENE]
setorder(gene_minp, lowest_p)
fwrite(gene_minp, file.path(path_audit, paste0(tissue, ".csv")))

# --- eGenes: genes with a cis signal, keep ALL their SNPs ---------------------
# All SNPs, not just significant ones: SuSiE/coloc need the full regional profile,
# so the cut is per gene, never per SNP. P_EGENE_DEFAULT is the shared eGene rule
# in code/lib/gwas_functions.R, used by metabrain/gtex/singlebrain alike.
genes_to_save <- gene_minp[kept == TRUE, GENE]
cat(sprintf("eGenes at min(P) < %g: %d of %d tested\n",
            P_EGENE, length(genes_to_save), nrow(gene_minp)))
if (length(genes_to_save) == 0) stop("No signal gene for ", tissue, " (unexpected)")
df <- df[GENE %in% genes_to_save]

# --- liftover hg38 -> hg19 ---------------------------------------------------
df <- lift_to_hg19(df, chain_hg38_to_hg19)
if (nrow(df) == 0) stop("Nothing lifted to hg19 for ", tissue)

# Per-(gene, chr) lowest P and hg19 pos range for the trait manifest
trait_summary <- df[, .(
  lowest_p = min(P, na.rm = TRUE),
  start = min(POS), end = max(POS),
  ENSG_id = .SD[which.min(P), ENSG_id]
), by = .(REGION, GENE, CHR)]
setorder(trait_summary, GENE, CHR)
fwrite(trait_summary, file.path(path_summary, paste0(tissue, ".csv")))

# --- per-gene harmonize + write ----------------------------------------------
# One file per gene, nested in a per-tissue subdir: clean/{tissue}/{ENSG}.tsv.gz
# (49 dirs of a few-k files each, vs ~250k files in one flat dir). Per-gene
# granularity is the coloc unit, so reads stay O(1) and skip-existing works per gene.
tissue_dir <- file.path(path_clean, tissue)
dir.create(tissue_dir, recursive = TRUE, showWarnings = FALSE)
for (gene in unique(df$GENE)) {
  file_out <- file.path(tissue_dir, paste0(gene, ".tsv.gz"))
  if (file.exists(file_out)) next
  df_gene <- df[GENE == gene, .(CHR, POS, SNP, EA, NEA, EA_FREQ, BETA, SE, P, N)]
  # harmonize: align alleles to 1KG ALT/REF, (re)assign rsID, drop palindromes, MAF filter
  # (pass freq directly, as metabrain_clean.R does: harmonize's setnames is idempotent and
  #  it never mutates freq rows, so no per-gene copy of the multi-million-row table is needed)
  df_gene <- round_gwas_clean(harmonize_to_1kg(df_gene, freq))
  if (nrow(df_gene) == 0) next
  write_atomic(df_gene, file_out)   # a killed job leaves no truncated file for the skip above
}

# --- fail-loud: assert outputs exist -----------------------------------------
produced <- list.files(tissue_dir, pattern = "^ENSG.*\\.tsv\\.gz$")
if (length(produced) == 0) stop("No clean files written for ", tissue)
message(sprintf("GTEx clean done: %s (%s) | signal genes=%d, gene files=%d",
                tissue, qtd, length(genes_to_save), length(produced)))
