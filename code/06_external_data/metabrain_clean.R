#!/usr/bin/env Rscript
library(data.table)
library(stringr)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
# Set fixed project root for stable batch execution
project_root <- "."
setwd(project_root)
source(file.path(project_root, "code/lib/gwas_functions.R"))

# Directories
path_MetaBrain <- path.expand("~/db/shared_data/gwas/metabrain")
path_raw <- file.path(path_MetaBrain, "raw")
path_clean <- file.path(path_MetaBrain, "clean")
path_info <- file.path(path_MetaBrain, "info")
# Manifest lives in project info/ for repo-tracked discoverability
path_info_repo <- file.path(project_root, "info", "metabrain")
dir.create(path_info_repo, recursive = TRUE, showWarnings = FALSE)

# reference afreq for harmonization (hg19)
file_freq <- path.expand("~/db/shared_data/ref_geno/1kg_v3/EUR.afreq")
freq <- fread(file_freq)

# liftover chain: hg38 -> hg19
chain_file <- path.expand("~/db/shared_data/ref_geno/liftover/hg38ToHg19.over.chain")
chain_hg38_to_hg19 <- rtracklayer::import(chain_file, format = "chain")

# Find all files
files <- list.files(path_raw)
files <- files[grepl('chr', files) & !grepl('\\.tbi$', files)]

P_EGENE <- P_EGENE_DEFAULT   # shared eGene rule, see code/lib/gwas_functions.R

trait_summary_list <- list()

for (file_in in files) {
  # Format: YYYY-MM-DD-{region}-...-{chr}.txt.gz
  match_result <- str_match(file_in, "^\\d{4}-\\d{2}-\\d{2}-([^-]+)-.*-(chr[^.]+)\\.txt\\.gz$")
  region <- match_result[1, 2]
  chr <- match_result[1, 3]
  cat("Processing", file_in, "-> region:", region, "chr:", chr, "\n")

  df_raw <- fread(file.path(path_raw, file_in), sep = "\t", colClasses = "character")
  if (nrow(df_raw) == 0) next
  df <- as.data.table(df_raw)
  df[, c("a1", "a2") := tstrsplit(SNPAlleles, "/", fixed = TRUE)]
  setnames(df, c("SNPChr", "SNPPos", "SNPEffectAllele", "SNPEffectAlleleFreq", "MetaBeta", "MetaSE", "MetaP", "MetaPN", "Gene"),
           c("CHR", "POS", "EA", "EA_FREQ", "BETA", "SE", "P", "N", "GENE"))
  df[, NEA := fifelse(a1 == EA, a2, fifelse(a2 == EA, a1, NA_character_))]
  # Use ENSG gene ID and remove version suffix (e.g. ENSG000001.10 -> ENSG000001)
  df[, GENE := sub("\\..*$", "", GENE)]
  df[, ENSG_id := GENE]
  df <- df[, .(CHR, POS, SNP, EA, NEA, EA_FREQ, BETA, SE, P, N, GENE, ENSG_id)]
  df[, `:=`(
    CHR = as.integer(CHR),
    POS = as.integer(POS),
    EA_FREQ = as.numeric(EA_FREQ),
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P),
    N = as.integer(N)
  )]
  df <- df[!is.na(CHR) & !is.na(POS) & !is.na(SNP) & SNP != "" & !is.na(GENE) & GENE != ""]
  if (nrow(df) == 0) next

  # liftover SNP positions from hg38 to hg19 before harmonization
  df <- lift_to_hg19(df, chain_hg38_to_hg19)
  if (nrow(df) == 0) next

  # Per (region, gene, chr) lowest P and pos range for trait summary
  df[, REGION := region]
  df_trait_min <- df[, .(
    lowest_p = min(P, na.rm = TRUE),
    start = min(POS),
    end = max(POS),
    ENSG_id = .SD[which.min(P), ENSG_id]
  ), by = .(REGION, GENE, CHR)]
  trait_summary_list[[length(trait_summary_list) + 1]] <- df_trait_min[, .(REGION, GENE, CHR, lowest_p, start, end, ENSG_id)]

  # Genes to save: min(P) < P_EGENE in this file; save ALL SNPs for those genes
  # (from this file, no per-SNP P filter). P_EGENE_DEFAULT is shared with
  # gtex_clean.R / singlebrain_clean.R via code/lib/gwas_functions.R.
  genes_to_save <- df[, .(min_p = min(P, na.rm = TRUE)), by = GENE][min_p < P_EGENE, GENE]
  if (length(genes_to_save) == 0) {
    cat(sprintf("  No gene with min(P) < %g, skip\n", P_EGENE))
    next
  }

  cat(sprintf("  Genes with min(P) < %g: %d, saving full eQTL for each (this chr)\n", P_EGENE, length(genes_to_save)))
  n_written <- 0L
  for (gene in genes_to_save) {
    # clean/{region}/{ENSG}.tsv.gz, the layout fine-mapping and coloc read (as gtex_clean.R).
    # Checked before harmonizing, and written atomically, so a rerun redoes only
    # missing genes and never keeps a half-written file.
    path_gz <- file.path(path_clean, region, paste0(gene, ".tsv.gz"))
    if (file.exists(path_gz)) next
    df_gene <- df[GENE == gene][, .(CHR, POS, SNP, EA, NEA, EA_FREQ, BETA, SE, P, N)]
    df_gene <- round_gwas_clean(harmonize_to_1kg(df_gene, freq))
    write_atomic(df_gene, path_gz)
    n_written <- n_written + 1L
  }
  cat(sprintf("  Wrote %d new gene files for %s %s\n", n_written, region, chr))
}

if (length(trait_summary_list) > 0) {
  trait_summary <- rbindlist(trait_summary_list, use.names = TRUE, fill = TRUE)
  trait_summary <- trait_summary[, .(
    lowest_p = min(lowest_p, na.rm = TRUE),
    start = min(start),
    end = max(end),
    ENSG_id = fcoalesce(ENSG_id[!is.na(ENSG_id) & ENSG_id != ""][1L], NA_character_)
  ), by = .(REGION, GENE, CHR)]
  setorder(trait_summary, REGION, GENE, CHR)
  fwrite(trait_summary, file.path(path_info_repo, "trait_lowest_p.csv"))
}

# ============================================================
# Check annotation results
# ============================================================
trait_summary <- fread(file.path(path_info_repo, "trait_lowest_p.csv"))
expected <- unique(trait_summary[lowest_p < P_EGENE, .(REGION, GENE)])
message(sprintf("Expected %d genes with min(P) < %g (of %d tested)",
                nrow(expected), P_EGENE, uniqueN(trait_summary[, .(REGION, GENE)])))