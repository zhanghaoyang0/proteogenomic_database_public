#!/usr/bin/env Rscript
# Clean SingleBrain single-nucleus cis-eQTL sumstats into the standard schema,
# one cell type per run (mirrors gtex_clean.R).
# Output: clean/{cell_type}/{ENSG}.tsv.gz, cols CHR POS SNP EA NEA EA_FREQ BETA SE P N (hg19).
#
# SingleBrain specifics:
#  - hg38 -> lifted to hg19 (same chain file as metabrain_clean.R / gtex_clean.R)
#  - variant_id is already an rsID -> SNP maps directly
#  - no allele frequency column -> imputed from 1kg EUR (harmonize impute_frq = TRUE)
#  - no per-SNP N -> constant n_donors from info/singlebrain/datasets.csv
#  - eGenes = shared min(Fixed_P) < P_EGENE_DEFAULT rule (see gwas_functions.R)
#  - BETA/SE/P come from the FIXED-effect meta (fixed_beta/fixed_sd/Fixed_P): the
#    random-effect meta reports only Z and P, and coloc/SuSiE need beta + se.
#
# Cell type chosen by CLI arg or SLURM_ARRAY_TASK_ID (1-based row of datasets.csv).
library(data.table)
library(GenomicRanges)
library(rtracklayer)

project_root <- "."
setwd(project_root)
source(file.path(project_root, "code/lib/gwas_functions.R"))

# Directories
path_sb    <- path.expand("~/db/shared_data/gwas/singlebrain")
path_raw   <- file.path(path_sb, "raw")
path_clean <- file.path(path_sb, "clean")
dir.create(path_clean, recursive = TRUE, showWarnings = FALSE)
path_info_repo <- file.path(project_root, "info", "singlebrain")
path_summary   <- file.path(path_info_repo, "trait_lowest_p")   # per-cell-type, array-safe
dir.create(path_summary, recursive = TRUE, showWarnings = FALSE)
path_audit     <- file.path(path_info_repo, "egene_audit")     # per-gene min P + kept flag
dir.create(path_audit, recursive = TRUE, showWarnings = FALSE)
path_tmp <- file.path(project_root, "tmp/gwas")                 # gene lists (gitignored)
dir.create(path_tmp, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(file.path(path_info_repo, "datasets.csv"))

# eGene rule: shared with metabrain_clean.R / gtex_clean.R via
# code/lib/gwas_functions.R, so the three eQTL resources stay comparable.
# Applied to the FIXED-effect P, which is the one coloc uses.
#
# Note this is stricter than SingleBrain's own qval < 0.05 eGene call (MG: 3,049
# genes at Fixed_P<1e-5, the threshold before 2026-10-04, vs 8,423 at qval<0.05;
# now P_EGENE_DEFAULT = 5e-8 keeps fewer still). Cross-resource comparability wins
# here: reporting a gene as cell-type-specific requires the bulk resources to have
# been screened by the same rule. info/singlebrain/egene_audit/ records every
# tested gene's P and qval regardless, so the dropped ones remain checkable.
P_EGENE <- P_EGENE_DEFAULT

# --- select cell type --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
slurm_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
if (length(args) >= 1 && nzchar(args[1])) {
  row <- manifest[cell_type == args[1]]
} else if (nzchar(slurm_id)) {
  row <- manifest[as.integer(slurm_id)]
} else {
  stop("Provide a cell_type as arg, or run under SLURM_ARRAY_TASK_ID")
}
if (nrow(row) != 1) stop("Could not resolve a single cell type from args/array id")
ct <- row$cell_type
n_donors <- as.integer(row$n_donors)
cat("Cell type:", ct, "|", row$label, "| n_donors:", n_donors, "\n")

file_top  <- file.path(path_raw, paste0(ct, "_eqtl_top_assoc.tsv.gz"))
file_full <- file.path(path_raw, paste0(ct, "_eqtl_full_assoc.tsv.gz"))
for (f in c(file_top, file_full)) if (!file.exists(f)) stop("Missing sumstats file: ", f)

# reference afreq (hg19) for harmonization + liftover chain hg38 -> hg19
freq <- fread(path.expand("~/db/shared_data/ref_geno/1kg_v3/EUR.afreq"))
chain_hg38_to_hg19 <- rtracklayer::import(
  path.expand("~/db/shared_data/ref_geno/liftover/hg38ToHg19.over.chain"), format = "chain")

# --- eGenes from the small top_assoc file ------------------------------------
# The full_assoc file is ~5 GB gzipped (~10^8 rows) and cannot be read whole, so
# the eGene set is decided here, on the one-row-per-feature top file, and used to
# prefilter the stream. qval lives only in top_assoc, which is why the call is made
# at this stage rather than after reading the full data.
# Prefilter deliberately looser (1e-3) than P_EGENE: top_assoc picks its row by the
# random-effect adjusted P, so the Fixed_P shown there is not necessarily the gene's
# minimum. A generous candidate set keeps the exact min(Fixed_P) < P_EGENE rule below
# authoritative -- the prefilter can only ever be looser, never tighter.
top <- fread(file_top, select = c("feature", "Fixed_P", "Random_P", "qval"))
cand <- top[qval < 0.05 | Random_P < 1e-3 | Fixed_P < 1e-3, unique(feature)]
cat(sprintf("Features tested: %d | candidates streamed: %d\n", nrow(top), length(cand)))
if (length(cand) == 0) stop("No candidate gene for ", ct, " (unexpected)")

# Audit over ALL tested genes (top_assoc is already one row per gene), written
# before any filtering so absence claims stay answerable. Same role as the
# metabrain / gtex manifests.
fwrite(top, file.path(path_audit, paste0(ct, "_all_tested.csv")))

file_cand <- file.path(path_tmp, paste0("singlebrain_cand_", ct, ".txt"))
writeLines(cand, file_cand)

# --- stream full_assoc: keep candidate features and the 10 needed columns -----
awk_prog <- paste0(
  "NR==FNR { g[$1]=1; next } ",
  "FNR==1 { for (i=1; i<=NF; i++) h[$i]=i; ",
  "print \"feature\",\"SNP\",\"CHR\",\"POS\",\"ref\",\"alt\",\"EA\",\"BETA\",\"SE\",\"P\"; next } ",
  "($h[\"feature\"] in g) && ($h[\"fixed_sd\"]+0) > 0 { ",
  "print $h[\"feature\"],$h[\"variant_id\"],$h[\"chr\"],$h[\"pos\"],",
  "$h[\"ref\"],$h[\"alt\"],$h[\"Allele\"],$h[\"fixed_beta\"],$h[\"fixed_sd\"],$h[\"Fixed_P\"] }")
# fixed_sd > 0 drops mmQTL's placeholder rows. Where a variant has too few
# contributing cohorts (1-2 of 4) mmQTL emits fixed_beta = fixed_sd = 0, Fixed_P = 1
# -- and Random_P = 1 on the very same rows, so these carry no information under
# either model, they are not a fixed-vs-random trade. Measured on Ast: 46.7% of all
# rows. They must go: SE = 0 makes z = beta/se undefined and would break SuSiE
# downstream. Dropping them in awk also halves what R has to hold. eGene selection
# is unaffected -- a placeholder has P = 1 and can never be a gene's minimum.
cmd <- sprintf("gzip -cd %s | awk -F'\t' -v OFS='\t' %s %s -",
               shQuote(file_full), shQuote(awk_prog), shQuote(file_cand))
cat("Streaming full_assoc for", length(cand), "eGenes ...\n")
df <- fread(cmd = cmd, showProgress = FALSE)
if (nrow(df) == 0) stop("awk prefilter returned no rows for ", ct,
                        " -- check column names in ", file_full)
cat("Rows after gene prefilter:", nrow(df), "\n")

# --- map to standard schema --------------------------------------------------
# Allele is the effect allele; NEA is whichever of ref/alt it is not.
df[, `:=`(
  CHR  = suppressWarnings(as.integer(sub("^chr", "", CHR))),   # X/Y -> NA, dropped below
  POS  = as.integer(POS),
  BETA = as.numeric(BETA),
  SE   = as.numeric(SE),
  P    = as.numeric(P),
  GENE = sub("\\..*$", "", feature),                           # ENSG00000186047.11 -> ENSG00000186047
  N    = n_donors
)]
df[, NEA := fifelse(EA == alt, ref, fifelse(EA == ref, alt, NA_character_))]
df[, ENSG_id := GENE]
df <- df[!is.na(CHR) & !is.na(POS) & !is.na(SNP) & SNP != "" &
         !is.na(EA) & !is.na(NEA) & !is.na(BETA) & !is.na(SE) & !is.na(P)]
# Belt-and-braces for the awk fixed_sd filter above: SE must be a usable positive
# number or SuSiE gets a non-finite z. Loud, because a non-zero count here means the
# awk condition stopped working.
n_bad_se <- df[!is.finite(SE) | SE <= 0, .N]
if (n_bad_se > 0) {
  warning("Dropping ", n_bad_se, " rows with SE <= 0 that the awk filter should have removed")
  df <- df[is.finite(SE) & SE > 0]
}
if (nrow(df) == 0) stop("No parseable autosomal rows for ", ct)

# --- eGenes: exact min(Fixed_P) < P_EGENE rule on the full data ---------------
# Retained genes keep ALL their cis SNPs, not just significant ones: SuSiE and
# coloc need the full regional profile, so the cut is per gene, never per SNP.
df[, REGION := ct]
gene_minp <- df[, .(lowest_p = min(P, na.rm = TRUE)), by = .(REGION, GENE)]
gene_minp[, kept := lowest_p < P_EGENE]
setorder(gene_minp, lowest_p)
fwrite(gene_minp, file.path(path_audit, paste0(ct, "_streamed.csv")))
genes_to_save <- gene_minp[kept == TRUE, GENE]
cat(sprintf("eGenes kept: %d of %d (P_EGENE = %g; dropped %d)\n",
            length(genes_to_save), nrow(gene_minp), P_EGENE,
            nrow(gene_minp) - length(genes_to_save)))
if (length(genes_to_save) == 0) stop("No signal gene for ", ct, " (unexpected)")
df <- df[GENE %in% genes_to_save]

# --- liftover hg38 -> hg19 ---------------------------------------------------
df <- lift_to_hg19(df, chain_hg38_to_hg19)
if (nrow(df) == 0) stop("Nothing lifted to hg19 for ", ct)

# Per-(gene, chr) lowest P and hg19 pos range -- same schema/build as the
# metabrain and gtex manifests, which prep_*_cis_loci.R reads as hg19.
trait_summary <- df[, .(
  lowest_p = min(P, na.rm = TRUE),
  start = min(POS), end = max(POS),
  ENSG_id = .SD[which.min(P), ENSG_id]
), by = .(REGION, GENE, CHR)]
setorder(trait_summary, GENE, CHR)
fwrite(trait_summary, file.path(path_summary, paste0(ct, ".csv")))

# --- per-gene harmonize + write ----------------------------------------------
# One file per gene under a per-cell-type subdir: clean/{cell_type}/{ENSG}.tsv.gz,
# the same layout metabrain/gtex use, so finemap + coloc read it unchanged.
ct_dir <- file.path(path_clean, ct)
dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

# Restrict the 1kg reference to each gene's own cis window before harmonizing.
# harmonize_to_1kg merges its input against the whole ~9M-row reference twice (once
# by rsID for the frequency, once by CHR:POS for allele alignment). At ~5,000 genes
# per cell type that was the dominant cost -- measured 44 gene files in 11 minutes,
# i.e. ~8 h for one cell type, uncomfortably close to the 12 h wall. A gene spans
# ~2 Mb, so a per-chromosome split plus a POS range subset leaves a few thousand
# reference rows and the merges become trivial. The result is identical: the
# reference rows dropped could never have matched this gene's variants.
setnames(freq, c("#CHROM", "ID", "ALT_FREQS"), c("CHR", "SNP", "ALT_FREQ"), skip_absent = TRUE)
freq[, `:=`(CHR = as.integer(CHR), POS = as.integer(POS))]
freq_by_chr <- split(freq, freq$CHR)

setorder(df, GENE)
for (gene in unique(df$GENE)) {
  file_out <- file.path(ct_dir, paste0(gene, ".tsv.gz"))
  if (file.exists(file_out)) next
  df_gene <- df[GENE == gene, .(CHR, POS, SNP, EA, NEA, BETA, SE, P, N)]
  fs <- freq_by_chr[[as.character(df_gene$CHR[1])]]
  if (is.null(fs)) next
  rng <- range(df_gene$POS)
  freq_sub <- fs[POS >= rng[1] & POS <= rng[2]]
  if (nrow(freq_sub) == 0) next
  # impute_frq = TRUE: SingleBrain ships no allele frequency, so EA_FREQ is taken
  # from 1kg EUR by rsID before allele alignment / palindrome + MAF filtering.
  df_gene <- round_gwas_clean(harmonize_to_1kg(df_gene, freq_sub, impute_frq = TRUE))
  if (nrow(df_gene) == 0) next
  write_atomic(df_gene, file_out)   # a killed job leaves no truncated file for the skip above
}

# --- fail-loud: assert outputs exist -----------------------------------------
produced <- list.files(ct_dir, pattern = "^ENSG.*\\.tsv\\.gz$")
if (length(produced) == 0) stop("No clean files written for ", ct)
message(sprintf("SingleBrain clean done: %s | signal genes=%d, gene files=%d",
                ct, length(genes_to_save), length(produced)))
