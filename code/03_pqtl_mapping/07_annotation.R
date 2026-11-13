#!/usr/bin/env Rscript
# Annotate GWAS signals with VEP consequence info.
# Signal set = clump lead SNPs UNION finemap credible-set (CS) SNPs.
#
# Role columns (per SNP, per trait):
#   is_clump_lead  TRUE if the SNP is a clumped lead SNP for the trait.
#   is_cs_lead     TRUE if the SNP is the top-PIP lead of its credible set.
#   cs_id / pip    credible-set id and PIP (NA when not in any CS).
#   => a "cs_member" is any row with cs_id present and is_cs_lead == FALSE.
#
# Section 1: parse clumped files -> per-trait summary + clump lead sets.
# Section 2: build unified SNP table (clump leads U finemap CS SNPs) -> VCF for VEP.
# Section 3: VEP, run by 06_annotation.sbatch between Section 2 and 4.
# Section 4: read VEP output, fill is_cis, write final tables.
#
# Env vars:
#   GWAS  = bf1 | bf2 | finngen | other   (default: bf2)
#   TYPE  = protein | olink | pheno | ""  (default: protein; "" = flat clump/finemap dirs)
#
# NOTE: this step now depends on finemap output -> run it AFTER finemapping.
# If finemap files are absent (e.g. pheno), it falls back to clump leads only.

library(data.table)
options(width = 200)

GWAS     <- Sys.getenv("GWAS",  unset = "bf2")
TYPE     <- Sys.getenv("TYPE",  unset = "protein")
source(path.expand("code/lib/common.R"))   # GWAS_BASE, CLUMP_SUBDIR, REPO_ROOT
# REPO_ROOT, not dirname() of --file=: that kept the "--file=" prefix and found the repo
# only when the script ran from its own directory.
wb_info_dir  <- file.path(REPO_ROOT, "info")
BASE_DIR <- GWAS_BASE
PVAL_THR <- "5e-8"   # names the clump input dir; the output dir is CLUMP_SUBDIR

gwas_dir <- file.path(BASE_DIR, GWAS)
type_dir <- if (TYPE == "olink") "protein" else TYPE

# Paths differ: bf/type have subfolders; finngen/other are flat
if (nchar(type_dir) > 0) {
  clump_base  <- file.path(gwas_dir, "clump", type_dir)
  finemap_dir <- file.path(gwas_dir, "finemapping", type_dir)
  anno_dir    <- file.path(gwas_dir, "annotation", type_dir)
} else {
  clump_base  <- file.path(gwas_dir, "clump")
  finemap_dir <- file.path(gwas_dir, "finemapping")
  anno_dir    <- file.path(gwas_dir, "annotation")
}
clump_out_dir   <- file.path(clump_base, CLUMP_SUBDIR)
clump_input_dir <- file.path(clump_base, sprintf("input_%sPthres", PVAL_THR))
info_dir        <- file.path(gwas_dir, "info")
dir.create(anno_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(info_dir,  recursive = TRUE, showWarnings = FALSE)

out_vcf <- file.path(anno_dir, "leadsnp_for_vep.vcf")
vep_out <- file.path(anno_dir, "leadsnp_vep_anno.txt")

message(sprintf("GWAS=%s  TYPE=%s", GWAS, TYPE))
message("Clump output: ", clump_out_dir)
message("Finemap dir:  ", finemap_dir)
message("Anno dir:     ", anno_dir)

# Load protein coords for protein GWAS only
use_protein_info <- (type_dir == "protein")
if (use_protein_info) {
  cis_coords_path <- file.path(wb_info_dir, "bf_protein/ensg/assay_ensg.csv")
  if (!file.exists(cis_coords_path)) stop("Missing ", cis_coords_path)
  cis_tbl <- fread(cis_coords_path)[, .(cohort, protein_assay, gene, gene_chr, gene_start, gene_end)]
}

# ============================================================
# Section 1: Parse clumped files; build summary + clump lead sets
# ============================================================
clumped_files <- list.files(clump_out_dir, pattern = "\\.clumped$", full.names = TRUE)
if (!length(clumped_files)) stop("No .clumped files under ", clump_out_dir)
traits <- sub("\\.clumped$", "", basename(clumped_files))

summary_list <- vector("list", length(clumped_files))
clump_leads  <- vector("list", length(clumped_files))   # trait -> lead SNP vector
names(clump_leads) <- traits
for (i in seq_along(clumped_files)) {
  trait <- traits[i]
  clump <- setDT(read.table(clumped_files[i], header = TRUE))
  clump_leads[[trait]] <- clump$SNP
  input_f <- file.path(clump_input_dir, paste0(trait, ".tsv"))
  # Stop, not NA: the same file supplies the leads outside a credible set (below).
  if (!file.exists(input_f)) stop("clump input missing for ", trait, ": ", input_f)
  n_sig <- nrow(fread(input_f, sep = "\t"))
  summary_list[[i]] <- data.table(trait = trait, n_sig_signal = n_sig, n_lead = nrow(clump))
}
summary <- rbindlist(summary_list)[order(-n_lead)]

if (use_protein_info) {
  setnames(summary, "trait", "protein_assay")
  # cis_tbl keys on (cohort, raw protein_assay); split "bf2__<assay>" for the join.
  # By assay alone, an assay measured in both cohorts got two rows.
  summary[, `:=`(cohort = sub("__.*$", "", protein_assay),
                 protein_assay_raw = sub("^bf[12]__", "", protein_assay))]
  n_before <- nrow(summary)
  summary <- merge(summary, cis_tbl, by.x = c("cohort", "protein_assay_raw"),
                   by.y = c("cohort", "protein_assay"), all.x = TRUE)
  if (nrow(summary) != n_before) stop("cis_tbl join changed the summary row count")
  summary[, c("cohort", "protein_assay_raw") := NULL]
  message("Proteins with no coords: ", nrow(summary[is.na(gene_chr)]), " / ", nrow(summary))
}

# ============================================================
# Section 2: Build unified SNP table (clump leads U finemap CS SNPs); write VCF
# ============================================================
# For each trait, take all finemap CS SNPs, then add any clump lead SNP that is
# not already present in a CS (so header-only / non-finemapped loci aren't lost).
build_trait_snps <- function(trait) {
  leads   <- clump_leads[[trait]]
  fm_file <- file.path(finemap_dir, paste0(trait, ".finemapped.tsv.gz"))

  fm_rows <- NULL
  if (file.exists(fm_file)) {
    fm <- fread(fm_file)
    if (nrow(fm)) {
      fm_rows <- fm[, .(trait,
                        SNP = snp, CHR = as.character(chr), POS = pos,
                        EA = ea, NEA = nea,
                        cs_id = cs_id, cs_purity = cs_purity, pip = pip,
                        is_cis = as.logical(is_cis),
                        is_cs_lead = (snp == lead_snp))]
    }
  }

  have    <- if (!is.null(fm_rows)) unique(fm_rows$SNP) else character(0)
  missing <- setdiff(leads, have)
  extra   <- NULL
  if (length(missing)) {
    input_f <- file.path(clump_input_dir, paste0(trait, ".tsv"))
    if (!file.exists(input_f)) stop("clump input missing for ", trait, ": ", input_f)
    pv    <- fread(input_f, sep = "\t")
    extra <- pv[SNP %in% missing,
                .(trait, SNP, CHR = as.character(CHR), POS, EA, NEA,
                  cs_id = NA_integer_, cs_purity = NA_real_, pip = NA_real_,
                  is_cis = NA, is_cs_lead = FALSE)]
  }

  out <- rbindlist(list(fm_rows, extra), use.names = TRUE, fill = TRUE)
  if (nrow(out)) out[, is_clump_lead := SNP %in% leads]
  out
}

snp_dt <- rbindlist(lapply(traits, build_trait_snps), use.names = TRUE, fill = TRUE)
if (!nrow(snp_dt)) stop("No SNPs to annotate (empty clump + finemap).")
message(sprintf("Unified SNPs: %d rows | clump_lead=%d  cs_lead=%d  cs_member=%d",
                nrow(snp_dt), sum(snp_dt$is_clump_lead),
                sum(snp_dt$is_cs_lead, na.rm = TRUE),
                sum(!is.na(snp_dt$cs_id) & !snp_dt$is_cs_lead)))

# One VCF record per unique variant (REF=NEA, ALT=EA), sorted by chr/pos.
vcf_dt <- unique(snp_dt[, .(CHR, POS, SNP, NEA, EA)])
vcf_dt[, chr_num := suppressWarnings(as.integer(CHR))]
vcf_dt[is.na(chr_num), chr_num := ifelse(CHR %in% c("X", "x"), 99L, 98L)]
setorder(vcf_dt, chr_num, POS)
vcf_dt[, chr_num := NULL]

vcf_body  <- vcf_dt[, paste(CHR, POS, SNP, NEA, EA, ".", ".", ".", sep = "\t")]
vcf_lines <- c("##fileformat=VCFv4.2",
               "##source=07_annotation.R",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
               vcf_body)

# Idempotent: the sbatch runs this script twice (before and after VEP). Only
# (re)write the VCF and drop the VEP output when the variant set actually changed,
# so the post-VEP pass keeps the fresh VEP result instead of nuking it. A changed
# VCF invalidates any old VEP output (else we'd silently annotate a stale set).
vcf_unchanged <- file.exists(out_vcf) && identical(readLines(out_vcf, warn = FALSE), vcf_lines)
if (vcf_unchanged) {
  message("VCF unchanged (", nrow(vcf_dt), " variants): ", out_vcf)
} else {
  writeLines(vcf_lines, out_vcf)
  message("Wrote ", nrow(vcf_dt), " unique variants to ", out_vcf)
  if (file.exists(vep_out)) {
    file.remove(vep_out)
    message("Removed stale VEP output (VCF changed): ", vep_out)
  }
}
message("Now run VEP (06_annotation.sbatch does), then re-run this script for Section 4.")

# ============================================================
# Section 4: Read VEP output; fill is_cis; write final tables
# ============================================================
if (!file.exists(vep_out)) {
  # The sbatch's post-VEP call sets REQUIRE_VEP=1: there a missing VEP output means
  # VEP was skipped or failed, and exiting 0 would leave the old tables in place.
  if (Sys.getenv("REQUIRE_VEP") == "1") stop("VEP output not found after the VEP step: ", vep_out)
  message("VEP output not found: ", vep_out, " — skipping Section 4.")
  quit(save = "no")
}

lines      <- readLines(vep_out, warn = FALSE)
header     <- sub("^#", "", tail(lines[grepl("^#[^#]", lines)], 1))
data_lines <- lines[!is.na(lines) & !grepl("^#", lines) & nzchar(lines)]
if (!length(header)) stop("Missing VEP header in ", vep_out)
vep_cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
if (length(data_lines)) {
  # header=FALSE: the column header was already stripped above; without it fread's
  # auto-detection can mistake the first data row for a header and drop one variant.
  vep_dt <- fread(text = paste(data_lines, collapse = "\n"),
                  header = FALSE, col.names = vep_cols)
  vep_dt <- unique(vep_dt)
} else {
  vep_dt <- as.data.table(setNames(rep(list(character()), length(vep_cols)), vep_cols))
}

# Lift key=value annotations out of the Extra field into dedicated columns so
# downstream consumers (e.g. projects/proteogenomic_database) get clean gene symbols. SYMBOL is
# the gene the variant sits in; NEAREST is the nearest gene (for intergenic hits).
get_extra <- function(x, key) {
  val <- sub(sprintf(".*(?:^|;)%s=([^;]*).*", key), "\\1", x, perl = TRUE)
  val[val == x] <- NA_character_   # key absent -> sub leaves string unchanged
  val
}
if ("Extra" %in% names(vep_dt) && nrow(vep_dt)) {
  vep_dt[, SYMBOL := get_extra(Extra, "SYMBOL")]
  vep_dt[, NEAREST := get_extra(Extra, "NEAREST")]
}

# Fill is_cis for clump-lead-only rows (NA, i.e. not from finemap). The clump-lead cis
# rule lives in clump/06_collect_pqtl.R; take its flag rather than a second copy.
if (use_protein_info && any(is.na(snp_dt$is_cis))) {
  pq <- fread(file.path(wb_info_dir, "bf_protein/clump_lead_pqtl.csv.gz"),
              select = c("cohort", "protein_assay", "SNP", "is_cis"))
  pq <- pq[, .(trait = paste0(cohort, "__", protein_assay), SNP, is_cis_pq = is_cis)]
  if (anyDuplicated(pq, by = c("trait", "SNP"))) stop("clump_lead_pqtl.csv.gz has duplicate (trait, SNP)")
  need <- snp_dt[is.na(is_cis), .(trait, SNP)]
  miss <- need[!pq, on = c("trait", "SNP")]
  if (nrow(miss)) stop(nrow(miss), " clump leads are not in clump_lead_pqtl.csv.gz, e.g. ",
                       miss$trait[1], " ", miss$SNP[1], " -- rerun clump/06_collect_pqtl.R")
  snp_dt[pq, on = c("trait", "SNP"), is_cis := fifelse(is.na(is_cis), i.is_cis_pq, is_cis)]
}

# Per-trait counts: credible sets, finemap SNPs, and lead cis/trans.
sig_counts <- snp_dt[, .(n_cs          = uniqueN(cs_id[!is.na(cs_id)]),
                         n_finemap_snp = uniqueN(SNP[!is.na(cs_id)])), by = trait]
if (use_protein_info) {
  lead_cis <- snp_dt[is_clump_lead == TRUE,
                     .(n_lead_cis   = sum(is_cis %in% TRUE),
                       n_lead_trans = sum(is_cis %in% FALSE)), by = trait]
  sig_counts <- merge(sig_counts, lead_cis, by = "trait", all = TRUE)
  message("Lead cis: ", sum(lead_cis$n_lead_cis), "  trans: ", sum(lead_cis$n_lead_trans))
}
by_col  <- if (use_protein_info) "protein_assay" else "trait"
summary <- merge(summary, sig_counts, by.x = by_col, by.y = "trait", all.x = TRUE)
cnt_cols <- intersect(c("n_cs", "n_finemap_snp", "n_lead_cis", "n_lead_trans"), names(summary))
summary[, (cnt_cols) := lapply(.SD, function(x) { x[is.na(x)] <- 0L; x }), .SDcols = cnt_cols]

if (use_protein_info) {
  # No LDSC h2 for pQTL. Measured on the previous run (3,915 bf1 protein GWAS):
  # median h2_se 0.867 against a median h2 of 0.275, 38% negative h2, 634 plasma
  # assays above 1, only 6% with h2/SE > 2, and median mean_chi2 1.001-1.007.
  # LDSC needs a polygenic background and N in the thousands; a pQTL puts nearly
  # all of its variance in one cis locus and these GWAS have per-variant N of
  # 525-792. The estimates were fitting noise. The calibration diagnostic worth
  # keeping (lambda_gc) is already computed per assay by clean_qc, so take it,
  # and the sample-size columns, from there instead.
  # Only the bf1 pQTL clean path passes QC_DIR, so bf2 has no protein_qc yet.
  # Skip the merge rather than failing: these are diagnostic columns, and the
  # h2 merge they replaced degraded the same way.
  qc_dir <- file.path(gwas_dir, "clean", paste0(type_dir, "_qc"))
  qc_files <- list.files(qc_dir, pattern = "\\.qc\\.tsv$", full.names = TRUE)
  if (!length(qc_files)) {
    message("No per-assay QC files in ", qc_dir,
            " -- summary will carry no lambda_gc/n_variants/n_median columns")
  } else {
    # Key on the file name (<cohort>__<assay>.qc.tsv), as summary is keyed: the trait column
    # inside the file carries no cohort prefix, so joining on it matched nothing (r014).
    qc_dt <- rbindlist(lapply(qc_files, function(f)
                         fread(f, colClasses = c(trait = "character"))[, trait := sub("\\.qc\\.tsv$", "", basename(f))]),
                       fill = TRUE)[, .(trait, n_variants, n_median, lambda_gc)]
    if (anyDuplicated(qc_dt$trait)) stop("duplicate traits in ", qc_dir)
    summary <- merge(summary, qc_dt, by.x = by_col, by.y = "trait", all.x = TRUE)
    message("Merged QC diagnostics for ", sum(!is.na(summary$lambda_gc)), "/",
            nrow(summary), " traits from ", length(qc_files), " QC files")
  }
} else if (type_dir == "pheno") {
  # Phenotypes are polygenic, so LDSC h2 stays meaningful for them. 05_ldsc_h2.sh
  # covers every clean pheno, a superset of the clumped `traits` annotated here.
  h2_file <- file.path(wb_info_dir, "h2", sprintf("%s_h2_%s.tsv", GWAS, type_dir))
  if (!file.exists(h2_file)) stop("No LDSC h2 table at ", h2_file,
                                  " -- finish 05_ldsc_h2.sh before annotation")
  h2_dt <- fread(h2_file, colClasses = c(trait = "character"))
  if (anyDuplicated(h2_dt$trait)) stop("Duplicate traits in ", h2_file)
  missing_h2 <- setdiff(traits, h2_dt$trait)
  if (length(missing_h2))
    stop(length(missing_h2), " annotated traits are missing from ", h2_file,
         ", e.g. ", missing_h2[1], " -- finish and recollect 05_ldsc_h2.sh")
  clean_trait_dir <- file.path(gwas_dir, "clean", type_dir)
  clean_files <- file.path(clean_trait_dir, paste0(traits, ".tsv.gz"))
  if (any(!file.exists(clean_files)))
    stop(sum(!file.exists(clean_files)), " clean GWAS files are missing for annotation traits")
  if (any(file.mtime(clean_files) > file.mtime(h2_file)))
    stop(sum(file.mtime(clean_files) > file.mtime(h2_file)),
         " clean GWAS files are newer than ", h2_file,
         " -- rerun LDSC before annotation")
  h2_dt <- h2_dt[, .(trait, lambda_gc, h2, h2_se)]
  summary <- merge(summary, h2_dt, by.x = by_col, by.y = "trait", all.x = TRUE)
} else {
  # finngen / other: LDSC is run on BF phenotypes only.
  message("No LDSC h2 for ", GWAS, " -- summary will carry no h2 columns")
}

fwrite(summary, file.path(anno_dir, "summary.tsv"), sep = "\t")
message("Wrote summary.tsv (", nrow(summary), " traits)")

# vep_dt may have multiple consequence rows per SNP; keep one per SNP to avoid
# row explosion when joining with snp_dt (many trait x SNP pairs).
vep_dt_uniq <- unique(vep_dt, by = "Uploaded_variation")
message("VEP rows: ", nrow(vep_dt), " -> unique SNPs: ", nrow(vep_dt_uniq))
# Left join that keeps snp_dt's row and column order: one VEP row per SNP (unique
# above), so match() cannot multiply rows.
snp_anno <- cbind(snp_dt, vep_dt_uniq[match(snp_dt$SNP, vep_dt_uniq$Uploaded_variation),
                                      !"Uploaded_variation"])
stopifnot(nrow(snp_anno) == nrow(snp_dt))
# Every variant needs a VEP row. A VEP run killed part-way left NA consequences in a
# table that was written anyway, and the sbatch never reran VEP once it existed.
setDT(snp_anno)
n_unannot <- uniqueN(snp_anno[is.na(Location), SNP])
if (n_unannot) stop(n_unannot, " variants have no VEP row; delete ", vep_out, " and rerun VEP")
fwrite(snp_anno, file.path(anno_dir, "variant_annotation.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
message("Wrote variant_annotation.tsv.gz (", nrow(snp_anno), " rows)")
file.remove(c(out_vcf, vep_out))
message("Deleted temporary VEP files.")
message("07_annotation.R Section 4 done.")
