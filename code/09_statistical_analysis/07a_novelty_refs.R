#!/usr/bin/env Rscript
# The pQTL novelty scan (docs/pipeline.md):
# turn the published supplementary workbooks into tidy reference tables.
#
#   Cruchaga CSF (Nat Genet 2024, PMID 39528825)  - SomaScan 7k, CSF, N=3,506
#   Western CSF  (Sci Transl Med 2026, PMID 42054495) - Olink Explore HT, CSF, N=2,524
#   UKB-PPP      (Nature 2023,   PMID 37794186)   - Olink 3k,    plasma, N=54,219
#
# Two products per reference:
#   *_panel.csv      - every protein the reference assayed. Defines tier C (a protein
#                      they never measured is coverage, not novelty).
#   *_sentinels.csv  - every association they report. Defines tier A (no matching
#                      sentinel = previously unreported).
#
# Build: UKB-PPP variant IDs already carry hg19 positions, so plasma needs no
# conversion. Cruchaga is GRCh38 only, so its sentinel *coordinates* are lifted to
# hg19 here and the lift is verified two ways (allele match in the 1000G EUR panel,
# and rsID agreement with their own Supplementary Table 2). Only sentinel positions
# cross the build boundary, never summary statistics.
#
# Usage: Rscript code/09_statistical_analysis/07a_novelty_refs.R

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

REF_DIR   <- path.expand("~/db/shared_data/pqtl_ref")
CR_XLSX   <- file.path(REF_DIR, "cruchaga_csf", "41588_2024_1972_MOESM4_ESM.xlsx")
UK_XLSX   <- file.path(REF_DIR, "ukbppp", "41586_2023_6592_MOESM3_ESM.xlsx")
CR_URL    <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-024-01972-8/MediaObjects/41588_2024_1972_MOESM4_ESM.xlsx"
UK_URL    <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06592-6/MediaObjects/41586_2023_6592_MOESM3_ESM.xlsx"
BIM       <- path.expand("~/db/shared_data/ref_geno/1kg_v3/EUR.bim")
CHAIN     <- path.expand("~/db/shared_data/ref_geno/liftover/hg38ToHg19.over.chain")
WO_DIR    <- file.path(REF_DIR, "western_csf_olink")
WO_S4     <- file.path(WO_DIR, "table_s4_olink_analytes.csv")
WO_S5     <- file.path(WO_DIR, "table_s5_olink_csf_index_pqtls.csv")
WO_GCST   <- file.path(WO_DIR, "gwas_catalog_studies.csv")
GC_API    <- "https://www.ebi.ac.uk/gwas/rest/api/v2/studies?pubmed_id=42054495&size=200&page="
OUT_DIR   <- file.path(REF_DIR, "parsed")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

fetch <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 1e6) return(invisible(dest))
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  message("downloading ", basename(dest))
  if (download.file(url, dest, mode = "wb", quiet = TRUE) != 0 || file.size(dest) < 1e6)
    stop("download failed: ", url)
  invisible(dest)
}
fetch(CR_URL, CR_XLSX)
fetch(UK_URL, UK_XLSX)
if (!file.exists(BIM)) stop("missing LD panel bim: ", BIM)

read_sheet <- function(xlsx, sheet, skip) {
  as.data.table(suppressMessages(
    read_excel(xlsx, sheet = sheet, skip = skip, .name_repair = "minimal")))
}

# ---- UKB-PPP plasma ---------------------------------------------------------
# ST3 = the Olink Explore 3k panel; ST9 = the 14,287 significant discovery pQTLs
# (p < 1.7e-11), with the hg19 position embedded in the variant ID.
uk_panel <- read_sheet(UK_XLSX, "ST3", skip = 1)
setnames(uk_panel,
         c("UniProt", "Gene symbol", "Olink ID", "Gene CHROM", "Gene start", "Gene end"),
         c("uniprot_id", "gene", "assay_id", "gene_chr", "gene_start", "gene_end"))
uk_panel <- unique(uk_panel[, .(uniprot_id, gene, assay_id, gene_chr, gene_start, gene_end)])

uk_sent <- read_sheet(UK_XLSX, "ST9", skip = 3)
setnames(uk_sent,
         c("Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)", "Target UniProt", "rsID",
           "cis/trans", "log10(p) (discovery)", "BETA (discovery, wrt. A1)"),
         c("variant_id", "uniprot_id", "rsid", "cis_trans", "log10p", "beta"))
# 167 associations sit on X, written as "X:" in the variant ID and 23 in CHROM;
# keep them as 23 so nothing is silently dropped (our data are autosomal, so they
# simply never match).
uk_sent[, chr := as.integer(sub("^X$", "23", tstrsplit(variant_id, ":", keep = 1)[[1]]))]
uk_sent[, pos_hg19 := as.integer(tstrsplit(variant_id, ":", keep = 2)[[1]])]
uk_sent <- uk_sent[, .(uniprot_id, rsid, chr, pos_hg19,
                       cis_trans = tolower(cis_trans), log10p, beta,
                       source = "ST9 primary")]
if (nrow(uk_sent) != 14287)
  stop("UKB-PPP ST9 should hold 14,287 associations, got ", nrow(uk_sent))
if (anyNA(uk_sent$pos_hg19)) stop("unparsed hg19 position in UKB-PPP ST9")

# ST9 keeps only the most significant association per locus, so on its own it makes any
# secondary signal of ours look unreported. ST16 carries their 29,420 independent
# (SuSiE) signals; both are used, so "unreported" means unreported in either.
uk_ind <- read_sheet(UK_XLSX, "ST16", skip = 3)
setnames(uk_ind, seq_len(12),
         c("protein_id", "variant_id", "rsid", "pip", "log10bf", "beta", "se",
           "log10p", "alt_freq", "cis_trans", "test_region", "ld_other_signals"))
uk_ind <- uk_ind[variant_id != "Variant ID"]   # ST16 repeats its header once
if (nrow(uk_ind) != 29420)
  stop("UKB-PPP ST16 should hold 29,420 independent signals, got ", nrow(uk_ind))
uk_ind[, uniprot_id := tstrsplit(protein_id, ":", keep = 2)[[1]]]
uk_ind[, chr := as.integer(sub("^X$", "23", tstrsplit(variant_id, ":", keep = 1)[[1]]))]
uk_ind[, pos_hg19 := as.integer(tstrsplit(variant_id, ":", keep = 2)[[1]])]
uk_ind <- uk_ind[, .(uniprot_id, rsid, chr, pos_hg19,
                     cis_trans = tolower(cis_trans), log10p, beta,
                     source = "ST16 independent")]
if (anyNA(uk_ind$pos_hg19)) stop("unparsed hg19 position in UKB-PPP ST16")
if (anyNA(uk_ind$uniprot_id)) stop("unparsed UniProt in UKB-PPP ST16")
message("UKB-PPP sentinels: ", nrow(uk_sent), " primary + ", nrow(uk_ind),
        " independent")
uk_sent <- unique(rbind(uk_sent, uk_ind))

# ---- Cruchaga CSF -----------------------------------------------------------
# ST1 = SomaScan 7k annotation; aptamers with a non-empty "Step Removed" were dropped
# before their GWAS, so the analysed panel is the 7,008 with no removal reason.
# ST8 = the 3,885 conditionally significant associations (the paper's headline count,
# a superset of the ST2 index pQTLs). ST8 has no rsID, ST2 does - used below to audit
# the lift.
cr_annot <- read_sheet(CR_XLSX, "Supplementary Table 1", skip = 1)
setnames(cr_annot, c("Analytes", "UniProt", "EntrezGeneSymbol", "Step Removed"),
         c("aptamer", "uniprot_id", "gene", "step_removed"))
cr_annot <- cr_annot[, .(aptamer, uniprot_id, gene, step_removed)]
cr_panel <- cr_annot[is.na(step_removed)][, step_removed := NULL]
if (nrow(cr_panel) != 7008)
  stop("Cruchaga analysed panel should hold 7,008 aptamers, got ", nrow(cr_panel))

cr_sent <- read_sheet(CR_XLSX, "Supplementary Table 8", skip = 1)
setnames(cr_sent, c("CHROM", "POS", "ID", "REF", "ALT", "SomaID", "pQTL Type",
                    "LOG10_P", "BETA"),
         c("chr", "pos_hg38", "variant_id", "ref", "alt", "aptamer", "cis_trans",
           "log10p", "beta"))
cr_sent <- cr_sent[, .(aptamer, chr = as.integer(chr), pos_hg38 = as.integer(pos_hg38),
                       variant_id, ref, alt, cis_trans = tolower(cis_trans),
                       log10p, beta, source = "ST8 conditional")]
if (nrow(cr_sent) != 3885)
  stop("Cruchaga ST8 should hold 3,885 associations, got ", nrow(cr_sent))
cr_sent <- merge(cr_sent, cr_panel[, .(aptamer, uniprot_id, gene)], by = "aptamer",
                 all.x = TRUE)
message("Cruchaga sentinels with no UniProt from the panel annotation: ",
        cr_sent[is.na(uniprot_id) | uniprot_id == "", .N], " of ", nrow(cr_sent))

# ---- lift GRCh38 sentinels to hg19 ------------------------------------------
# Cruchaga and Western 2026 are GRCh38. lift_hg19() adds pos_hg19 and, from the LD panel
# we will compute r2 in, the rsID at that position with the same alleles: a sentinel
# that is not in the panel cannot be LD-matched anyway, so panel membership is the
# useful definition of "found". audit_lift() then checks the rsIDs against the study's own.
lift_bin <- Sys.which("liftOver")
if (lift_bin == "") lift_bin <- path.expand("~/miniforge3/envs/gwas/bin/liftOver")
if (!file.exists(lift_bin)) stop("liftOver binary not found")
bim <- fread(BIM, header = FALSE, select = c(1, 2, 4, 5, 6),
             col.names = c("chr", "rsid", "pos_hg19", "a1", "a2"))
bim[, `:=`(al1 = pmin(a1, a2), al2 = pmax(a1, a2))]

lift_hg19 <- function(d, label) {   # d: chr, pos_hg38, variant_id, ref, alt
  tmpd <- tempfile("liftover"); dir.create(tmpd)
  bed_in  <- file.path(tmpd, "in.bed")
  bed_out <- file.path(tmpd, "out.bed")
  lift_in <- unique(d[, .(chr, pos_hg38, variant_id)])
  fwrite(lift_in[, .(paste0("chr", chr), pos_hg38 - 1L, pos_hg38, variant_id)],
         bed_in, sep = "\t", col.names = FALSE)
  system2(lift_bin, c(bed_in, CHAIN, bed_out, file.path(tmpd, "unmapped.bed")),
          stdout = NULL, stderr = NULL)
  if (!file.exists(bed_out)) stop("liftOver produced no output")
  lift <- fread(bed_out, header = FALSE,
                col.names = c("lift_chr", "start", "pos_hg19", "variant_id"))
  lift[, lift_chr := sub("^chr", "", lift_chr)]
  d <- merge(d, lift[, .(variant_id, lift_chr, pos_hg19)], by = "variant_id", all.x = TRUE)
  d[lift_chr != as.character(chr), pos_hg19 := NA_integer_]  # cross-chromosome lift
  message(label, " sentinels lifted to hg19: ", d[!is.na(pos_hg19), .N], " of ", nrow(d))
  d[, `:=`(al1 = pmin(ref, alt), al2 = pmax(ref, alt))]
  d <- merge(d, bim[, .(chr, pos_hg19, al1, al2, rsid)],
             by = c("chr", "pos_hg19", "al1", "al2"), all.x = TRUE)
  d[, c("al1", "al2") := NULL]
  message(label, " sentinels matched in the 1000G EUR panel by position+alleles: ",
          d[!is.na(rsid), .N], " of ", nrow(d))
  d
}

# rsid_paper: variant_id -> the rsID the study itself reports.
audit_lift <- function(d, rsid_paper, label) {
  audit <- unique(merge(d[!is.na(rsid), .(variant_id, rsid)], rsid_paper, by = "variant_id"))
  message(label, " lift audit vs their own rsIDs: ", audit[rsid == rsid_paper, .N], " agree / ",
          audit[rsid != rsid_paper, .N], " disagree, of ", nrow(audit), " checkable")
  if (nrow(audit) > 0 && audit[rsid == rsid_paper, .N] / nrow(audit) < 0.98)
    stop("liftOver of the ", label, " sentinels does not reproduce their rsIDs - stop and check")
}

cr_sent <- lift_hg19(cr_sent, "Cruchaga")
# audit against ST2, the index pQTLs with the rsIDs Cruchaga report (ST8 has none)
cr_idx <- read_sheet(CR_XLSX, "Supplementary Table 2", skip = 1)
audit_lift(cr_sent, unique(cr_idx[, .(variant_id = ID, rsid_paper = rsID)])[!is.na(rsid_paper)],
           "Cruchaga")

setcolorder(cr_sent, c("aptamer", "uniprot_id", "gene", "chr", "pos_hg19", "pos_hg38",
                       "rsid", "variant_id", "ref", "alt", "cis_trans", "log10p",
                       "beta", "source"))

# ---- Western CSF Olink ------------------------------------------------------
# Olink Explore HT, the platform family of our CSF, so a CSF protein that SomaScan did
# not show can still be reported here. The supplement is an .xlsb behind the journal's
# paywall, saved by hand to WO_DIR (2026-10-30); readxl has no .xlsb reader, so tables S4
# (5,420 analytes) and S5 (1,062 index pQTLs; no conditional signals) were exported beside
# it to CSV, unchanged, with pandas + pyxlsb.
for (f in c(WO_S4, WO_S5)) if (!file.exists(f)) stop("missing Western 2026 table: ", f)
wo_annot <- fread(WO_S4)
if (nrow(wo_annot) != 5420) stop("Western S4 should hold 5,420 analytes, got ", nrow(wo_annot))

# Analysed panel: S4 flags LOD and CV but does not say which analytes were analysed. The
# GWAS Catalog holds one study per analysed analyte, traits named "Cerebrospinal fluid
# protein <gene> levels", so the panel is the S4 analytes with a study (all but 90, every
# one of which fails LOD). The list is fetched once to WO_GCST; one study has no trait in
# the API and takes it from its metadata file.
if (!file.exists(WO_GCST)) {
  message("fetching the Western 2026 study list from the GWAS Catalog (~10 min)")
  get_page <- function(p) {
    for (i in 1:3) {
      x <- tryCatch(jsonlite::fromJSON(paste0(GC_API, p)), error = function(e) NULL)
      if (!is.null(x$`_embedded`$studies)) return(x)
    }
    stop("GWAS Catalog page ", p, " failed three times")
  }
  first <- get_page(0)
  pages <- c(list(first), lapply(seq_len(first$page$totalPages - 1), get_page))
  gc <- rbindlist(lapply(pages, function(x) {
    s <- x$`_embedded`$studies
    data.table(accession = s$accession_id, ftp_dir = s$full_summary_stats,
               trait = if (is.null(s$disease_trait)) NA_character_ else s$disease_trait)
  }))
  if (nrow(gc) != first$page$totalElements || anyDuplicated(gc$accession))
    stop("GWAS Catalog list incomplete: ", nrow(gc), " of ", first$page$totalElements)
  for (i in gc[is.na(trait), which = TRUE]) {
    y <- readLines(sprintf("%s/%s.tsv.gz-meta.yaml", gc$ftp_dir[i], gc$accession[i]))
    set(gc, i, "trait", trimws(sub("^\\s*-", "", y[grep("^trait_description:", y) + 1L])))
  }
  gc[, gene := sub("^Cerebrospinal fluid protein (.*) levels$", "\\1", trait)]
  if (gc[is.na(gene) | gene == trait, .N]) stop("unparsed GWAS Catalog trait")
  fwrite(gc, WO_GCST)
}
gc <- fread(WO_GCST)
if (nrow(gc) != 5330) stop("expected 5,330 Western 2026 GWAS Catalog studies, got ", nrow(gc))
if (!all(gc$gene %in% wo_annot$Assay))
  stop("GWAS Catalog studies with no S4 analyte: ", paste(setdiff(gc$gene, wo_annot$Assay), collapse = ", "))
# GBP1 and MAP2K1 are three analytes each, with three studies each, so matching on the
# name keeps all of them.
wo_panel <- wo_annot[Assay %in% gc$gene,
                     .(olink_id = OlinkID, uniprot_id = gsub("_", "|", UniProt), gene = Assay)]
if (nrow(wo_panel) != nrow(gc)) stop("Western analysed panel: ", nrow(wo_panel), " analytes for ", nrow(gc), " studies")

# Sentinels. type is per gene of the assay ("cis,trans" for a two-gene complex); 07b reads
# any "cis" as a cis-pQTL of the protein.
wo_sent <- fread(WO_S5)
if (nrow(wo_sent) != 1062) stop("Western S5 should hold 1,062 index pQTLs, got ", nrow(wo_sent))
setnames(wo_sent, c("CHROM", "POS", "ID", "REF", "ALT", "rsID", "OlinkID", "type", "LOG10_P", "BETA"),
         c("chr", "pos_hg38", "variant_id", "ref", "alt", "rsid_paper", "olink_id", "cis_trans",
           "log10p", "beta"))
wo_sent <- wo_sent[, .(olink_id, chr = as.integer(chr), pos_hg38 = as.integer(pos_hg38), variant_id,
                       ref, alt, rsid_paper, cis_trans, log10p, beta, source = "S5 index")]
n_s5 <- nrow(wo_sent)
wo_sent <- merge(wo_sent, wo_panel[, .(olink_id, uniprot_id, gene)], by = "olink_id")
if (nrow(wo_sent) != n_s5) stop("Western S5 pQTLs of analytes outside the analysed panel: ", n_s5 - nrow(wo_sent))
wo_sent <- lift_hg19(wo_sent, "Western Olink")
audit_lift(wo_sent, unique(wo_sent[grepl("^rs", rsid_paper), .(variant_id, rsid_paper)]), "Western Olink")
wo_sent[, rsid_paper := NULL]
setcolorder(wo_sent, c("olink_id", "uniprot_id", "gene", "chr", "pos_hg19", "pos_hg38",
                       "rsid", "variant_id", "ref", "alt", "cis_trans", "log10p",
                       "beta", "source"))

# ---- write ------------------------------------------------------------------
fwrite(uk_panel, file.path(OUT_DIR, "ukbppp_panel.csv"))
fwrite(uk_sent,  file.path(OUT_DIR, "ukbppp_sentinels.csv"))
fwrite(cr_panel, file.path(OUT_DIR, "cruchaga_panel.csv"))
fwrite(cr_sent,  file.path(OUT_DIR, "cruchaga_sentinels.csv"))
fwrite(wo_panel, file.path(OUT_DIR, "western_olink_panel.csv"))
fwrite(wo_sent,  file.path(OUT_DIR, "western_olink_sentinels.csv"))

for (f in c("ukbppp_panel.csv", "ukbppp_sentinels.csv",
            "cruchaga_panel.csv", "cruchaga_sentinels.csv",
            "western_olink_panel.csv", "western_olink_sentinels.csv")) {
  p <- file.path(OUT_DIR, f)
  if (!file.exists(p)) stop("expected output missing: ", p)
  message(sprintf("wrote %-26s %6d rows", f, nrow(fread(p, nrows = Inf))))
}
message("07a done")
