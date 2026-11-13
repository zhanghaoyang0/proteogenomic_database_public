# Build BF1/BF2 genetic sample ID mapping: sid -> FID, IID
# Output:
#   ~/db/secure_data/bf/pheno/clean/bf1_sid_iid.csv
#   ~/db/secure_data/bf/pheno/clean/bf2_sid_iid.csv
# Usage: Rscript 03_make_genetic_key.R
#
# A row is kept only if its IID is a sample in the RAW array delivery (before
# imputation): that is what "has array genotype" means. The sample xlsx also
# carries placeholders in the IID column for samples without usable genotypes
# (BF1: "QC FAIL" on 134 sids, 2026-10-04); written through, they made one fake
# "IID" shared by many sids. Dropped rows are counted, and placeholder strings
# (any IID value shared by several sids) are printed -- they are not sample ids.
#
# The kept IIDs must also be the ids of the imputed pfiles that downstream
# joins against (01_vcf_import renames via the genetic key); the overlap is
# checked and a low one stops the script rather than writing a key no
# pipeline can use.

library(data.table)
library(readxl)

genome_dir <- path.expand(Sys.getenv("BF_ROOT", "~/db/secure_data/bf"))
pheno_clean_dir <- file.path(genome_dir, "pheno/clean")
dir.create(pheno_clean_dir, showWarnings = FALSE, recursive = TRUE)
MIN_IMPUTED_OVERLAP <- 0.9   # share of kept IIDs that must be in the imputed psam

norm_id <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "<NA>")] <- NA_character_
  x
}

pick_first <- function(nm, rx) {
  hit <- grep(rx, nm, value = TRUE, ignore.case = TRUE)[1]
  if (is.na(hit)) stop("xlsx: no column matching ", rx)
  hit
}

# IIDs are compared without leading zeros: the BF2 raw .fam writes them with an
# extra leading 0 (8 digits) that the xlsx and the imputed pfiles do not carry
# (7 digits) -- same people, same numbers (checked 2026-10-04: all 1817 xlsx IIDs
# match the .fam once zeros are stripped). 06_qc_harmonize_1kg.sh compares the
# same way. The key keeps the xlsx spelling, which is what the pfiles use.
zstrip <- function(x) sub("^0+(?=.)", "", x, perl = TRUE)

# sid comparison key: case and -/_/space ignored. The BF1 xlsx also writes 11 MLDS
# sids as MLDSMIX### where the database has MLDSMI### (same number; the 10 of them
# in the 2024 key array/bf1/archive/BF1_genetic_key.csv have the database spelling;
# its other pairs all recur here except 2 Parkinson sids no longer in the database,
# 2026-11-03), and no
# database sid starts MLDSMIX, so that prefix is read as MLDSMI.
skey <- function(x) sub("^MLDSMIX", "MLDSMI", toupper(gsub("[-_[:space:]]", "", x)))

# Phenotype-database spelling of every sid, looked up by skey. A skey with several
# spellings in the table is ambiguous; build_key stops if a key sid lands on one.
# Read from the id columns of the PET export (code/01_genotype_qc/02_extract_sample_ids.py),
# not the 1.5 GB export itself; the PET export is a superset of the main export
# (checked 2026-10-09) and what process_pheno/02 and code/02_protein_qc/02 read.
pheno_main <- file.path(genome_dir, "pheno/raw/haoyang__20260227_112046_PET/haoyang__20260227_112046_x.csv")
pheno_ids <- file.path(genome_dir, "pheno/clean/haoyang__20260227_112046_PET_ids.csv")
for (f in c(pheno_main, pheno_ids)) if (!file.exists(f)) stop("missing ", f)
if (file.mtime(pheno_main) > file.mtime(pheno_ids)) {
  stop(pheno_ids, " is older than the PET export -- rerun code/01_genotype_qc/02_extract_sample_ids.py")
}
# Sample-level ids, one row per (sid, Visit): a key sid that is not a person sid
# may be one of these (BF1, 2026-10-13: array ids that are an OBARCODE or a
# visit's blood_sample_id, found by the WGS-array identity check), so a key sid
# not found among the sids is looked up here and replaced by the row's sid.
# csf_sample_id is left out: no key sid is one, and its values are reused across
# sids (up to 15 sids per value, 00_extract_ids.log). blood_sample_id is reused
# too (up to 18); a key sid landing on a reused value stops the script below.
SAMPLE_ID_COLS <- c("OBARCODE", "blood_sample_id")
pet <- fread(pheno_ids, select = c("sid", SAMPLE_ID_COLS), colClasses = "character")
if (!all(c("sid", SAMPLE_ID_COLS) %in% names(pet))) {
  stop("id table lacks: ", paste(setdiff(c("sid", SAMPLE_ID_COLS), names(pet)), collapse = ", "))
}
for (j in names(pet)) set(pet, j = j, value = norm_id(pet[[j]]))
pheno_sids <- unique(na.omit(pet$sid))
pheno_ambiguous <- unique(skey(pheno_sids)[duplicated(skey(pheno_sids))])
pheno_sid <- setNames(pheno_sids, skey(pheno_sids))[!duplicated(skey(pheno_sids))]
message("phenotype table: ", length(pheno_sids), " sids (", length(pheno_ambiguous), " ambiguous skeys)")

sample_alias <- unique(rbindlist(lapply(SAMPLE_ID_COLS, function(j)
  pet[!is.na(get(j)) & !is.na(sid), .(skey = skey(get(j)), sid, col = j)])), by = c("skey", "sid", "col"))
alias_ambiguous <- sample_alias[, .(n = uniqueN(sid)), by = skey][n > 1, skey]
alias_first <- sample_alias[!duplicated(skey)]
alias_sid <- setNames(alias_first$sid, alias_first$skey)
alias_col <- setNames(alias_first$col, alias_first$skey)
message("phenotype table sample ids: ", length(alias_sid), " (", length(alias_ambiguous),
        " on more than one sid)")
rm(pet)

# digit-masked format shape, safe to print
id_shape <- function(x) gsub("[A-Za-z]", "A", gsub("[0-9]", "#", x))

# IIDs of a plink sample file: .psam (header "#IID" or "#FID IID") or .fam (no header, IID = col 2)
read_iids <- function(path) {
  if (!file.exists(path)) stop("missing ", path)
  if (grepl("\\.fam$", path)) {
    iid <- fread(path, header = FALSE, select = 2, colClasses = "character")[[1]]
  } else {
    d <- fread(path, colClasses = "character")
    col <- intersect(c("IID", "#IID"), names(d))[1]
    if (is.na(col)) stop("no IID or #IID column in ", path)
    iid <- d[[col]]
  }
  unique(norm_id(iid))
}

build_key <- function(label, xlsx, sid_rx, raw_samples, imputed_psam, out) {
  message("Processing ", label, "...")
  genotyped <- read_iids(raw_samples)
  message("  samples in raw array delivery: ", length(genotyped))
  geno_z <- zstrip(genotyped)

  # Both xlsx have duplicated column names (two FID/IID blocks; unique_quiet
  # renames them quietly). Use the IID column holding the most genotyped
  # samples, with the FID of the same block: in BF2 the lab overwrote the first
  # IID column with "Removed in QC" for 127 samples (2026-10-09, see
  # ../docs/issues.md) while the second keeps every IID. QC is not decided here.
  raw <- as.data.table(read_excel(xlsx, col_types = "text", .name_repair = "unique_quiet"))
  iid_cols <- grep("^IID", names(raw), value = TRUE)
  fid_cols <- grep("^FID", names(raw), value = TRUE)
  if (!length(iid_cols)) stop(label, ": xlsx has no IID column")
  hits <- sapply(iid_cols, function(j) sum(zstrip(norm_id(raw[[j]])) %in% geno_z))
  k <- which.max(hits)
  message("  genotyped IIDs per xlsx IID column: ",
          paste(sprintf("%s %d", iid_cols, hits), collapse = ", "), " -> using ", iid_cols[k])
  if (length(fid_cols) != length(iid_cols)) stop(label, ": ", length(fid_cols), " FID vs ",
                                                 length(iid_cols), " IID columns")
  cols <- c(fid_cols[k], iid_cols[k], pick_first(names(raw), sid_rx))
  raw <- raw[, ..cols]
  setnames(raw, c("FID", "IID", "sid"))
  for (j in names(raw)) set(raw, j = j, value = norm_id(raw[[j]]))
  raw <- raw[!is.na(IID) & !is.na(sid)]
  message("  xlsx rows with sid and IID: ", nrow(raw))

  # Placeholders: one value on several sids AND not a genotyped sample. Only
  # these are printed; a real IID shared by sids is caught (count only) below.
  shared <- raw[, .(n_sid = uniqueN(sid)), by = IID][n_sid > 1 & !(zstrip(IID) %in% geno_z)]
  if (nrow(shared)) {
    message("  placeholder IID values (several sids, not in raw delivery; printed as-is):")
    for (i in seq_len(nrow(shared))) message(sprintf("    %-20s on %d sids", shared$IID[i], shared$n_sid[i]))
  }
  key <- raw[zstrip(IID) %in% geno_z]
  if (!nrow(key)) stop(label, ": no xlsx IID is in ", raw_samples, " -- id format mismatch; not writing ", out)
  message("  dropped, IID not in raw delivery: ", nrow(raw) - nrow(key),
          " (of which placeholders: ", sum(raw$IID %in% shared$IID), ")")

  # The xlsx spells some sids differently from the phenotype database, and every
  # downstream join on sid is exact (BF1, 2026-10-09: 313 differ in case, 281 lack
  # the dash, so 594 genotyped people lost age/sex and dropped out of the GWAS).
  # Write the database spelling when the sid matches it ignoring case and -/_/space.
  k_xlsx <- skey(key$sid)
  hit <- k_xlsx %in% names(pheno_sid)
  if (any(hit & k_xlsx %in% pheno_ambiguous)) {
    stop(label, ": ", sum(hit & k_xlsx %in% pheno_ambiguous),
         " sids match several phenotype-table spellings")
  }
  respelled <- hit & key$sid != pheno_sid[k_xlsx]
  key[hit, sid := pheno_sid[k_xlsx[hit]]]
  message("  sids in phenotype table: ", sum(hit), " / ", nrow(key),
          " (respelled to its spelling: ", sum(respelled), "; not in it: ", sum(!hit), ")")

  # Not a person sid: try the sample-level ids (see SAMPLE_ID_COLS above).
  via <- !hit & k_xlsx %in% names(alias_sid)
  if (any(via & k_xlsx %in% alias_ambiguous)) {
    stop(label, ": ", sum(via & k_xlsx %in% alias_ambiguous),
         " key sids are a sample id of more than one phenotype-table sid")
  }
  if (any(via)) {
    by_col <- table(alias_col[k_xlsx[via]])
    message("  resolved through a sample id to the person's sid: ", sum(via), " (",
            paste(sprintf("%s %d", names(by_col), by_col), collapse = ", "), ")")
    key[via, sid := alias_sid[k_xlsx[via]]]
  }
  left <- !hit & !via
  if (any(left)) {
    sh <- sort(table(id_shape(key$sid[left])), decreasing = TRUE)
    message("  still not in the phenotype table: ", sum(left), " -- shapes: ",
            paste(sprintf("%s %d", names(sh), sh), collapse = ", "))
  }
  # Their IIDs leave in Stage 2 with the sex-check exclusions: they have no phenotype,
  # and left in they only take part in het/KING, where --king-cutoff dropped the
  # phenotyped member of 3 BF1 pairs instead of them (2026-11-03; BF1: 75 are TauPET
  # sids, a sub-study outside the export). Written even when empty, next to the key
  # and only after its checks pass, so the two files always come from the same run.
  no_pheno_iid <- unique(key$IID[left])

  key <- unique(key, by = c("sid", "IID"))
  if (anyDuplicated(key$sid)) stop(label, ": ", sum(duplicated(key$sid)), " sids map to more than one IID",
                                   if (any(via)) " -- possibly a person genotyped twice, found by the sample-id lookup" else "")
  if (anyDuplicated(zstrip(key$IID))) stop(label, ": ", sum(duplicated(zstrip(key$IID))), " IIDs map to more than one sid")

  imputed <- read_iids(imputed_psam)
  in_imp <- zstrip(key$IID) %in% zstrip(imputed)
  overlap <- mean(in_imp)
  message(sprintf("  kept IIDs also in imputed psam: %d / %d (%.1f%%)",
                  sum(in_imp), nrow(key), 100 * overlap))
  if (overlap < MIN_IMPUTED_OVERLAP) {
    stop(label, ": only ", round(100 * overlap, 1), "% of kept IIDs are in ", imputed_psam,
         " -- raw and imputed IIDs are not the same id space; not writing ", out)
  }

  no_pheno_out <- sub("_sid_iid\\.csv$", "_iid_no_pheno.txt", out)
  writeLines(no_pheno_iid, no_pheno_out)
  message("  IIDs without a phenotype-table sid: ", length(no_pheno_iid), " -> ", no_pheno_out)
  fwrite(key[, .(sid, FID, IID)], out)
  message("  ", label, " key rows: ", nrow(key), " -> ", out)
  invisible(key)
}

# BF1 (since 2026-11-02): the raw 1680-sample bfile, whose QC FAIL samples carry
# their ids only in the xlsx's second IID column; the overlap is checked against
# the pre-impute psam, whose IIDs are the ones sent to TOPMed and exists before
# Stage 2 needs this key.
build_key("BF1",
          xlsx         = file.path(genome_dir, "array/bf1/raw/Biofinder_Main_GWAS_DATA.xlsx"),
          sid_rx       = "Clinical_ID|^sid",
          raw_samples  = file.path(genome_dir, "array/bf1/raw/BioFINDER1RAWDATA/GSP_2018-054-biofinder.fam"),
          imputed_psam = file.path(genome_dir, "array/bf1/for_impute/bf1_hg19.psam"),
          out          = file.path(pheno_clean_dir, "bf1_sid_iid.csv"))

build_key("BF2",
          xlsx         = file.path(genome_dir, "array/bf2/raw/GWAS_Full_Data_BF-2.xlsx"),
          sid_rx       = "^sid",
          raw_samples  = file.path(genome_dir, "array/bf2/raw/UB-2844_240229_PLINK_PCF_PLUS/1_BF-2_Genetic_Data.fam"),
          imputed_psam = file.path(genome_dir, "array/bf2/imputed_qc/bf2_hg19_minimalqc.psam"),
          out          = file.path(pheno_clean_dir, "bf2_sid_iid.csv"))

message("\n✓ Complete: bf1_sid_iid.csv and bf2_sid_iid.csv generated")
