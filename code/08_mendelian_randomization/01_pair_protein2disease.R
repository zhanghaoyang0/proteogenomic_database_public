#!/usr/bin/env Rscript

# Generate protein -> disease MR pairs.
# Exposure: bf2 proteins with >=1 cis clump lead (instruments = clump cis-pQTL).
# Outcome: FinnGen + Other (disease GWAS).
# Output: mr_pair.csv.gz with full pair list, mr_pair_todo.csv.gz with pending pairs.

library(data.table)

SCRIPT_DIR <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)[1L]
  if (is.na(fa)) getwd() else dirname(normalizePath(sub("^--file=", "", fa)))
})()
source(file.path(SCRIPT_DIR, "functions_db.R"))
source(file.path(SCRIPT_DIR, "..", "lib", "mhc.R"))


INFO_DIR     <- path.expand("./info")
MR_DIR       <- file.path(MR_ROOT, "protein2disease")
MR_INFO_DIR  <- file.path(MR_DIR, "info")
TMP_DIR      <- file.path(MR_DIR, "tmp_chunks")
DB_FILE      <- file.path(MR_DIR, "mr.db")
RESULT       <- file.path(MR_DIR, "result", "mr_protein2disease.tsv.gz")  # gzipped; fread reads transparently
RESULT_OTH   <- file.path(MR_DIR, "result", "mr_protein2disease_other.tsv")
BF2_CLUMP    <- file.path(GWAS_BASE, "bf2/clump/protein", CLUMP_SUBDIR)
FINNGEN_CLUMP<- file.path(GWAS_BASE, "finngen/clump", CLUMP_SUBDIR)
OTHER_CLUMP  <- file.path(GWAS_BASE, "other/clump", CLUMP_SUBDIR)
assert_clump_dirs(BF2_CLUMP, FINNGEN_CLUMP, OTHER_CLUMP)

dir.create(MR_INFO_DIR, showWarnings = FALSE, recursive = TRUE)

t0  <- proc.time()[[3L]]

# ---- load trait_info ----
ti <- fread(file.path(INFO_DIR, "trait_info.csv"))
msg("trait_info: ", nrow(ti), " rows")

# Per-assay cis windows, one per subunit (gene +/-500kb); windows of a gene in the MHC
# are dropped (../lib/mhc.R).
source(file.path(REPO_ROOT, "code/lib/cis_windows.R"))
cis_win <- read_cis_windows(REPO_ROOT)
cis_win <- cis_win[!mhc_cis_window(chr, cis_start, cis_end)]

# ---- helper: protein clump file has ≥1 lead SNP inside the cis window ----
# Clump-based cis instruments (follow Eldjarn et al.): only proteins with a cis
# clump lead are valid MR exposures (cis-only avoids trans pleiotropy). bf2 only
# (bf1 protein dropped).
clump_has_cis_lead <- function(dir) {
  files <- list.files(dir, pattern = "\\.clumped$", full.names = TRUE)
  ids   <- sub("\\.clumped$", "", basename(files))
  keep <- vapply(seq_along(files), function(j) {
    if (file.size(files[j]) == 0) return(FALSE)
    # Stop on a read error: FALSE would drop the protein from the pairs, and the collect
    # would then delete its results.
    cl <- tryCatch(fread(files[j], select = c("CHR", "BP")),
                   error = function(e) stop("cannot read ", files[j], ": ", conditionMessage(e)))
    if (!nrow(cl)) return(FALSE)
    pa <- sub("^bf[12]__", "", ids[j])
    # No MR in the MHC: neither a window whose gene lies there (dropped above) nor an
    # instrument inside it (../lib/mhc.R). 02b applies the same rule to the instruments.
    any(!mhc_variant(cl$CHR, cl$BP) & in_cis_window(rep(pa, nrow(cl)), cl$CHR, cl$BP, cis_win))
  }, logical(1L))
  ids[keep]
}
clump_cis_ids <- clump_has_cis_lead(BF2_CLUMP)
msg("protein with ≥1 cis clump lead (bf2): ", length(clump_cis_ids))

# ---- outcomes: FinnGen + Other with ≥1 sig SNP (clump lead) ----
disease_sig_ids <- c(clump_has_lead(FINNGEN_CLUMP), clump_has_lead(OTHER_CLUMP))
msg("diseases with ≥1 sig SNP (clump lead): ", length(disease_sig_ids))

# ---- exposures: bf2 protein with ≥1 cis clump lead ----
is_bf_protein <- ti$source %in% c("BioFINDER1", "BioFINDER2") &
                 ti$category %in% c("Protein CSF", "Protein plasma")
expo <- ti[is_bf_protein & trait_id %in% clump_cis_ids]
# trait_id is unique in trait_info (trait_info/collect_trait_info.R stops otherwise).
msg("exposures (bf2 protein with cis clump lead): ", nrow(expo))

# ---- outcomes: FinnGen + Other with ≥1 sig SNP ----
is_disease <- ti$source %in% c("FinnGen meta-analysis", "Other")
outc <- ti[is_disease & trait_id %in% disease_sig_ids]
msg("outcomes (FinnGen + Other with sig SNP): ", nrow(outc))

# ---- cross-join ----
expo_sub <- expo[, .(
  exposure_source      = source,
  exposure_id          = trait_id,
  exposure_name        = trait_name,
  exposure_category    = category,
  exposure_n_lead_snp  = n_lead_snp,
  exposure_N           = N
)]
outc_sub <- outc[, .(
  outcome_source       = source,
  outcome_id           = trait_id,
  outcome_name         = trait_name,
  outcome_category     = category,
  outcome_n_lead_snp   = n_lead_snp,
  outcome_N            = N
)]
expo_sub[, k := 1L]
outc_sub[, k := 1L]
pair <- merge(expo_sub, outc_sub, by = "k", allow.cartesian = TRUE)[, k := NULL]
pair[, mr_direction := "protein2disease"]
# Sort by outcome then exposure so chunks share outcome reads (I/O cache friendly).
setorder(pair, outcome_id, exposure_id)
setcolorder(pair, c("mr_direction",
                    grep("^exposure_", names(pair), value = TRUE),
                    grep("^outcome_",  names(pair), value = TRUE)))
msg("total pairs: ", nrow(pair))

# ---- write mr_pair.csv.gz ----
out_pair <- file.path(MR_INFO_DIR, "mr_pair.csv.gz")
fwrite_gz(pair, out_pair)
msg("wrote ", nrow(pair), " rows -> ", out_pair)

# ---- todo = pairs not yet in the ledger ----
write_pair_todo(pair, DB_FILE, TMP_DIR, MR_INFO_DIR)

# ---- breakdown ----
msg("")
msg("--- breakdown ---")
msg("exposures by source:")
print(expo[, .N, by = source])
msg("exposures by category:")
print(expo[, .N, by = category])
msg("outcomes by source:")
print(outc[, .N, by = source])
msg("outcomes by category (top 15):")
print(outc[, .N, by = category][order(-N)][1:min(15, .N)])

msg("done in ", round(proc.time()[[3L]] - t0, 1), "s")
