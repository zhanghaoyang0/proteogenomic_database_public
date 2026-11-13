#!/usr/bin/env Rscript

# Per-trait metadata from GWAS clean sumstats only (no mr_pair.csv).
# Reads bf2/other/finngen clean/; for FinnGen prefers meta_fin_ukbb_mvp > meta_fin_ukbb > meta_fin
# per phenocode; writes source as BioFINDER2, Other, FinnGen-MVP-UKBB, FinnGen-UKBB.
# BF2 n_lead_snp from cis clump dirs (..._cis), same as 02_clump.sh CLUMP_CIS=1.

library(data.table)

source(path.expand("code/lib/common.R"))   # GWAS_BASE, CLUMP_SUBDIR, msg
BF1_DIR          <- file.path(GWAS_BASE, "bf1")
BF2_DIR          <- file.path(GWAS_BASE, "bf2")
OTHER_DIR        <- file.path(GWAS_BASE, "other")
FIN_DIR          <- file.path(GWAS_BASE, "finngen")
INFO_DIR         <- path.expand("./info")
BF_N_INFO_DIR    <- INFO_DIR
PQTL_FILE        <- file.path(INFO_DIR, "bf_protein/clump_lead_pqtl.csv.gz")

CLEAN_OTHER      <- file.path(OTHER_DIR, "clean")
CLEAN_FIN        <- file.path(FIN_DIR, "clean")

FIN_INFO_UKBB    <- file.path(FIN_DIR, "info", "finngen_R12_meta_analysis_mapping_with_definitions.tsv")
FIN_INFO_MVP     <- file.path(FIN_DIR, "info", "FinnGen_R12_MVP_UKBB_manifest.tsv")
FIN_INFO_R12     <- file.path(FIN_DIR, "info", "finngen_R12_manifest.tsv")
# Single source of truth for FinnGen category -> short disease class (editable
# data, not buried code). Includes semantic remaps (e.g. Diabetes/GI/ILD
# endpoint-groups folded into their parent ICD chapter).
FIN_CAT_MAP_FILE <- file.path(INFO_DIR, "finngen_category_map.csv")

dir.create(INFO_DIR, showWarnings = FALSE, recursive = TRUE)

t0 <- proc.time()[[3L]]

# Read the first per-variant N for sources where N is constant.
first_n_from_gwas <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  con <- tryCatch(gzfile(path, "rt"), error = function(e) NULL)
  if (is.null(con)) return(NA_integer_)
  on.exit(close(con), add = TRUE)
  hdr <- readLines(con, n = 1L, warn = FALSE)
  row1 <- readLines(con, n = 1L, warn = FALSE)
  if (!length(hdr) || !nzchar(hdr) || !length(row1)) return(NA_integer_)
  cols <- strsplit(hdr, "\t", fixed = TRUE)[[1L]]
  j <- match("N", cols)
  if (is.na(j)) return(NA_integer_)
  vals <- strsplit(row1, "\t", fixed = TRUE)[[1L]]
  if (length(vals) < j) return(NA_integer_)
  v <- suppressWarnings(as.numeric(vals[j]))
  if (!is.finite(v)) return(NA_integer_)
  as.integer(round(v))
}

# External "Other" meta-analyses can have variant-specific N. Use the maximum
# available sample size as the trait-level summary reported in trait_info.csv.
max_n_from_gwas <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  n <- tryCatch(
    fread(path, select = "N", showProgress = FALSE)$N,
    error = function(e) numeric()
  )
  n <- suppressWarnings(as.numeric(n))
  n <- n[is.finite(n)]
  if (!length(n)) return(NA_integer_)
  as.integer(round(max(n)))
}

# Count clumped lead SNP rows using fread instead of wc -l.
# Some PLINK .clumped files wrap long SP2 values across physical lines, so
# line counts can overestimate the number of lead SNPs.
clump_n_lead_vec <- function(paths, chunk = 200L, msg = NULL) {
  n <- length(paths)
  ans <- integer(n)
  fi <- file.info(paths)
  ok <- !is.na(fi$size) & fi$size > 0L
  ans[!ok] <- 0L
  n_miss <- sum(!ok)
  if (is.function(msg)) {
    msg("clump: ", n, " paths, ", sum(ok), " non-empty .clumped, ", n_miss, " missing/empty")
  }
  po <- paths[ok]
  if (!length(po)) return(ans)
  io <- which(ok)
  res_ok <- integer(length(po))
  n_chunk <- ceiling(length(po) / chunk)
  ch_i <- 0L
  for (s in seq(1L, length(po), by = chunk)) {
    ch_i <- ch_i + 1L
    e <- min(s + chunk - 1L, length(po))
    if (is.function(msg)) msg("clump fread: chunk ", ch_i, "/", n_chunk, " (", s, ":", e, " / ", length(po), ")")
    for (k in s:e) {
      # Stop, not 0: an unreadable file is not a trait without leads.
      dt <- tryCatch(fread(po[k], showProgress = FALSE),
                     error = function(e) stop("cannot read ", po[k], ": ", conditionMessage(e)))
      res_ok[k] <- nrow(dt)
    }
  }
  ans[io] <- res_ok
  ans
}

resolve_clump_file <- function(gwas_root, trait_id, clump_subdir, clump_type = NA_character_) {
  base <- if (!is.na(clump_type) && nzchar(clump_type))
    file.path(gwas_root, "clump", clump_type)
  else
    file.path(gwas_root, "clump")
  d1 <- file.path(base, clump_subdir)
  alt_subdir <- sub("Pthres", "pthres", clump_subdir, fixed = TRUE)
  d2 <- file.path(base, alt_subdir)
  # A missing dir (a renamed CLUMP_SUBDIR) would read as 0 leads for every trait.
  if (!dir.exists(d2)) assert_clump_dirs(d1)
  p1 <- file.path(d1, paste0(trait_id, ".clumped"))
  if (file.exists(p1)) return(p1)
  p2 <- file.path(d2, paste0(trait_id, ".clumped"))
  if (file.exists(p2)) return(p2)
  p1
}

# Load a FinnGen manifest (R12 / MVP / UKBB) into (fg_phenotype, outcome_name).
# Accepts either "phenocode" or "fg_phenotype" as the key column; picks first
# non-empty name from candidates.
read_fin_manifest <- function(path, msg) {
  if (!file.exists(path)) { msg("FinnGen lut: skip missing ", basename(path)); return(NULL) }
  d <- tryCatch(fread(path, showProgress = FALSE, quote = ""), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) { msg("FinnGen lut: skip empty ", basename(path)); return(NULL) }
  nm <- names(d); nm_l <- tolower(trimws(nm))
  ik <- match(TRUE, nm_l %in% c("phenocode", "fg_phenotype"))
  if (is.na(ik)) { msg("FinnGen lut: no key col in ", basename(path)); return(NULL) }
  name_cands <- c("phenotype", "name", "long_name", "definition")
  in_cols <- nm[nm_l %in% name_cands]
  if (!length(in_cols)) { msg("FinnGen lut: no name col in ", basename(path)); return(NULL) }
  name <- rep(NA_character_, nrow(d))
  for (cc in in_cols) {
    v <- trimws(as.character(d[[cc]]))
    hit <- is.na(name) | !nzchar(name)
    name[hit] <- ifelse(is.na(v[hit]) | !nzchar(v[hit]), NA_character_, v[hit])
  }
  data.table(fg_phenotype = trimws(as.character(d[[nm[ik]]])), outcome_name = name)
}

load_fin_name_lut <- function(msg) {
  rows <- Filter(Negate(is.null),
                 lapply(c(FIN_INFO_R12, FIN_INFO_MVP, FIN_INFO_UKBB), read_fin_manifest, msg = msg))
  if (!length(rows)) return(data.table(fg_phenotype = character(), outcome_name = character()))
  lut <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  lut[, .(outcome_name = {
    ok <- which(!is.na(outcome_name) & nzchar(outcome_name))
    if (length(ok)) outcome_name[ok[1L]] else outcome_name[1L]
  }), by = fg_phenotype]
}

# Phenocode for manifest match: same idea as finngen_clean.R (fg_phenotype has no finngen_R12_ prefix).
fg_pheno_from_trait_id <- function(trait_id) {
  x <- as.character(trait_id)
  x <- sub("_meta_fin_ukbb_mvp$", "", x)
  x <- sub("_meta_fin_ukbb$", "", x)
  x <- sub("_meta_fin$", "", x)
  sub("^finngen_R12_", "", x)
}

fg_source_rank <- function(trait_id) {
  x <- as.character(trait_id)
  fcase(
    grepl("_meta_fin_ukbb_mvp$", x), 1L,
    grepl("_meta_fin_ukbb$", x), 2L,
    grepl("_meta_fin$", x), 3L,
    default = 9L
  )
}

prefer_finngen_meta <- function(d, msg = NULL) {
  if (!nrow(d)) return(d)
  d <- copy(d)
  d[, fg_phenocode := vapply(trait_id, fg_pheno_from_trait_id, character(1))]
  d[, fg_rank := fg_source_rank(trait_id)]
  dup_n <- d[, .N, by = fg_phenocode][N > 1L, sum(N)]
  setorder(d, fg_phenocode, fg_rank, trait_id)
  keep <- d[, .SD[1L], by = fg_phenocode]
  if (is.function(msg)) {
    msg(
      "FinnGen prefer-meta: kept ", nrow(keep), "/", nrow(d),
      " rows after phenocode dedup",
      if (!is.na(dup_n) && dup_n > 0L) paste0(" (", dup_n, " duplicate-version rows seen)") else ""
    )
  }
  keep[, c("fg_phenocode", "fg_rank") := NULL]
  keep
}

# Shorten FinnGen phenotype strings: drop ", excluding ..." and rest, and a
# comma-less trailing "excluding ..." outside parentheses ("Autoimmune diseases
# excluding thyroid diseases"); "(controls excluding all cancers)" is left to
# strip_finngen_parenthetical_suffix.
strip_finngen_excluding <- function(ch) {
  trimws(sub("\\s+[Ee]xcluding[^()]*$", "", sub(",\\s*[Ee]xcluding.*$", "", ch)))
}

# Drop a trailing definition qualifier: "Gout, strict definition" -> "Gout".
# Endpoints that then share a name are reduced to one further down.
strip_finngen_definition <- function(ch) {
  trimws(sub(paste0(",\\s*(wide definition|strict definition|strict|",
                    "definitions combined|broad definition|including avohilmo)$"),
             "", ch, ignore.case = TRUE))
}

# Drop trailing parenthetical qualifiers, e.g. "(more control exclusions)".
# Vectorized: loops until no more trailing "(...)" remain.
strip_finngen_parenthetical_suffix <- function(ch) {
  x <- trimws(as.character(ch))
  repeat {
    y <- trimws(sub("\\s*\\([^)]*\\)\\s*$", "", x))
    if (identical(y, x)) return(y)
    x <- y
  }
}

# Rows: trait_id, gwas_path, source, gwas_root (for clump path)
from_clean_dir <- function(clean_dir, source_label, gwas_root) {
  if (!dir.exists(clean_dir)) {
    message("Skip (missing dir): ", clean_dir)
    return(data.table())
  }
  f <- list.files(clean_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
  if (!length(f)) {
    message("No .tsv.gz under ", clean_dir)
    return(data.table())
  }
  message("clean [", source_label, "]: ", length(f), " .tsv.gz -> ", clean_dir)
  tid <- sub("\\.tsv\\.gz$", "", basename(f))
  data.table(
    trait_id = tid,
    gwas_path = f,
    source = source_label,
    gwas_root = gwas_root
  )
}

msg("start trait_info: MR info dir -> ", INFO_DIR)
msg("clump subdir: ", CLUMP_SUBDIR, " (BF2 + _cis)")

fin_lut <- load_fin_name_lut(msg)
if (!nrow(fin_lut)) message("Note: FinnGen name manifest missing; fin trait_name := trait_id.")

# FinnGen raw category lookup, from the union of the three R12 manifests: the
# analysed outcomes include FinnGen-UKBB and FinnGen-UKBB-MVP meta-analyses whose
# phenocodes are not all in the base manifest (8 dental endpoints were, and fell to
# "Drug/other"). A phenocode with two different categories across manifests stops.
read_fin_category <- function(path) {
  if (!file.exists(path)) stop("Missing FinnGen manifest: ", path)
  d <- fread(path, showProgress = FALSE, quote = "")
  nl <- tolower(trimws(names(d)))
  k <- names(d)[match(TRUE, nl %in% c("phenocode", "fg_phenotype"))]
  cc <- names(d)[nl == "category"]
  if (is.na(k) || length(cc) != 1L) stop(basename(path), ": no phenocode/fg_phenotype or category column")
  d <- data.table(fg_phenotype = trimws(as.character(d[[k]])),
                  category = trimws(as.character(d[[cc]])))
  d[!is.na(fg_phenotype) & nzchar(fg_phenotype) & !is.na(category) & nzchar(category)]
}
fin_cat_lut <- unique(rbindlist(lapply(c(FIN_INFO_R12, FIN_INFO_MVP, FIN_INFO_UKBB), read_fin_category)))
cat_conflict <- fin_cat_lut[, .(n = uniqueN(category)), by = fg_phenotype][n > 1L]
if (nrow(cat_conflict)) stop(nrow(cat_conflict), " FinnGen phenocode(s) with conflicting categories across manifests, e.g. ",
                             paste(head(cat_conflict$fg_phenotype, 5L), collapse = ", "))
fin_cat_lut <- unique(fin_cat_lut, by = "fg_phenotype")
msg("FinnGen category lookup: ", nrow(fin_cat_lut), " phenocodes from 3 manifests")

add_type <- function(d, type) { if (nrow(d)) d[, clump_type := type]; d }
fin_part <- add_type(from_clean_dir(CLEAN_FIN, "fin", FIN_DIR), NA_character_)
fin_part <- prefer_finngen_meta(fin_part, msg)

parts <- list(
  add_type(from_clean_dir(file.path(BF1_DIR, "clean/protein"), "bf1", BF1_DIR), "protein"),
  add_type(from_clean_dir(file.path(BF1_DIR, "clean/pheno"), "bf1", BF1_DIR), "pheno"),
  add_type(from_clean_dir(file.path(BF2_DIR, "clean/protein"), "bf2", BF2_DIR), "protein"),
  add_type(from_clean_dir(file.path(BF2_DIR, "clean/pheno"), "bf2", BF2_DIR), "pheno"),
  add_type(from_clean_dir(CLEAN_OTHER, "other", OTHER_DIR), NA_character_),
  fin_part
)
tab <- rbindlist(parts, use.names = TRUE, fill = TRUE)
if (!nrow(tab)) stop("No clean .tsv.gz found")

msg("trait_name: protein = ENSG IDs; pheno = bf_pheno_info.csv; FinnGen = manifest; other = other_info.csv")
# Raw trait_id (strip "bf1__"/"bf2__" prefix) for joining lookup tables that key on unprefixed ids.
strip_bf_prefix <- function(x) sub("^bf[12]__", "", as.character(x))
tab[, raw_trait_id := strip_bf_prefix(trait_id)]
tab[, trait_name := as.character(trait_id)]

# protein: use ENSG ID(s) joined by "_" from bf_protein/ensg/assay_ensg.csv
ensg_lut <- fread(file.path(BF_N_INFO_DIR, "bf_protein/ensg/assay_ensg.csv"),
                  select = c("protein_assay", "ensembl_gene_id"), showProgress = FALSE)
ensg_lut[, ensg_label := gsub(",", "_", ensembl_gene_id, fixed = TRUE)]
is_protein <- tab$clump_type == "protein" & tab$source %in% c("bf1", "bf2")
protein_ids   <- as.character(tab$raw_trait_id[is_protein])
protein_ensg  <- ensg_lut$ensg_label[match(protein_ids, ensg_lut$protein_assay)]
protein_fluid <- ifelse(grepl("csf$", protein_ids, ignore.case = TRUE), "CSF", "Plasma")
protein_ensg[is.na(protein_ensg)] <- sub("__.*$", "", protein_ids[is.na(protein_ensg)])
protein_panel <- sub("^.+?__(.+)_o[12](csf|plasma)$", "\\1", protein_ids, ignore.case = TRUE)
protein_panel[protein_panel == protein_ids] <- ""
# Fluid omitted from trait_name — captured in `category`.
# CSF/plasma pairs will share trait_name; they're disambiguated by category.
protein_name <- trimws(paste(protein_ensg, protein_panel))
# Disambiguate only within the same fluid (still need unique label per category).
for (fl in unique(protein_fluid)) {
  idx <- which(protein_fluid == fl)
  nm <- protein_name[idx]
  dup <- duplicated(nm) | duplicated(nm, fromLast = TRUE)
  if (any(dup)) {
    pfx <- sub("__.*$", "", protein_ids[idx])
    protein_name[idx[dup]] <- trimws(paste(pfx[dup], protein_ensg[idx[dup]], protein_panel[idx[dup]]))
  }
}
tab[is_protein, trait_name := protein_name]

# other: use label from other_info.csv if available
other_info_path <- file.path(BF_N_INFO_DIR, "other_info.csv")
if (file.exists(other_info_path)) {
  other_lut <- fread(other_info_path, select = c("id", "label"), showProgress = FALSE)
  tab[source == "other", trait_name := {
    lbl <- other_lut$label[match(trait_id, other_lut$id)]
    ifelse(is.na(lbl) | !nzchar(lbl), as.character(trait_id), lbl)
  }]
}
pheno_info_path <- file.path(BF_N_INFO_DIR, "bf_pheno/bf_pheno_info.csv")
if (file.exists(pheno_info_path)) {
  pheno_lut <- fread(pheno_info_path, select = c("id", "label"), showProgress = FALSE)
  is_pheno <- tab$clump_type == "pheno"
  tab[is_pheno, trait_name := {
    lbl <- pheno_lut$label[match(raw_trait_id, pheno_lut$id)]
    ifelse(is.na(lbl) | !nzchar(lbl), as.character(trait_id), lbl)
  }]
}
tab[source == "fin", trait_name := {
  fg <- vapply(trait_id, fg_pheno_from_trait_id, character(1))
  key <- trimws(fg)
  nm <- fin_lut$outcome_name[match(key, fin_lut$fg_phenotype)]
  miss <- is.na(nm) | nm == ""
  if (any(miss)) {
    ii <- which(miss)
    nm[ii] <- fin_lut$outcome_name[match(toupper(key[ii]), toupper(fin_lut$fg_phenotype))]
    miss <- is.na(nm) | nm == ""
  }
  out <- ifelse(miss, as.character(trait_id), nm)
  ok <- !miss
  if (any(ok)) {
    out[ok] <- strip_finngen_definition(
      strip_finngen_parenthetical_suffix(strip_finngen_excluding(out[ok])))
  }
  out
}]
n_fin <- sum(tab$source == "fin")
n_named <- nrow(tab[source == "fin" & trait_name != as.character(trait_id)])
msg("FinnGen rows: ", n_fin, "; with manifest name: ", n_named, "; still trait_id as name: ", n_fin - n_named)

miss_fin <- tab[source == "fin" & trait_name == as.character(trait_id)]
if (nrow(miss_fin)) {
  miss_fin[, phenocode_key := vapply(trait_id, fg_pheno_from_trait_id, character(1))]
  miss_fin[, lut_name := fin_lut$outcome_name[match(phenocode_key, fin_lut$fg_phenotype)]]
  miss_fin[, in_lut := phenocode_key %in% fin_lut$fg_phenotype]
  miss_fin[, in_lut_ci := toupper(trimws(phenocode_key)) %in% toupper(fin_lut$fg_phenotype)]
  out_miss <- file.path(INFO_DIR, "trait_info_finngen_unmatched.tsv")
  fwrite(
    miss_fin[, .(trait_id, phenocode_key, lut_name, in_lut, in_lut_ci)],
    out_miss,
    sep = "\t"
  )
  msg(
    "FinnGen MISSING name: ", nrow(miss_fin), " traits -> ", out_miss,
    " (check phenocode_key vs finngen_R12_manifest phenocode)"
  )
  show <- min(15L, nrow(miss_fin))
  if (show > 0L) {
    msg("  first ", show, " examples:")
    for (i in seq_len(show)) {
      message("    ", miss_fin$trait_id[i], "  -> key=", miss_fin$phenocode_key[i])
    }
  }
}

# Formal source labels
tab[, source := fcase(
  source == "bf1", "BioFINDER1",
  source == "bf2", "BioFINDER2",
  source == "other", "Other",
  source == "fin",  "FinnGen meta-analysis",
  default = source
)]

cnt <- tab[, .N, by = source]
msg(
  "traits total: ", nrow(tab), " (by source: ",
  paste(paste0(cnt$source, "=", cnt$N), collapse = ", "), ")"
)

is_bf_protein <- tab$source %in% c("BioFINDER1", "BioFINDER2") & tab$clump_type == "protein"

# BF protein lead SNPs from clump_lead_pqtl.csv.gz, split cis / trans.
#
# n_lead_snp stays CIS-ONLY for backward compatibility -- 01_pair_protein2disease.R
# reads it as the instrument count. The two explicit
# columns exist because "n_lead_snp = 0" next to a fine-mapping plot is
# unreadable: fine-mapping runs genome-wide (finemap/01_finemapping.R), so an
# assay with a trans-only signal legitimately has 0 cis leads.
if (!file.exists(PQTL_FILE)) stop("Missing pQTL file: ", PQTL_FILE)
pqtl <- fread(PQTL_FILE, showProgress = FALSE)
pqtl_n <- pqtl[, .(n_cis = sum(is_cis %in% c(TRUE, "TRUE")),
                   n_trans = sum(!(is_cis %in% c(TRUE, "TRUE")))),
               by = .(cohort, protein_assay)]
pqtl_key <- paste(pqtl_n$cohort, pqtl_n$protein_assay)
tab[, `:=`(n_lead_snp = 0L, n_lead_cis = 0L, n_lead_trans = 0L)]
tab[is_bf_protein, c("n_lead_cis", "n_lead_trans") := {
  coh <- fifelse(source == "BioFINDER1", "bf1", "bf2")
  m <- match(paste(coh, raw_trait_id), pqtl_key)
  # No row in the clump table means no genome-wide significant lead at all.
  list(as.integer(fifelse(is.na(pqtl_n$n_cis[m]), 0L, pqtl_n$n_cis[m])),
       as.integer(fifelse(is.na(pqtl_n$n_trans[m]), 0L, pqtl_n$n_trans[m])))
}]
tab[is_bf_protein, n_lead_snp := n_lead_cis]
# The table is BF2 only: BF1 proteins are never clumped (replication looks BF2 signals
# up in them), so their lead counts are unknown, not 0.
if (any(pqtl_n$cohort != "bf2")) stop("clump_lead_pqtl now holds ", paste(setdiff(pqtl_n$cohort, "bf2"), collapse = ", "), "; revisit the BF1 NA below")
tab[is_bf_protein & source == "BioFINDER1", c("n_lead_snp", "n_lead_cis", "n_lead_trans") := NA_integer_]

# Non-BF-protein traits: n_lead_snp from genome-wide clump files
tab[!is_bf_protein, cf := mapply(
  resolve_clump_file, gwas_root, trait_id, CLUMP_SUBDIR, clump_type
)]
msg("counting lead SNPs for non-protein-BF traits ...")
if (any(!is_bf_protein))
  tab$n_lead_snp[!is_bf_protein] <- clump_n_lead_vec(tab$cf[!is_bf_protein], msg = msg)
msg("n_lead_snp summary: min=", min(tab$n_lead_snp, na.rm = TRUE), " max=", max(tab$n_lead_snp, na.rm = TRUE),
    " median=", median(tab$n_lead_snp, na.rm = TRUE), " NA (BF1 protein)=", sum(is.na(tab$n_lead_snp)))

# N for bf traits: fast lookup from info files (keyed on raw trait_id, unprefixed).
msg("reading N: bf from lookup, Other=max variant N, finngen=constant N ...")
# Protein N: from bf_protein/ensg/assay_ensg.csv (single source of truth; collapse multi-gene rows).
bf_protein_n <- fread(file.path(BF_N_INFO_DIR, "bf_protein/ensg/assay_ensg.csv"),
                      select = c("cohort", "protein_assay", "N"))
bf_protein_n <- unique(bf_protein_n, by = c("cohort", "protein_assay"))

# assay_ensg.csv used to carry a THIRD cohort, bf1_bf2_olink_ht -- the
# network-analysis bf1+bf2 bridged matrices, which have no GWAS at all
# (~/db/shared_data/gwas holds only bf1/bf2, and qc_panel_summary.csv lists only
# those two). The fifelse below is binary, so those rows were labelled BioFINDER2
# and collided with the real bf2 assays on (source, raw_trait_id); the merge then
# fanned every affected trait into two rows differing only in N (2,125 of them,
# all Protein plasma, e.g. bf2__SORT1__ExploreHT_o7plasma with N=1616 and
# N=2235). Keep only cohorts that actually have GWAS.
#
# 03_annot_ensg.R no longer emits that cohort, so this is now a transition
# guard: it keeps this script working against a copy of assay_ensg.csv written
# before that change. Drop bf1_bf2_olink_ht from the vector once the table has
# been regenerated everywhere; the stop() below is the part worth keeping.
GWAS_COHORTS     <- c("bf1", "bf2")
NON_GWAS_COHORTS <- c("bf1_bf2_olink_ht")
extra_cohorts <- setdiff(unique(bf_protein_n$cohort), c(GWAS_COHORTS, NON_GWAS_COHORTS))
if (length(extra_cohorts)) {
  stop("assay_ensg.csv has unexpected cohort(s): ", paste(extra_cohorts, collapse = ", "),
       " -- add to GWAS_COHORTS if it has GWAS, or to NON_GWAS_COHORTS if not. ",
       "The fifelse below is binary and would silently label it BioFINDER2.")
}
n_dropped <- nrow(bf_protein_n[!cohort %in% GWAS_COHORTS])
if (n_dropped) msg("dropping ", n_dropped, " assay_ensg rows from non-GWAS cohort(s): ",
                   paste(intersect(unique(bf_protein_n$cohort), NON_GWAS_COHORTS), collapse = ", "))
bf_protein_n <- bf_protein_n[cohort %in% GWAS_COHORTS]

bf_var_n     <- fread(file.path(BF_N_INFO_DIR, "bf_pheno/bf_var_n.csv"))  # source, var, n
bf_n_lut <- rbindlist(list(
  bf_protein_n[, .(source = fifelse(cohort == "bf1", "BioFINDER1", "BioFINDER2"), raw_trait_id = protein_assay, N = as.integer(N))],
  bf_var_n[,     .(source = fifelse(source == "bf1", "BioFINDER1", "BioFINDER2"), raw_trait_id = var, N = as.integer(n))]
))
tab <- merge(tab, bf_n_lut, by = c("source", "raw_trait_id"), all.x = TRUE)

# trait_id is the key every consumer relies on: the MR/coloc pair builders
# join on it, and the Shiny trait picker builds its choices straight from it.
# A duplicate here shows up as doubled dropdown entries, so assert the invariant
# instead of letting it leak downstream again.
if (anyDuplicated(tab$trait_id)) {
  dup_ids <- unique(tab$trait_id[duplicated(tab$trait_id)])
  stop("trait_id is not unique after the N merge: ", length(dup_ids), " duplicated id(s), e.g. ",
       paste(utils::head(dup_ids, 3), collapse = ", "),
       ". Check for another cohort or source collapsing onto an existing key.")
}

is_bf <- tab$source %in% c("BioFINDER1", "BioFINDER2")
n_bf_miss <- sum(is_bf & is.na(tab$N))
if (n_bf_miss) msg("bf lookup miss for ", n_bf_miss, " traits; will read N from gwas file")

need_n_idx <- which(!is_bf | is.na(tab$N))
msg("reading N from gwas file for ", length(need_n_idx), " traits (non-bf + bf with no lookup) ...")
gp <- tab$gwas_path
for (i in need_n_idx) {
  tab$N[i] <- if (tab$source[i] == "Other") {
    max_n_from_gwas(gp[i])
  } else {
    first_n_from_gwas(gp[i])
  }
}
na_n <- sum(is.na(tab$N))
if (na_n) msg("warning: N is NA for ", na_n, " traits")
tab[source %in% c("BioFINDER1", "BioFINDER2") & clump_type == "protein",
    category := ifelse(grepl("csf$", trait_id, ignore.case = TRUE), "Protein CSF", "Protein plasma")]
tab[source %in% c("BioFINDER1", "BioFINDER2") & clump_type == "pheno", category := "Phenotype"]
tab[source == "Other", category := "Published GWAS"]
# FinnGen raw category -> short class from the external map. The map is generated
# (with a fail-loud manifest-coverage check) by code/06_external_data/finngen_reclass.R --
# edit the taxonomy there and re-run it, not this file. A trait with no category in
# any manifest, or a raw category missing from the map, stops the run below.
if (!file.exists(FIN_CAT_MAP_FILE)) stop("Missing FinnGen category map: ", FIN_CAT_MAP_FILE)
fin_cat_map_dt <- fread(FIN_CAT_MAP_FILE, showProgress = FALSE)
FIN_CAT_MAP <- setNames(fin_cat_map_dt$class, fin_cat_map_dt$finngen_category)
msg("FinnGen category map: ", length(FIN_CAT_MAP), " raw -> ",
    uniqueN(fin_cat_map_dt$class), " short classes (", basename(FIN_CAT_MAP_FILE), ")")

# Source was relabeled to "FinnGen meta-analysis" above (line ~344); the old
# split labels never exist here, which silently left FinnGen category empty.
is_fin <- tab$source == "FinnGen meta-analysis"
n_fin0 <- sum(is_fin)
tab[, finngen_category_raw := NA_character_]
tab[is_fin, finngen_category_raw :=
      fin_cat_lut$category[match(vapply(trait_id, fg_pheno_from_trait_id, character(1)),
                                 fin_cat_lut$fg_phenotype)]]
unmatched <- tab[is_fin & is.na(finngen_category_raw), trait_id]
if (length(unmatched)) stop(length(unmatched), " FinnGen trait(s) with no category in any manifest: ",
                            paste(head(unmatched, 10L), collapse = ", "))
tab[is_fin, category := unname(FIN_CAT_MAP[finngen_category_raw])]
unmapped <- unique(tab[is_fin & (is.na(category) | !nzchar(category)), finngen_category_raw])
if (length(unmapped)) stop("FinnGen raw categories missing from ", basename(FIN_CAT_MAP_FILE), ": ",
                           paste(unmapped, collapse = " | "), " -- add them in code/06_external_data/finngen_reclass.R")
# No count-based folding: a class with few analysed traits keeps its name. A figure
# that needs fewer rows groups them at plot time, not here.
stopifnot(sum(tab$source == "FinnGen meta-analysis") == n_fin0)
msg("FinnGen classes (", n_fin0, " traits):")
print(tab[is_fin, .N, by = category][order(-N)])

# --- one definition per FinnGen disease ---------------------------------------
# All analysed FinnGen traits are FinnGen core endpoints, but the core set still
# holds several definitions of some diseases (Asthma vs Asthma as main diagnosis,
# T2D combined vs wide). Their names are equal once the definition qualifier is
# stripped above, so among endpoints sharing a name keep the one with the most
# cases, summed over the cohorts of the meta-analysis used. R12: 824 -> 812.
num0 <- function(x) { x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x }
mvp_m  <- fread(FIN_INFO_MVP)
ukbb_m <- fread(FIN_INFO_UKBB)
fin_cases <- rbind(
  mvp_m[, .(fg = fg_phenotype, meta = "mvp",
            cases = num0(fg_n_cases) + num0(ukbb_n_cases) + num0(MVP_AFR_n_cases) +
                    num0(MVP_EUR_n_cases) + num0(MVP_AMR_n_cases))],
  ukbb_m[, .(fg = fg_phenotype, meta = "ukbb",
             cases = num0(fg_n_cases) + num0(ukbb_n_cases))])
if (anyDuplicated(fin_cases, by = c("fg", "meta"))) stop("FinnGen manifest: duplicated phenocode")
fin <- tab[is_fin, .(trait_id, trait_name,
                     fg   = vapply(trait_id, fg_pheno_from_trait_id, character(1)),
                     meta = fifelse(grepl("_mvp$", trait_id), "mvp", "ukbb"))]
fin <- merge(fin, fin_cases, by = c("fg", "meta"), all.x = TRUE)
stopifnot(nrow(fin) == n_fin0)
# "Alzheimer's disease" and "Alzheimer disease" are one name.
fin[, key := tolower(gsub("['\u2019]s\\b", "", trait_name))]
fin[, n_def := .N, by = key]
if (anyNA(fin[n_def > 1L, cases]) || any(fin[n_def > 1L, cases] <= 0))
  stop("FinnGen: no case count for a trait with several definitions")
setorder(fin, key, -cases)
fin_drop <- fin[n_def > 1L & duplicated(key)]
msg("FinnGen definitions: ", uniqueN(fin[n_def > 1L, key]), " diseases with several; dropping ",
    nrow(fin_drop), " (fewer cases):")
print(fin[n_def > 1L, .(keep = !duplicated(key), trait_name, fg, cases)])
tab <- tab[!trait_id %chin% fin_drop$trait_id]
is_fin <- tab$source == "FinnGen meta-analysis"
stopifnot(sum(is_fin) == n_fin0 - nrow(fin_drop))
tab[, c("gwas_path", "gwas_root", "cf", "clump_type") := NULL]

# --- decompose the BF protein trait_id ---------------------------------------
# `bf2__GFAP__Oncology_o2csf` encodes cohort, assay and panel, and the panel token
# encodes platform and fluid. The Shiny app used to re-derive all five by string
# surgery in four separate places, which is how a display string once came back
# as a fluid and how NuLISA kept a redundant "_CSF" suffix. Store the fields
# instead, so consumers read columns rather than decode an id.
#
# Panel tokens follow three upstream conventions, all ending in the fluid:
#   Olink Explore   Oncology_o2csf          ExploreHT_o7plasma
#   Olink Target    NEUROLOGY_v8001_o1csf   INFLAMMATION_v3004_o1plasma
#   NULISA          CNS_nulisa2_CSF
#   SomaScan        (no panel token at all)
# Version tokens (v8001, nulisa2) are kept -- they name the panel, not the fluid.
# NOT is_bf_protein: that vector was built at the top of the file, and the
# bf_n_lut merge above re-sorts `tab`, so it no longer lines up row for row
# (289 rows off). Rebuild the mask from `category`, which was just assigned and
# is what every consumer filters on.
is_prot <- grepl("^Protein", tab$category)
tab[, c("cohort", "assay", "panel", "platform", "fluid") :=
      list(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)]

bf_parts <- strsplit(tab$trait_id[is_prot], "__", fixed = TRUE)
bf_nth <- function(k) vapply(bf_parts, function(p)
  if (length(p) >= k) p[k] else NA_character_, character(1))

raw_panel <- bf_nth(3L)
raw_panel[is.na(raw_panel)] <- ""
tab[is_prot, `:=`(
  cohort   = fifelse(source == "BioFINDER1", "BF1", "BF2"),
  assay    = bf_nth(2L),
  panel    = sub("_(o[0-9]+)?(csf|plasma)$", "", raw_panel, ignore.case = TRUE),
  platform = fifelse(!nzchar(raw_panel), "SomaScan",
             fifelse(grepl("nulisa", raw_panel, ignore.case = TRUE), "NULISA", "Olink")),
  fluid    = fifelse(grepl("CSF", category, ignore.case = TRUE), "CSF", "Plasma")
)]

# Fail loudly: a new panel naming convention must not silently produce a label
# that still carries its fluid suffix, which is what the app shipped for months.
bad_assay <- tab[is_prot & (is.na(assay) | !nzchar(assay)), trait_id]
if (length(bad_assay))
  stop("cannot read assay from ", length(bad_assay), " trait_id(s), e.g. ",
       paste(head(bad_assay, 3L), collapse = ", "))
bad_panel <- tab[is_prot & grepl("(csf|plasma)$", panel, ignore.case = TRUE), trait_id]
if (length(bad_panel))
  stop("panel still ends in a fluid suffix for ", length(bad_panel),
       " trait_id(s) -- add the new convention to the sub() above, e.g. ",
       paste(head(bad_panel, 3L), collapse = ", "))

msg("BF protein decomposition: ", sum(is_prot), " assays | platforms ",
    paste(sort(unique(na.omit(tab$platform))), collapse = "/"), " | panels ",
    uniqueN(tab$panel[is_prot]))

trait_info <- tab[, .(trait_id, trait_name, source, category, cohort, assay, panel,
                      platform, fluid, N, n_lead_snp, n_lead_cis, n_lead_trans,
                      finngen_category_raw)]
if (anyDuplicated(trait_info$trait_id)) stop("trait_id not unique in trait_info")
setorder(trait_info, source, trait_id)

# Uniqueness check on (trait_name, category) — label doesn't need to be
# unique on its own; trait_id is the key. (trait_name, category) catches
# real collisions while ignoring expected CSF/plasma or BF1/BF2 sharing.
dup_keys <- trait_info[, .N, by = .(trait_name, category)][N > 1L]
if (nrow(dup_keys)) {
  msg("warning: (trait_name, category) duplicated for ", nrow(dup_keys), " combo(s):")
  print(dup_keys)
} else {
  msg("trait_name: unique per category")
}

out_tsv <- file.path(INFO_DIR, "trait_info.csv")
fwrite(trait_info, out_tsv)
msg("wrote ", nrow(trait_info), " rows -> ", out_tsv)

msg("done; total elapsed ", round(proc.time()[[3L]] - t0, 1), "s")
