#!/usr/bin/env Rscript
# RUN LOCATION: LOCAL (Mac) — reads ~/db/secure_data/ (individual-level data).
# Outputs: info/bf_protein/ensg/{unmatched,assay_ensg,assay_cis_windows}.csv  (rename.csv is hand-curated input)
#
# Build protein -> ENSG/cis-window lookup for BF protein GWAS phenotypes.
# Source protein names from local PLINK .phe headers, then map gene symbols using
# offline GRCh37 Ensembl GTF, HGNC aliases, plus optional manual overrides.
library(data.table)

# ============================================================
# Step 0: set up
# ============================================================
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
WB_DIR <- normalizePath(file.path(script_dir, "../.."))

# Read from post-QC matrices (one column per protein) rather than .phe files.
# .phe is filtered to N>1000 non-missing, which drops every protein when a panel
# has <1000 baseline samples (bf1 Olink panels), so .phe-derived coverage is
# incomplete. QC matrices contain the full post-QC protein set.
PROTEIN_QC_DIR <- path.expand("~/db/secure_data/bf/protein/qc")
WB_INFO_DIR     <- file.path(WB_DIR, "info")
PROTEIN_INFO_DIR <- file.path(WB_INFO_DIR, "bf_protein")
ENSG_DIR <- file.path(PROTEIN_INFO_DIR, "ensg")
source(file.path(WB_DIR, "code/lib/cis_windows.R"))   # CIS_WINDOW_BP

GTF_GRCH37_GZ <- path.expand("~/db/shared_data/ref_geno/ensembl/Homo_sapiens.GRCh37.87.chr.gtf.gz")
HGNC_TSV <- path.expand("~/db/shared_data/ref_geno/ensembl/hgnc_complete_set.txt")
SOMA_META <- path.expand(
  "~/db/secure_data/bf/pheno/raw/haoyang__20260429_234900_somalogic/plasma_somalogic4.1_metadata_Jan2024.csv"
)
UNMATCH_CSV <- file.path(ENSG_DIR, "unmatched.csv")
RENAME_CSV <- file.path(ENSG_DIR, "rename.csv")

dir.create(WB_INFO_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PROTEIN_INFO_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ENSG_DIR, recursive = TRUE, showWarnings = FALSE)

out_all       <- file.path(ENSG_DIR, "assay_ensg.csv")
out_windows   <- file.path(ENSG_DIR, "assay_cis_windows.csv")
out_unmatched <- UNMATCH_CSV

if (!file.exists(GTF_GRCH37_GZ)) stop("GRCh37 GTF not found: ", GTF_GRCH37_GZ)
if (!file.exists(HGNC_TSV)) stop("HGNC complete set not found: ", HGNC_TSV)
if (!dir.exists(PROTEIN_QC_DIR)) stop("Protein QC dir not found: ", PROTEIN_QC_DIR)

if (!file.exists(RENAME_CSV)) {
  fwrite(
    data.table(source = character(), original_name = character(), rename = character()),
    RENAME_CSV
  )
  message("Created empty rename CSV: ", RENAME_CSV)
}

# ============================================================
# Step 1: process reference lookups
# ============================================================
extract_attr <- function(attr, key) {
  m <- regexpr(paste0(key, " \"([^\"]+)\""), attr, perl = TRUE)
  out <- rep(NA_character_, length(attr))
  ok <- m > 0
  if (any(ok)) {
    out[ok] <- sub(paste0(".*", key, " \"([^\"]+)\".*"), "\\1", attr[ok], perl = TRUE)
  }
  out
}

read_gene_gtf <- function(path) {
  gtf_cols <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  g <- fread(
    cmd = paste("gzip -cd", shQuote(path)),
    sep = "\t", header = FALSE, quote = "",
    col.names = gtf_cols,
    showProgress = FALSE
  )
  g <- g[!startsWith(seqname, "#")]
  g <- g[feature == "gene"]
  g[, gene := extract_attr(attribute, "gene_name")]
  g[, ensembl_gene_id := extract_attr(attribute, "gene_id")]
  g <- g[nzchar(gene) & nzchar(ensembl_gene_id)]
  g[, ensembl_gene_id := sub("\\..*$", "", ensembl_gene_id)]
  unique(g[, .(
    seqname = sub("^chr", "", as.character(seqname), ignore.case = TRUE),
    start = as.integer(start),
    end = as.integer(end),
    strand = as.character(strand),
    gene,
    gene_upper = toupper(gene),
    ensembl_gene_id
  )])
}

gtf <- read_gene_gtf(GTF_GRCH37_GZ)
message("Loaded GRCh37 GTF: ", GTF_GRCH37_GZ)

match_key <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
}

split_hgnc_field <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || !nzchar(x)) return(character())
  vals <- unlist(strsplit(x, "\\|"))
  vals <- trimws(vals)
  vals[nzchar(vals)]
}

hgnc <- fread(HGNC_TSV, sep = "\t", quote = "", fill = TRUE, showProgress = FALSE)
required_hgnc_cols <- c("symbol", "prev_symbol", "alias_symbol", "entrez_id", "ensembl_gene_id")
if (!all(required_hgnc_cols %in% names(hgnc))) {
  stop("HGNC complete set needs columns: ", paste(required_hgnc_cols, collapse = ", "))
}
hgnc[, ensembl_gene_id := sub("\\..*$", "", trimws(as.character(ensembl_gene_id)))]
hgnc <- hgnc[nzchar(ensembl_gene_id) & nzchar(symbol)]

hgnc_alias <- rbindlist(lapply(seq_len(nrow(hgnc)), function(i) {
  aliases <- unique(c(
    trimws(as.character(hgnc$symbol[i])),
    split_hgnc_field(hgnc$prev_symbol[i]),
    split_hgnc_field(hgnc$alias_symbol[i])
  ))
  aliases <- aliases[nzchar(aliases)]
  data.table(
    alias = aliases,
    alias_key = match_key(aliases),
    hgnc_symbol = hgnc$symbol[i],
    ensembl_gene_id = hgnc$ensembl_gene_id[i]
  )
}), use.names = TRUE, fill = TRUE)
hgnc_alias <- unique(hgnc_alias[nzchar(alias_key)])
message("Loaded HGNC aliases: ", HGNC_TSV)

hgnc_entrez <- unique(hgnc[nzchar(entrez_id), .(
  entrez_id = trimws(as.character(entrez_id)),
  hgnc_symbol = symbol,
  ensembl_gene_id
)])

hgnc_uniprot <- unique(hgnc[nzchar(uniprot_ids), .(
  uniprot_id = trimws(as.character(uniprot_ids)),
  hgnc_symbol = symbol,
  ensembl_gene_id
)])

missing_symbol_values <- c("", "NA", "NAN", "NONE", "NULL", "MISSING")
valid_symbol <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) & !(toupper(x) %in% missing_symbol_values)
}

read_commentable_csv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  if (!length(lines)) return(data.table())
  fread(
    text = paste(lines, collapse = "\n"),
    fill = TRUE,
    blank.lines.skip = TRUE,
    showProgress = FALSE
  )
}

canonical_gene <- function(x) {
  x <- trimws(as.character(x))
  out <- x
  for (i in seq_along(x)) {
    if (is.na(x[i]) || !nzchar(x[i])) next
    hit <- gtf[gene_upper == toupper(x[i]), unique(gene)]
    if (length(hit) == 1L) out[i] <- hit[1]
  }
  out
}

# ============================================================
# Step 2: collect proteins from .phe headers
# ============================================================
platform_from_tag <- function(tag) {
  if (grepl("nulisa", tag, ignore.case = TRUE)) return("nulisa")
  if (grepl("somalogic", tag, ignore.case = TRUE)) return("somalogic")
  if (grepl("olink|_o[0-9]+(csf|plasma)$", tag, ignore.case = TRUE)) return("olink")
  "protein"
}

source_from_tag <- platform_from_tag

qc_files <- list.files(
  PROTEIN_QC_DIR,
  pattern = "^bf[12]_.+_protein\\.csv$",
  full.names = TRUE
)
qc_files <- qc_files[!grepl("/\\._", qc_files)]
if (!length(qc_files)) stop("No protein QC files found in: ", PROTEIN_QC_DIR)

protein_list <- lapply(qc_files, function(path) {
  base <- basename(path)
  cohort <- sub("^(bf[12])_.*$", "\\1", base)
  tag <- sub("^bf[12]_", "", sub("_protein\\.csv$", "", base))
  dt <- fread(path, sep = ",", showProgress = FALSE, na.strings = c("NA", "-9", ""))
  # QC matrix has `sid` (sample id) + one column per protein assay.
  protein_cols <- setdiff(names(dt), c("sid", "mean_protein"))
  if (!length(protein_cols)) return(NULL)
  n_vals <- vapply(protein_cols, function(col) sum(!is.na(dt[[col]])), integer(1))
  data.table(
    cohort = cohort,
    tag = tag,
    source = source_from_tag(tag),
    platform = platform_from_tag(tag),
    protein_assay = protein_cols,
    N = n_vals
  )
})

# This used to also scan projects/proteogenomic_network/knn_imputed/ and emit a third
# cohort, bf1_bf2_olink_ht, so that the bridged bf1+bf2 Olink HT matrices got
# annotated. Those 3,162 rows were redundant: annotation is a property of the
# assay, not of a cohort or a QC pass, and bf1 + bf2 already cover all 2,225
# assays in the bridged matrices (the only mismatches were four hyphen/underscore
# spellings -- H1-1, HLA-A, INS-CPEPTIDE, NT-proBNP). The project now looks its
# own columns up in this table; see projects/proteogenomic_network/code/03_format_input.py.
summary <- unique(rbindlist(protein_list, use.names = TRUE, fill = TRUE))

# QC matrices use `{assay}__{tag}` columns (Olink/NULISA) or bare `seq_{seqid}`
# (SomaLogic) -- the platform metric suffix (NPX/NPQ/ANML) is dropped in 01.
# Strip everything from the first `__` onward to recover the bare assay/seq id.
summary[, original_name := sub("__.*$", "", protein_assay)]

# SomaLogic assay IDs are not gene names; use official metadata first.
summary[, soma_seqid := NA_character_]
summary[platform == "somalogic", soma_seqid := sub("^seq_", "", original_name)]
summary[platform == "somalogic", soma_seqid := sub("_", "-", soma_seqid)]

if (file.exists(SOMA_META)) {
  soma_meta <- fread(SOMA_META, showProgress = FALSE)
  if (all(c("SeqId", "EntrezGeneSymbol", "EntrezGeneID", "Target") %in% names(soma_meta))) {
    soma_meta <- unique(soma_meta[, .(
      soma_seqid = as.character(SeqId),
      soma_gene = as.character(EntrezGeneSymbol),
      soma_entrez_id = trimws(as.character(EntrezGeneID)),
      soma_target = trimws(as.character(Target)),
      soma_uniprot = trimws(as.character(UniProt))
    )])
    summary <- merge(summary, soma_meta, by = "soma_seqid", all.x = TRUE, sort = FALSE)
  } else {
    summary[, soma_gene := NA_character_]
    summary[, soma_entrez_id := NA_character_]
    summary[, soma_target := NA_character_]
    summary[, soma_uniprot := NA_character_]
    warning("SomaLogic metadata lacks SeqId/EntrezGeneSymbol/EntrezGeneID/Target columns: ", SOMA_META)
  }
} else {
  summary[, soma_gene := NA_character_]
  summary[, soma_entrez_id := NA_character_]
  summary[, soma_target := NA_character_]
  summary[, soma_uniprot := NA_character_]
  warning("SomaLogic metadata not found: ", SOMA_META)
}

summary[, gene_raw := fifelse(
  platform == "somalogic" & valid_symbol(soma_gene),
  soma_gene,
  fifelse(platform == "somalogic" & valid_symbol(soma_target), soma_target, original_name)
)]

summary[, gene_from_entrez := NA_character_]
summary[
  platform == "somalogic" & !valid_symbol(soma_gene) & !is.na(soma_entrez_id) & nzchar(soma_entrez_id),
  gene_from_entrez := hgnc_entrez$ensembl_gene_id[match(soma_entrez_id, hgnc_entrez$entrez_id)]
]

summary[, gene_from_uniprot := NA_character_]
summary[
  platform == "somalogic" & !is.na(soma_uniprot) & nzchar(soma_uniprot),
  gene_from_uniprot := hgnc_uniprot$ensembl_gene_id[match(soma_uniprot, hgnc_uniprot$uniprot_id)]
]

# ============================================================
# Step 3: automatic and manual gene-symbol normalization
# ============================================================
map_us <- c(
  "HLA_A" = "HLA-A",
  "HLA_DRA" = "HLA-DRA"
)

map_tok <- c(
  "HLAA"    = "HLA-A",
  "HLADRA"  = "HLA-DRA",
  "CXCL8"   = "IL8"
)

map_multi <- c(
  "AMY1AAMY1BAMY1C" = "AMY1A_AMY1B_AMY1C",
  "CKMT1ACKMT1B"    = "CKMT1A_CKMT1B",
  "DEFA1DEFA1B"      = "DEFA1_DEFA1B",
  "DEFB4ADEFB4B"     = "DEFB4A_DEFB4B",
  "EBI3IL27"         = "EBI3_IL27",
  "FUT3FUT5"        = "FUT3_FUT5",
  "IFNA1IFNA13"      = "IFNA1_IFNA13",
  "IL12AIL12B"      = "IL12A_IL12B",
  "LGALS7LGALS7B"   = "LGALS7_LGALS7B",
  "MICBMICA"        = "MICB_MICA",
  "SPACA5SPACA5B"   = "SPACA5_SPACA5B"
)

normalize_gene <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  normalize_one <- function(g) {
    if (is.na(g) || !nzchar(g)) return(NA_character_)
    if (g %in% names(map_us)) {
      key <- map_us[g]
    } else {
      key <- toupper(gsub("_", "", g, fixed = TRUE))
    }
    if (key %in% names(map_tok)) key <- map_tok[key]
    if (key %in% names(map_multi)) key <- map_multi[key]
    key
  }
  for (i in seq_along(x)) {
    g <- x[i]
    if (is.na(g) || !nzchar(g)) next
    if (grepl("\\|", g)) {
      key <- paste(vapply(strsplit(g, "\\|")[[1]], normalize_one, character(1)), collapse = "_")
    } else {
      key <- normalize_one(g)
    }
    parts <- strsplit(key, "_", fixed = TRUE)[[1]]
    parts <- canonical_gene(parts)
    out[i] <- paste(parts[nzchar(parts)], collapse = "_")
  }
  out
}

summary[, gene := normalize_gene(gene_raw)]
summary[!is.na(gene_from_entrez) & nzchar(gene_from_entrez), gene := gene_from_entrez]
summary[!is.na(gene_from_uniprot) & nzchar(gene_from_uniprot) & (is.na(gene) | !nzchar(gene)), gene := gene_from_uniprot]
summary[, override := FALSE]

rename_match_key <- match_key

overrides <- read_commentable_csv(RENAME_CSV)
required_override_cols <- c("source", "original_name", "rename")
if (!all(required_override_cols %in% names(overrides))) {
  stop("Rename CSV needs columns: ", paste(required_override_cols, collapse = ", "))
}
overrides <- overrides[
  !is.na(rename) & nzchar(trimws(rename)) &
    !is.na(original_name) & nzchar(trimws(original_name))
]
overrides[, source := trimws(as.character(source))]
overrides[, original_name := trimws(as.character(original_name))]
overrides[, rename := trimws(as.character(rename))]

if (nrow(overrides)) {
  summary[, original_name_key := rename_match_key(original_name)]
  summary[, protein_assay_key := rename_match_key(protein_assay)]
  summary[, gene_raw_key := rename_match_key(gene_raw)]
  for (i in seq_len(nrow(overrides))) {
    src <- overrides$source[i]
    orig <- overrides$original_name[i]
    orig_key <- rename_match_key(orig)
    idx <- summary$original_name == orig |
      summary$gene_raw == orig |
      summary$protein_assay == orig |
      (!is.na(summary$soma_uniprot) & summary$soma_uniprot == orig) |
      summary$original_name_key == orig_key |
      summary$gene_raw_key == orig_key |
      summary$protein_assay_key == orig_key
    if (!is.na(src) && nzchar(src) && !(src %in% c("*", "all", "ALL"))) {
      idx <- idx & summary$source == src
    }
    if (any(idx)) {
      summary[idx, gene := trimws(overrides$rename[i])]
      summary[idx, override := TRUE]
    } else {
      message("Override did not match any assay: source=", src, ", original_name=", orig)
    }
  }
  summary[, c("original_name_key", "protein_assay_key", "gene_raw_key") := NULL]
}

# Manual coordinates for genes absent from this GRCh37 GTF.
manual_missing <- data.table(
  gene = c("PECAM1", "PRSS2"),
  chr = c("17", "7"),
  start = c(62396775L, 142468247L),
  end = c(62404856L, 142475461L),
  # NA rather than a guess: these two are absent from this GTF, so nothing here
  # states their orientation. Consumers that need a strand see it is unknown.
  strand = c(NA_character_, NA_character_),
  ensembl_gene_id = c("ENSG00000261371", "ENSG00000275896")
)

collapse_ensg <- function(x) {
  x <- paste(stats::na.omit(as.character(x)), collapse = ",")
  x <- unlist(strsplit(x, ",", fixed = TRUE))
  x <- unique(trimws(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(sort(x), collapse = ",")
}

genes_needed <- unique(unlist(strsplit(summary$gene[nzchar(summary$gene)], "_", fixed = TRUE)))
gene_coords <- gtf[gene %in% genes_needed, .(
  gene,
  chr = as.character(seqname),
  start = as.integer(start),
  end = as.integer(end),
  strand = as.character(strand),
  ensembl_gene_id
)]
gene_coords <- rbindlist(list(
  gene_coords[!gene %in% manual_missing$gene],
  manual_missing[gene %in% genes_needed]
), use.names = TRUE, fill = TRUE)

primary_chr <- c(as.character(1:22), "X", "Y", "MT")
gene_coords[, chr := sub("^chr", "", as.character(chr), ignore.case = TRUE)]
gene_coords <- gene_coords[, {
  sub <- .SD
  on_primary <- sub$chr %in% primary_chr
  if (any(on_primary)) sub <- sub[on_primary]
  # A name with several ENSG (e.g. CKS1B, LSP1) is placed at its first ENSG in sorted
  # order, merging only the others within that one's cis window; min/max over all of
  # them gave windows of up to 94 Mb. All ENSG are kept, as for multi-part assays below.
  ensg_all <- collapse_ensg(sub$ensembl_gene_id)
  anchor <- if (is.na(ensg_all)) sub[1] else sub[grepl(sub(",.*$", "", ensg_all), ensembl_gene_id, fixed = TRUE)][1]
  near <- sub[chr == anchor$chr & end >= anchor$start - CIS_WINDOW_BP & start <= anchor$end + CIS_WINDOW_BP]
  .(
    chr = anchor$chr,
    start = min(near$start),
    end = max(near$end),
    strand = anchor$strand,
    ensembl_gene_id = ensg_all
  )
}, by = gene]

gtf_by_ensg <- gtf[seqname %in% primary_chr, .(
  gene = gene[1],
  chr = as.character(seqname[1]),
  start = min(as.integer(start)),
  end = max(as.integer(end)),
  strand = as.character(strand[1]),
  ensembl_gene_id = collapse_ensg(ensembl_gene_id)
), by = ensembl_gene_id]

resolve_gene_part <- function(part) {
  part <- trimws(as.character(part))
  if (is.na(part) || !nzchar(part)) return(NULL)

  if (grepl("^ENSG", part)) {
    ensg <- sub("\\..*$", "", part)
    hit <- gtf_by_ensg[ensembl_gene_id == ensg]
    if (nrow(hit)) {
      return(list(coord = hit[, .(chr, start, end, strand, ensembl_gene_id)], source = "ensg_grch37"))
    }
    return(NULL)
  }

  hit <- gene_coords[gene == part]
  if (nrow(hit)) {
    return(list(coord = hit[, .(chr, start, end, strand, ensembl_gene_id)], source = "name_grch37"))
  }

  hit_hgnc <- hgnc_alias[alias_key == match_key(part), unique(ensembl_gene_id)]
  hit_hgnc <- hit_hgnc[nzchar(hit_hgnc)]
  if (length(hit_hgnc)) {
    hit37 <- gtf_by_ensg[ensembl_gene_id %in% hit_hgnc]
    if (nrow(hit37) == 1L) {
      return(list(
        coord = hit37[, .(chr, start, end, strand, ensembl_gene_id)],
        source = "hgnc_alias_ensg_grch37"
      ))
    }
  }

  NULL
}

# ============================================================
# Step 4: attach coordinates and write outputs
# ============================================================
coord_by_row <- lapply(summary$gene, function(g) {
  genes <- strsplit(g, "_", fixed = TRUE)[[1]]
  resolved <- lapply(genes, resolve_gene_part)
  if (any(vapply(resolved, is.null, logical(1)))) {
    sub <- data.table()
    sources <- character()
  } else {
    sub <- rbindlist(lapply(seq_along(resolved), function(i) resolved[[i]]$coord[, part := genes[i]]),
                     use.names = TRUE, fill = TRUE)
    sources <- unique(vapply(resolved, function(r) r$source, character(1)))
  }
  if (nrow(sub) == 0 || any(is.na(sub$chr))) {
    return(list(
      gene_chr = NA_character_, gene_start = NA_integer_, gene_end = NA_integer_,
      gene_strand = NA_character_,
      cis_start = NA_integer_, cis_end = NA_integer_,
      ensembl_gene_id = NA_character_,
      lookup_source = NA_character_, windows = NULL
    ))
  }
  # cis windows: one per part (gene +/- CIS_WINDOW_BP), overlapping ones merged, so a
  # complex such as EBI3_IL27 has one window per subunit and a pQTL at either gene is cis.
  # They go to assay_cis_windows.csv; is_cis is decided there (code/lib/cis_windows.R).
  w <- sub[, .(chr, s = pmax(0, start - CIS_WINDOW_BP), e = end + CIS_WINDOW_BP, part)][order(chr, s)]
  w[, grp := cumsum(c(TRUE, s[-1] > cummax(e)[-.N])), by = chr]
  windows <- w[, .(cis_start = as.integer(min(s)), cis_end = as.integer(max(e)),
                   genes = paste(unique(part), collapse = ";")), by = .(chr, grp)][, grp := NULL]
  # The single gene_* / cis_* columns of assay_ensg.csv describe the assay's named gene only
  # (display, region extraction). A multi-part assay is anchored on the part that gives it
  # its gene name below (the first ENSG in sorted order), merging only parts within that
  # part's window; min/max over all parts gave windows of up to 190 Mb.
  ensg_all <- collapse_ensg(sub$ensembl_gene_id)
  anchor <- if (is.na(ensg_all)) sub[1] else sub[grepl(sub(",.*$", "", ensg_all), ensembl_gene_id, fixed = TRUE)][1]
  if (is.na(anchor$chr)) stop("no anchor part for ", g)
  near <- sub[chr == anchor$chr & end >= anchor$start - CIS_WINDOW_BP & start <= anchor$end + CIS_WINDOW_BP]
  gs <- min(near$start)
  ge <- max(near$end)
  list(
    gene_chr = anchor$chr,
    gene_start = gs,
    gene_end = ge,
    # "+"/"-" from the GTF, NA when unknown; the anchor part's orientation.
    gene_strand = as.character(anchor$strand),
    cis_start = as.integer(max(0, gs - CIS_WINDOW_BP)),
    cis_end = as.integer(ge + CIS_WINDOW_BP),
    # every part's ENSG, so UniProt matching (04) still sees the whole complex
    ensembl_gene_id = ensg_all,
    lookup_source = paste(sort(sources), collapse = ";"),
    windows = windows
  )
})

summary[, gene_chr := sapply(coord_by_row, function(r) r$gene_chr)]
summary[, gene_start := sapply(coord_by_row, function(r) r$gene_start)]
summary[, gene_end := sapply(coord_by_row, function(r) r$gene_end)]
summary[, gene_strand := sapply(coord_by_row, function(r) r$gene_strand)]
summary[, cis_start := sapply(coord_by_row, function(r) r$cis_start)]
summary[, cis_end := sapply(coord_by_row, function(r) r$cis_end)]
summary[, ensembl_gene_id := sapply(coord_by_row, function(r) r$ensembl_gene_id)]
summary[, lookup_source := sapply(coord_by_row, function(r) r$lookup_source)]
cis_windows <- rbindlist(lapply(seq_along(coord_by_row), function(i) {
  w <- coord_by_row[[i]]$windows
  if (is.null(w)) return(NULL)
  cbind(summary[i, .(cohort, tag, protein_assay)], w)
}))
stopifnot(nrow(cis_windows) >= sum(!is.na(summary$gene_chr)),
          !anyDuplicated(cis_windows[, .(cohort, tag, protein_assay, chr, cis_start)]))


# Replace gene with the canonical Ensembl gene name looked up from the primary ENSG.
summary[, gene := {
  first_ensg <- sub(",.*$", "", ensembl_gene_id)
  canonical <- gtf_by_ensg$gene[match(first_ensg, gtf_by_ensg$ensembl_gene_id)]
  fifelse(!is.na(canonical), canonical, gene)
}]

all_cols <- c(
  "cohort", "tag", "protein_assay", "gene", "ensembl_gene_id",
  "gene_chr", "gene_start", "gene_end", "gene_strand", "cis_start", "cis_end", "N"
)

# Compare against the previous run before overwriting. The tag is derived from
# the QC matrix filename, so renaming or losing a matrix silently empties that
# tag instead of failing, with no error. Shout about it here.
#
# The first run after the bf1_bf2_olink_ht rows were dropped will legitimately
# report those two (cohort, tag) pairs as gone. That is the intended change, not
# a missing matrix.
if (file.exists(out_all)) {
  prev_n <- fread(out_all, select = c("cohort", "tag"))[, .(n_prev = .N), by = .(cohort, tag)]
  new_n <- summary[, .(n_new = .N), by = .(cohort, tag)]
  cmp <- merge(prev_n, new_n, by = c("cohort", "tag"), all.x = TRUE)

  gone <- cmp[is.na(n_new)][order(-n_prev)]
  if (nrow(gone) > 0) {
    message("WARNING: ", nrow(gone), " (cohort, tag) pair(s) present last run are now absent.")
    message("  A renamed or missing QC matrix looks exactly like this. Check before trusting the output:")
    print(gone[, .(cohort, tag, n_prev)])
  }

  shrunk <- cmp[!is.na(n_new) & n_new < 0.8 * n_prev][order(n_new / n_prev)]
  if (nrow(shrunk) > 0) {
    message("WARNING: ", nrow(shrunk), " tag(s) lost more than 20% of their assays:")
    print(shrunk)
  }
}

fwrite(summary[, ..all_cols], out_all)
message("Wrote ", nrow(summary), " rows to ", out_all)
fwrite(cis_windows, out_windows)
message("Wrote ", nrow(cis_windows), " cis windows (", cis_windows[, uniqueN(paste(cohort, tag, protein_assay))],
        " assays) to ", out_windows)

# assay_annot.csv is assay_ensg.csv plus a uniprot_id column, written by
# 04_annot_uniprot.py. Leaving the old copy in place after rewriting assay_ensg
# lets the two disagree silently -- postgwas/ reads assay_ensg while
# projects/ read assay_annot, so they would report different assay sets. Remove
# it so any consumer fails loudly until 04 regenerates it.
annot_out <- file.path(PROTEIN_INFO_DIR, "assay_annot.csv")
if (file.exists(annot_out)) {
  if (!file.remove(annot_out)) {
    stop("Could not remove stale ", annot_out,
         " -- it would silently disagree with the assay_ensg.csv just written.")
  }
  message("Removed stale ", annot_out, " -- run 04_annot_uniprot.py to regenerate it.")
}

unmatched <- unique(summary[is.na(gene_chr), .(source, gene_raw, gene, soma_uniprot)])[order(gene_raw, source)]
unmatched[, original_name := fifelse(
  source == "somalogic" & !is.na(soma_uniprot) & nzchar(soma_uniprot),
  soma_uniprot,
  gene_raw
)]
unmatched[, soma_uniprot := NULL]
unmatched[, original_name_key := match_key(original_name)]
unmatched <- unmatched[!duplicated(original_name_key)]
unmatched <- unmatched[order(source, original_name), .(source, original_name, rename = "")]
fwrite(unmatched, out_unmatched)
message("No coords: ", nrow(unmatched), " / ", nrow(summary), ". Wrote ", out_unmatched)
if (nrow(unmatched) > 0) {
  print(unmatched)
}
