# Shared helpers for code/09_statistical_analysis. Scripts run from the repo root
# and source "code/09_statistical_analysis/lib/common.R".

# Coloc: the strong-PP4 bar and the reported-row rule are the pipeline's
# (code/lib/common.R), read here rather than restated -- four copies of the
# row rule and two of 0.8 existed before 2026-10-25.
local({
  f <- file.path(c(".", "..", "../..", "../../.."), "code/lib/common.R")
  f <- f[file.exists(f)][1L]
  if (is.na(f)) stop("code/lib/common.R not found from ", getwd())
  pg <- new.env(); sys.source(f, envir = pg)
  COLOC_PP4      <<- pg$PP4_STRONG
  coloc_reported <<- pg$coloc_reported
  FM_WINDOW_BP   <<- pg$FM_WINDOW_BP   # fine-mapping half-window around a clump lead
  GWAS_BASE      <<- pg$GWAS_BASE      # data roots
  COLOC_ROOT     <<- pg$COLOC_ROOT
  MR_ROOT        <<- pg$MR_ROOT
  LD_REF         <<- pg$LD_REF         # 1000 Genomes EUR bfile prefix
})

# Assay platform and fluid from the assay_ensg `tag` (e.g. "Neurology_o2csf").
# These two rules were copied into seven places before 2026-10-06.
assay_platform <- function(tag) data.table::fcase(
  grepl("nulisa", tolower(tag)), "NuLISA",
  grepl("somalogic", tolower(tag)), "SomaLogic",
  default = "Olink"
)
assay_fluid <- function(tag) data.table::fcase(
  grepl("csf", tolower(tag)), "CSF",
  grepl("plasma", tolower(tag)), "Plasma",
  default = NA_character_
)

# trans-pQTL hotspot: a locus acting on >= this many distinct trans-target proteins
# within one fluid (01b_hotspot_data.R defines the loci, lib/hotspot_panels.R draws them).
HOTSPOT_MIN <- 5

# APOE (hg19 chr19:45,409,011-45,412,650) +- 1 Mb; rs429358 sits at 45,411,941. A signal
# here can take an effect from e4, which drives both AD pathology and many proteins.
APOE_REGION <- list(chr = 19L, start = 44409011L, end = 46412650L)
in_apoe_region <- function(chr, pos)
  chr == APOE_REGION$chr & pos >= APOE_REGION$start & pos <= APOE_REGION$end

# GRCh37 chromosome lengths in Mb, for genome-wide x axes.
CHR_LEN_HG19 <- c(
  "1" = 249.25, "2" = 243.19, "3" = 198.02, "4" = 191.15,
  "5" = 180.92, "6" = 171.12, "7" = 159.14, "8" = 146.36,
  "9" = 141.21, "10" = 135.53, "11" = 135.01, "12" = 133.85,
  "13" = 115.17, "14" = 107.35, "15" = 102.53, "16" = 90.35,
  "17" = 81.20, "18" = 78.08, "19" = 59.13, "20" = 63.03,
  "21" = 48.13, "22" = 51.30, "X" = 155.27
)

# One assay's credible-set members from its fine-mapping output (NULL when the file is
# absent or empty; check_finemap_cs() decides whether that is allowed). Used by 03 and
# 06, which need SNP-level members that no info/ table carries.
read_finemap_cs <- function(cohort, assay, fm_dir, select = NULL) {
  f <- file.path(fm_dir, paste0(cohort, "__", assay, ".finemapped.tsv.gz"))
  if (!file.exists(f)) return(NULL)
  d <- if (is.null(select)) data.table::fread(f, showProgress = FALSE)
       else data.table::fread(f, select = select, showProgress = FALSE)
  if (!nrow(d)) return(NULL)
  d[, `:=`(cohort = cohort, protein_assay = assay)]
  d[]
}

# The per-assay outputs (read above) and info/bf_protein (read by 01a/01b/02) are two
# views of one fine-mapping run. Stop when they disagree, instead of quietly treating a
# missing output as "no pQTL":
#   - a kept assay with a clump lead has no output -> a failed or skipped run;
#   - the credible sets read differ from info/bf_protein/credible_sets.tsv.gz ->
#     fine-mapping was re-run and finemap/05 was not.
check_finemap_cs <- function(cs, kept, fm_dir,
                             lead_file = "info/bf_protein/clump_lead_pqtl.csv.gz",
                             cs_file   = "info/bf_protein/credible_sets.tsv.gz") {
  lead <- unique(data.table::fread(lead_file, select = c("cohort", "protein_assay")))
  need <- merge(unique(kept[, .(cohort, protein_assay)]), lead, by = c("cohort", "protein_assay"))
  need[, trait := paste0(cohort, "__", protein_assay)]
  miss <- need[!file.exists(file.path(fm_dir, paste0(trait, ".finemapped.tsv.gz")))]
  if (nrow(miss)) stop(nrow(miss), " kept assays with a clump lead have no fine-mapping output, e.g. ",
                       miss$trait[1], call. = FALSE)
  got <- unique(cs[, paste0(cohort, "__", protein_assay, "|", cs_id)])
  ref <- data.table::fread(cs_file, select = c("trait", "cs_id"))
  ref <- ref[trait %in% need$trait, unique(paste0(trait, "|", cs_id))]
  if (!setequal(got, ref)) stop("credible sets read (", length(got), ") differ from ", cs_file, " (",
                                length(ref), "): ", length(setdiff(got, ref)), " extra, ",
                                length(setdiff(ref, got)), " missing -- re-run finemap/05", call. = FALSE)
  invisible(TRUE)
}

# S1 with 07b's per-credible-set novelty call (07_novelty_credible_sets.csv.gz) as columns,
# in S1's final row order. Both 06 (when it rebuilds S1) and 07b (after its scan) fold it
# in, so running 06 alone keeps the columns (r019 lost them). NULL when the file is
# missing or leaves a credible-set lead without a call: the credible sets changed since
# 07b ran, so 07b has to run again.
add_novelty <- function(s1, novelty_file) {
  if (!file.exists(novelty_file)) return(NULL)
  add <- data.table::fread(novelty_file, na.strings = c("", "NA"))[, .(
    cs_key, novelty, novelty_resolution, novelty_ref_study = ref_study,
    novelty_ref_rsid = ref_best_rsid, novelty_ref_r2 = ref_best_r2,
    novelty_primary_cis)]   # primary cis association call; NA for secondary cis and trans
  s1 <- data.table::copy(s1)
  stale <- intersect(names(add)[-1], names(s1))
  if (length(stale)) s1[, (stale) := NULL]
  n_before <- nrow(s1)
  s1 <- merge(s1, add, by = "cs_key", all.x = TRUE, sort = FALSE)
  if (nrow(s1) != n_before)
    stop("merging novelty into S1 changed the row count: ", n_before, " -> ", nrow(s1))
  if (s1[is_cs_lead == TRUE & is.na(novelty), .N] > 0) return(NULL)
  data.table::setorder(s1, cohort, fluid, protein_assay, cs_id, -pip)
  s1[]
}

# S1 with 07c's external replication in Western 2026 CSF Olink (07c_western_replication.csv)
# as columns, per credible set like add_novelty and folded in by the same two scripts (06,
# 07c). Only the distinct CSF signals were tested, so plasma sets, duplicates of a distinct
# set and trans sets above the study-wide threshold stay NA. western_snp is the variant
# tested (the lead, or the next credible-set variant Western has); western_beta is aligned
# to S1's effect allele. NULL when the file is missing or names a credible set S1 lacks:
# the credible sets changed since 07c ran, so 07c has to run again.
add_western <- function(s1, western_file) {
  if (!file.exists(western_file)) return(NULL)
  add <- data.table::fread(western_file, na.strings = c("", "NA"))[, .(
    cs_key, western_status, western_snp = tested_snp, western_beta = beta_western,
    western_p = p_western, western_cs_min_p)]   # smallest Western P over the credible set
  if (anyDuplicated(add$cs_key)) stop("07c lists a credible set twice")
  if (!all(add$cs_key %in% s1$cs_key)) return(NULL)
  s1 <- data.table::copy(s1)
  stale <- intersect(names(add)[-1], names(s1))
  if (length(stale)) s1[, (stale) := NULL]
  n_before <- nrow(s1)
  s1 <- merge(s1, add, by = "cs_key", all.x = TRUE, sort = FALSE)
  if (nrow(s1) != n_before)
    stop("merging Western replication into S1 changed the row count: ", n_before, " -> ", nrow(s1))
  data.table::setorder(s1, cohort, fluid, protein_assay, cs_id, -pip)
  s1[]
}

# VEP consequence (code/03_pqtl_mapping) -> one class per variant. Used by 03b and
# 09e, so a protein-altering share means the same thing in both.
# Most severe first; a variant takes the first class any of its terms falls in.
CSQ_CLASS <- list(
  "Protein-altering" = c("stop_gained", "stop_lost", "start_lost", "frameshift_variant",
                         "inframe_insertion", "inframe_deletion", "missense_variant",
                         "protein_altering_variant", "splice_donor_variant",
                         "splice_acceptor_variant"),
  "Splice region"    = c("splice_region_variant", "splice_donor_region_variant",
                         "splice_donor_5th_base_variant", "splice_polypyrimidine_tract_variant"),
  "5' UTR"           = "5_prime_UTR_variant",
  "3' UTR"           = "3_prime_UTR_variant",
  "Synonymous"       = c("synonymous_variant", "stop_retained_variant", "start_retained_variant"),
  "Non-coding exon"  = c("non_coding_transcript_exon_variant", "mature_miRNA_variant"),
  "Intronic"         = "intron_variant",
  "Up/downstream"    = c("upstream_gene_variant", "downstream_gene_variant"),
  "Intergenic"       = c("intergenic_variant", "regulatory_region_variant", "TF_binding_site_variant")
)
csq_class <- function(csq) {
  terms <- strsplit(csq, ",", fixed = TRUE)
  vapply(terms, function(t) {
    for (k in names(CSQ_CLASS)) if (any(t %in% CSQ_CLASS[[k]])) return(k)
    "Other"
  }, "")
}

# One spelling for a UniProt ID list: "P29459;P29460", "P29459_P29460" and
# "Q96QH8;Q96QH8" become accessions deduplicated, sorted and joined by "_" (07b, 07c).
norm_ids <- function(x)
  vapply(strsplit(as.character(x), "[_;|]"),
         function(p) paste(sort(unique(p[p != ""])), collapse = "_"), character(1))

# Brain eQTL colocalization of a CSF protein's own gene, from 04's tissue table (09d,
# 09e): tested = a cognate coloc ran (status ok; no_sig_coloc means the eQTL had no
# signal, so not tested); hit = a retained row with PP4 > COLOC_PP4.
brain_coloc_sets <- function(tis_file) {
  tis <- fread(tis_file)[fluid == "CSF" & tissue_class == "Brain"]
  list(tis    = tis,
       tested = tis[status == "ok" & cognate == TRUE, unique(protein_assay)],
       hit    = tis[coloc_retained == TRUE & PP4 > COLOC_PP4, unique(protein_assay)])
}
