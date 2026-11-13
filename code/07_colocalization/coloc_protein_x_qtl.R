#!/usr/bin/env Rscript
# Shared body of the three "protein x eQTL" coloc runners (GTEx, SingleBrain,
# MetaBrain). All three are COGNATE-scoped: a protein's cis-pQTL against the eQTL
# of ITS OWN gene only, across tissues / cell types. (MetaBrain was overlap-scoped
# -- every credible set against every gene in range -- until 2026-10-01; the
# switches for that were removed on 2026-10-06.)
#
# The "x disease" chains keep their own loop (coloc_x_disease.R): they iterate
# (exposure, outcome) pairs and write one chunk file per row range, whereas
# these iterate proteins and append to staging/. What is kept from each comparison
# is shared (coloc_acc_* in functions_coloc.R).
#
# coloc.susie only: both sides must be fine-mapped (>= 1 credible set). Per
# protein locus x overlapping QTL locus, run coloc.susie on the two saved SuSiE
# fits, aligned by rsid. No coloc.abf fallback -- a status row is recorded
# instead. Rows with PP4 >= PP4_KEEP are "ok"; a protein tested only below it gets
# one "no_sig_coloc" row carrying its highest-PP4 comparison (since 2026-10-24, as the
# x_disease chains have since 2026-10-23).
#
# cfg:
#   coloc_dir     output/input root under ~/db/shared_data/coloc
#   qtl_info      (tissue, ENSG) manifest csv with CHR/start/end/GENE
#   qtl_label     name used in the manifest log line
#   qtl_finemap   function(tissue, ensg) -> finemap tsv path
#   qtl_rds       function(tissue, ensg, lead) -> saved SuSiE fit path
#   no_cs_status  status string when no locus pair had a CS on both sides
run_coloc_protein_x_qtl <- function(cfg) {
  # Every row is stamped with the job's start; 01_pair redoes a protein whose fine-maps
  # are newer (db_stale_keys in functions_db.R).
  RUN_STARTED <- as.integer(Sys.time())
  COLOC_DIR <- path.expand(cfg$coloc_dir)
  # Overridable so a run can be compared against another build without touching
  # the real staging directory.
  STAGE_DIR <- Sys.getenv("COLOC_STAGE_DIR", file.path(COLOC_DIR, "staging"))
  # COLOC_TODO: this submission's own copy (submit_coloc_chain.sh).
  PROT_LIST <- Sys.getenv("COLOC_TODO", Sys.getenv("COLOC_PROT_LIST", file.path(COLOC_DIR, "info", "protein_list_todo.csv")))
  ASSAY_ENSG <- path.expand("./info/bf_protein/ensg/assay_ensg.csv")

  # Fixed column set for the staging tsv. Status-only rows are padded with NA so
  # the appended file stays rectangular; the ledger inherits this schema.
  STAGE_COLS <- c("exposure_id", "status", "cs_id", "coloc_method",
                  "tissue", "ENSG", "GENE", "lead_snp", "o_lead", "chr", "pos",
                  "n_snps", "e_hit", "o_hit", "e_idx", "o_idx",
                  "PP0", "PP1", "PP2", "PP3", "PP4", "computed_at")

  ROW_START <- as.integer(Sys.getenv("PROT_ROW_START", "1"))
  ROW_END   <- as.integer(Sys.getenv("PROT_ROW_END",   "0"))
  if (ROW_END < ROW_START) ROW_END <- ROW_START
  STAGE_FILE <- stage_file_path(STAGE_DIR, fallback = sprintf("rows%s_%s", ROW_START, ROW_END))

  t0  <- proc.time()[[3L]]
  msg <- function(...) message(sprintf("[%.1fs] ", proc.time()[[3L]] - t0), ...)

  prot <- fread(PROT_LIST)
  if (ROW_START > nrow(prot)) stop("rows ", ROW_START, "-", ROW_END, " are past the ", nrow(prot), " rows of ", PROT_LIST)
  prot <- prot[ROW_START:min(ROW_END, nrow(prot))]
  if (anyNA(prot$exposure_id)) stop("NA exposure_id in ", PROT_LIST)
  msg("rows ", ROW_START, "-", ROW_END, " -> ", nrow(prot), " proteins")
  msg("staging -> ", STAGE_FILE)

  # QTL candidate prefilter: (tissue, ENSG) cis range (gene +-window) by chr.
  qtl_manifest <- fread(path.expand(cfg$qtl_info))
  setnames(qtl_manifest, c("REGION", "ENSG_id"), c("tissue", "ENSG"))
  qtl_manifest[, CHR := as.character(CHR)]
  msg(cfg$qtl_label, " manifest: ", nrow(qtl_manifest), " (tissue,ENSG) entries")

  # Cognate map: exposure_id (<cohort>__<protein_assay>) -> the protein's own
  # gene ENSG(s).
  ae <- fread(ASSAY_ENSG)
  ae[, exposure_id := paste0(cohort, "__", protein_assay)]
  # ensembl_gene_id can list several genes ("ENSG1,ENSG2"): one row per gene.
  ae <- ae[!is.na(ensembl_gene_id), .(ENSG = sub("\\..*$", "", trimws(unlist(strsplit(ensembl_gene_id, ","))))),
           by = exposure_id]
  cog_map <- ae[nzchar(ENSG), .(ENSG = unique(ENSG)), by = exposure_id]
  setkey(cog_map, exposure_id)
  msg("cognate map: ", uniqueN(cog_map$exposure_id), " proteins -> gene(s)")

  # Protein finemap loci (clump_lead_snp + chr/pos). Keyed on the clump lead to
  # locate the saved protein SuSiE fit.
  ereg_cache <- list()
  get_ereg <- function(id, src) {
    if (is.null(ereg_cache[[id]]))
      ereg_cache[[id]] <<- finemap_loci(finemap_path(src, id), cis_only = TRUE)
    ereg_cache[[id]]
  }
  # QTL finemap loci per (tissue, ENSG) (one locus each); cached.
  qtl_loci_cache <- list()
  get_qtl_loci <- function(tissue, ensg) {
    key <- paste0(tissue, "_", ensg)
    if (is.null(qtl_loci_cache[[key]]))
      qtl_loci_cache[[key]] <<- finemap_loci(cfg$qtl_finemap(tissue, ensg), cis_only = FALSE)
    qtl_loci_cache[[key]]
  }

  n_keep <- 0L
  for (i in seq_len(nrow(prot))) {
    p <- prot[i]
    cogs <- cog_map[.(p$exposure_id), ENSG, nomatch = NULL]
    if (!length(cogs)) {
      stage_append(STAGE_FILE, data.table(exposure_id = p$exposure_id, status = "no_cognate")[, computed_at := RUN_STARTED], STAGE_COLS)
      next
    }
    regs <- get_ereg(p$exposure_id, p$exposure_source)
    if (is.null(regs) || !nrow(regs)) {
      stage_append(STAGE_FILE, data.table(exposure_id = p$exposure_id, status = "no_region")[, computed_at := RUN_STARTED], STAGE_COLS)
      next
    }

    acc <- coloc_acc()
    for (k in seq_len(nrow(regs))) {
      r     <- regs[k]
      rds_e <- protein_susie_path(p$exposure_source, p$exposure_id, r$clump_lead_snp)
      # The locus has credible-set rows, so its SuSiE fit must exist; skipping it
      # would record the pair as done with no coloc.
      if (!file.exists(rds_e)) stop("finemapped locus without its SuSiE fit: ", rds_e)

      # Candidates: (tissue, ENSG) whose cis range overlaps the protein locus span +-
      # window, restricted to the protein's own gene(s).
      cand <- qtl_manifest[CHR == as.character(r$chr) &
                           start <= r$hi + CIS_HALF_BP &
                           end   >= r$lo - CIS_HALF_BP]
      cand <- cand[ENSG %in% cogs]
      if (!nrow(cand)) next

      # Exposure fit is constant across the candidate loop below -- load once per region.
      fe_k <- rekey_susie(readRDS(rds_e))
      if (is.null(fe_k)) next

      for (j in seq_len(nrow(cand))) {
        ol <- get_qtl_loci(cand$tissue[j], cand$ENSG[j])
        if (is.null(ol) || !nrow(ol)) next
        ol <- loci_near(ol, r)
        if (!nrow(ol)) next
        for (m in seq_len(nrow(ol))) {
          rds_o <- cfg$qtl_rds(cand$tissue[j], cand$ENSG[j], ol$clump_lead_snp[m])
          sres  <- run_coloc_susie(rds_e, rds_o, fe = fe_k)
          annot <- function(d, status)
            d[, `:=`(exposure_id = p$exposure_id, cs_id = k, status = status,
                     coloc_method = "susie", lead_snp = r$clump_lead_snp,
                     o_lead = ol$clump_lead_snp[m], chr = r$chr, pos = r$pos,
                     tissue = cand$tissue[j], ENSG = cand$ENSG[j], GENE = cand$GENE[j])]
          n_keep <- n_keep + coloc_acc_add(acc, sres, annot)
        }
      }
    }

    out <- coloc_acc_rows(acc, cfg$no_cs_status,
                          function(st) data.table(exposure_id = p$exposure_id, status = st),
                          unit = p$exposure_id)
    stage_append(STAGE_FILE, out[, computed_at := RUN_STARTED], STAGE_COLS)
    msg("done ", i, "/", nrow(prot), " ", p$exposure_id, " (rows=", nrow(out), ", PP4>",
        PP4_KEEP, " rows kept so far: ", n_keep, ")")
  }
  coloc_failed_report()
}
