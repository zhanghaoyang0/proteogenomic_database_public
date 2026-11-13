#' Harmonize GWAS to reference alleles from afreq (CHR, POS, ID/SNP, REF, ALT).
#' @param impute_frq if TRUE, fill EA_FREQ from freq (REF/ALT/ALT_FREQ or ALT_FREQS) by SNP before harmonizing.
#' @param collapse_by_snp if TRUE, keep one row per SNP (lowest P).
#' @param maf_min if not NULL, keep rows with MAF = min(EA_FREQ, 1-EA_FREQ) >= maf_min (default 0.01).
#' @param pal_maf_max if not NULL, drop palindromic SNPs (A/T or C/G) with MAF > pal_maf_max
#'   (default 0.4). Mirrors post-impute QC behavior — strand of high-MAF palindromic
#'   SNPs cannot be inferred reliably from frequency alone, so flip/flip_swap decisions
#'   on them are unsafe.
# --- eQTL cleaning: shared eGene rule ----------------------------------------
# ONE threshold for every eQTL resource (metabrain, gtex, singlebrain) so the
# cleaned gene sets stay comparable across them. A gene passes if its best cis P
# beats this; passing genes keep ALL their cis SNPs, because SuSiE and coloc need
# the full regional profile -- the cut is per gene, never per SNP.
#
# 5e-8, the genome-wide threshold every other trait uses in post-GWAS (clump,
# MR instruments, and the credible-set gate in 01_finemapping.R). That gate already
# dropped every credible set without a P < 5e-8 variant, so a gene between 5e-8 and
# the former 1e-5 was fine-mapped but could never reach coloc; screening at 5e-8
# drops it up front instead. Not resource-specific, and needs no permutation/qvalue
# file (GTEx via eQTL Catalogue ships none locally). Was 1e-5 until 2026-10-04;
# clean trees built before then still hold the 1e-5 genes, so
# code/04_finemapping/prep_eqtl_cis_loci.R applies this value again.
#
# Regardless of this value, every clean script keeps a record of ALL tested genes, so
# "was this gene tested, and how strong was it?" stays answerable without re-cleaning:
# metabrain_clean.R writes its trait manifest before the filter; gtex_clean.R and
# singlebrain_clean.R write theirs after it (eGenes only) and keep the full list in
# info/<resource>/egene_audit/. That record, not the threshold, is what supports an
# absence claim ("no eQTL in bulk cortex") -- check it before writing one, since a
# gene may have been dropped at P = 1e-7 rather than truly null.
P_EGENE_DEFAULT <- 5e-8

harmonize_to_1kg <- function(gwas, freq, impute_frq = FALSE, collapse_by_snp = TRUE, maf_min = 0.01, pal_maf_max = 0.4) {
  comp <- function(x) {
    x <- toupper(x)
    out <- ifelse(is.na(x), NA_character_, chartr("ACGT","TGCA", x))
    out
  }
  data.table::setDT(gwas)
  data.table::setDT(freq)
  setnames(freq, c("#CHROM", "ID", "ALT_FREQS"), c("CHR", "SNP", "ALT_FREQ"), skip_absent = TRUE)

  if (impute_frq) {
    cols <- intersect(c("SNP", "REF", "ALT", "ALT_FREQ"), names(freq))
    n_in <- nrow(gwas)
    gwas <- merge(gwas, freq[, ..cols], by = "SNP", all.x = TRUE)
    # A left join multiplies rows on a duplicate reference key (CLAUDE.md); 1kg_v3 has none.
    if (nrow(gwas) != n_in) stop("harmonize_to_1kg: duplicate SNP IDs in the reference added ",
                                 nrow(gwas) - n_in, " rows")
    comp_i <- function(x) chartr("ACGT", "TGCA", toupper(x))
    gwas[, EA_c := comp_i(EA)]
    gwas[, EA_FREQ := data.table::fcase(
      EA == ALT, ALT_FREQ,
      EA == REF, 1 - ALT_FREQ,
      EA_c == ALT, ALT_FREQ,
      EA_c == REF, 1 - ALT_FREQ,
      default = NA_real_
    )]
    gwas[, c("EA_c", "REF", "ALT", "ALT_FREQ") := NULL]
  }

  ref <- freq[, .(CHR = as.integer(CHR), POS = as.integer(POS), SNP_bim = SNP, REF, ALT)]

  ## merge (keep SNP name from ref as SNP_bim)
  n_in <- nrow(gwas)
  d <- merge(
    gwas,
    ref,
    by = c("CHR", "POS"),
    all.x = TRUE
  )
  if (nrow(d) != n_in) stop("harmonize_to_1kg: duplicate CHR:POS in the reference added ",
                            nrow(d) - n_in, " rows")

  ## complement
  d[, `:=`(
    EA_c  = comp(EA),
    NEA_c = comp(NEA)
  )]

  ## classify: target is EA=ALT(1kg), NEA=REF(1kg)
  d[, case := data.table::fifelse(
    !is.na(REF) & !is.na(ALT) & EA == ALT & NEA == REF, "ok",
    data.table::fifelse(
      !is.na(REF) & !is.na(ALT) & EA == REF & NEA == ALT, "swap",
      data.table::fifelse(
        !is.na(REF) & !is.na(ALT) & EA_c == ALT & NEA_c == REF, "flip",
        data.table::fifelse(
          !is.na(REF) & !is.na(ALT) & EA_c == REF & NEA_c == ALT, "flip_swap",
          "unmatched"
        )
      )
    )
  )]

  ## apply fixes: set EA=ALT, NEA=REF; swap/flip_swap flip beta and freq
  d[case != "unmatched",
    `:=`(
      EA = ALT,
      NEA = REF,
      BETA = data.table::fifelse(case %in% c("swap","flip_swap"), -BETA, BETA),
      EA_FREQ = data.table::fifelse(case %in% c("swap","flip_swap"), 1 - EA_FREQ, EA_FREQ)
    )
  ]

  ## drop high-MAF palindromic SNPs (ambiguous strand). A/T or C/G pairs are their
  ## own complement, so flip/flip_swap decisions on them rely on MAF being far from
  ## 0.5 — drop those where MAF is too close to 0.5 to tell strands apart. EA_FREQ
  ## has been case-aligned to ALT above, so MAF = min(EA_FREQ, 1-EA_FREQ).
  if (!is.null(pal_maf_max) && is.finite(pal_maf_max) && pal_maf_max > 0) {
    pal <- paste(d$REF, d$ALT, sep = "/") %in% c("A/T", "T/A", "C/G", "G/C")
    maf <- pmin(d$EA_FREQ, 1 - d$EA_FREQ)
    drop <- pal & !is.na(maf) & maf > pal_maf_max
    message("palindromic SNPs with MAF > ", pal_maf_max, " dropped: ", sum(drop))
    d <- d[!drop]
  }

  ## for matched variants, assign SNP name from bim
  d[case != "unmatched" & !is.na(SNP_bim), SNP := SNP_bim]
  ## An unmatched variant that came without an ID (empty rsid in GTEx and ad_uffelmann)
  ## still has none, and nothing downstream can join on it: a GTEx gene's top eQTL like
  ## that gave 7 reported colocalizations an empty o_lead.
  no_id <- is.na(d$SNP) | !nzchar(d$SNP)
  if (any(no_id)) {
    message("variants with no ID dropped: ", sum(no_id))
    d <- d[!no_id]
  }

  matched_prop <- d[, sum(case != "unmatched") / .N]
  message("matched proportion: ", round(matched_prop, 2))

  ## only keep original GWAS columns (drop others in place to avoid copying large table)
  keep_cols <- c("CHR", "POS", "SNP", "EA", "NEA", "EA_FREQ", "BETA", "SE", "P", "N")
  drop_cols <- setdiff(names(d), keep_cols)
  if (length(drop_cols) > 0) {
    d[, (drop_cols) := NULL]
  }
  data.table::setcolorder(d, keep_cols)

  ## strict complete cases for analysis columns
  d <- d[complete.cases(EA, NEA, EA_FREQ, BETA, SE, P, N)]

  ## MAF filter (biallelic: minor allele frequency = min(p, 1-p))
  if (!is.null(maf_min) && is.finite(maf_min) && maf_min > 0) {
    n_maf <- nrow(d)
    d <- d[pmin(EA_FREQ, 1 - EA_FREQ) >= maf_min]
    message("MAF >= ", maf_min, ": ", n_maf, " \u2192 ", nrow(d), " rows")
  }

  ## optionally keep one row per SNP (lowest P)
  if (isTRUE(collapse_by_snp) && "SNP" %in% names(d)) {
    data.table::setorder(d, SNP, P)
    d <- d[!duplicated(SNP)]
  }

  before_n <- nrow(gwas)
  after_n <- nrow(d)
  message("before / after GWAS rows: ", before_n, " \u2192 ", after_n)

  return(d)
}

lift_to_hg19 <- function(dt, chain) {
  chr <- as.character(dt$CHR)
  chr <- ifelse(grepl("^chr", chr, ignore.case = TRUE), chr, paste0("chr", chr))
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = dt$POS, end = dt$POS)
  )
  lifted <- rtracklayer::liftOver(gr, chain)
  keep <- lengths(lifted) == 1L
  if (!any(keep)) return(dt[0])
  lifted_gr <- unlist(lifted[keep])
  # A variant can land on an alt or unplaced contig (chr6_cox_hap2): as.integer() made
  # its CHR NA and harmonize_to_1kg kept it, 140 rows per FinnGen file. chr1-22 only.
  lifted_chr <- sub("^chr", "", as.character(GenomeInfoDb::seqnames(lifted_gr)))
  auto <- lifted_chr %in% as.character(1:22)
  dt2 <- dt[keep][auto]
  dt2[, CHR := as.integer(lifted_chr[auto])]
  dt2[, POS := as.integer(IRanges::start(lifted_gr))[auto]]
  message("liftover unique-mapped variants rate: ", round(mean(keep), 3),
          "; landed off chr1-22, dropped: ", sum(!auto))
  return(dt2)
}

# --- Manhattan + QQ (CMplot): trait title on top; square QQ panel on the right ---

clip_p_gwas <- function(p) {
  p <- as.numeric(p)
  p[p <= 0] <- .Machine$double.xmin
  p
}

normalize_chr_gwas <- function(chr) sub("^chr", "", as.character(chr), ignore.case = TRUE)

#' Genomic inflation lambda from P (1 df): median(chi^2) / median under null.
lambda_gc_from_p <- function(p) {
  p <- pmax(pmin(as.numeric(p), 1 - 1e-300), .Machine$double.xmin)
  chisq <- stats::qchisq(1 - p, df = 1L)
  stats::median(chisq, na.rm = TRUE) / stats::qchisq(0.5, df = 1L)
}

#' Read harmonized GWAS (SNP, CHR, POS, P), build CMplot frame, filter to autosomes.
# One pass over the vectors instead of six subsetting copies of a 5.5M-row
# table: 8.0 s -> 1.6 s on a typical trait, same rows out. The regex in
# normalize_chr_gwas() is the single most expensive step, so it is skipped
# when CHR already came off fread as a number.
prep_manhattan_dt <- function(dt, chr_levels = as.character(1:22)) {
  data.table::setDT(dt)
  chr <- if (is.numeric(dt$CHR)) dt$CHR
         else suppressWarnings(as.numeric(normalize_chr_gwas(dt$CHR)))
  pos <- as.numeric(dt$POS)
  p   <- as.numeric(dt$P)
  keep <- !is.na(chr) & chr %in% as.numeric(chr_levels) &
          is.finite(pos) & pos > 0 & is.finite(p)
  data.table::data.table(
    SNP = dt$SNP[keep], Chromosome = chr[keep],
    Position = pos[keep], P = clip_p_gwas(p[keep]))
}

#' Manhattan (left) + square QQ (right); requires CMplot.
#' @param df_manhattan data.table with SNP, Chromosome, Position, P (from prep_manhattan_dt).
#' @param trait_title title above plots.
#' @param out_png output PNG path.
#' @param cis_info optional list(chr, start, end) to add a cis-zoom panel below.
#' @param zoom_height_in height in inches of the cis-zoom panel (used only if cis_info provided).
plot_gwas_manhattan_qq <- function(
    df_manhattan,
    out_manhattan,
    out_qq,
    p_gws = 5e-8,
    p_suggestive = NULL,
    width_in = 18,
    height_in = 4.5,
    qq_size_in = 5,
    dpi = 300,
    chr_levels = as.character(1:22),
    cis_info = NULL,
    zoom_height_in = 3.5,
    label_snps = NULL) {
  if (!requireNamespace("CMplot", quietly = TRUE)) {
    stop("Package 'CMplot' is required for plot_gwas_manhattan_qq()", call. = FALSE)
  }
  if (nrow(df_manhattan) == 0L) {
    message("plot_gwas_manhattan_qq: no rows, skip")
    return(invisible(FALSE))
  }
  gws_snps  <- df_manhattan[P < p_gws, SNP]
  lam       <- lambda_gc_from_p(df_manhattan$P)
  thresholds <- if (is.null(p_suggestive)) p_gws else c(p_gws, p_suggestive)
  hl_text   <- if (!is.null(label_snps) && length(gws_snps) > 0L)
    ifelse(gws_snps %in% label_snps, gws_snps, "") else NULL
  present_chrs <- as.character(sort(unique(df_manhattan$Chromosome)))

  # --- Manhattan PNG ---
  if (!is.null(out_manhattan)) {
    has_zoom     <- !is.null(cis_info)
    total_height <- if (has_zoom) height_in + zoom_height_in else height_in
    zoom_frac    <- zoom_height_in / total_height
    y_lo <- if (has_zoom) zoom_frac + 0.08 * (1 - zoom_frac) else 0.06
    y_hi <- if (has_zoom) zoom_frac + 0.94 * (1 - zoom_frac) else 0.96

    grDevices::png(out_manhattan, width = width_in, height = total_height, units = "in", res = dpi, bg = "white")
    # The bottom outer margin exists for the zoom panel's x axis; without a zoom
    # the same 2.2 lines are an empty band under the chromosome labels.
    graphics::par(oma = c(if (has_zoom) 2.2 else 0.3, 0, 0, 0))
    
    # Save device state before CMplot
    old_fig <- graphics::par("fig")
    old_mar <- graphics::par("mar")
    
    graphics::par(fig = c(0.022, 0.98, y_lo, y_hi), new = FALSE)
    # mgp[2] = 0.95, not the usual 0.55: it sets the gap between the y tick
    # numbers and the axis. CMplot draws the chromosome labels with its own
    # mgp = c(3, xticks.pos, 0), so this moves the y side only.
    graphics::par(mar = c(2.2, 2.9, 0.4, 0.5), mgp = c(1.4, 0.95, 0), xaxs = "i", yaxs = "i")
    CMplot::CMplot(df_manhattan,
      plot.type = "m", LOG10 = TRUE,
      # Same y label as the coloc and fine-mapping panels.
      ylab = expression(-log[10](italic(P))),
      ylab.pos = 2.0,
      axis.cex = 0.9,
      lab.cex = 1.1,
      lab.font = 1,     # CMplot defaults to 2; the app wants no bold anywhere
      chr.labels = present_chrs,
      threshold = thresholds,
      threshold.lty = if (length(thresholds) == 1L) 1 else c(1, 2),
      threshold.col = if (length(thresholds) == 1L) "grey60" else c("grey30", "grey60"),
      threshold.lwd = if (length(thresholds) == 1L) 2 else c(2, 1),
      col = c("#8DA0CB", "#66C2A5"),
      highlight = gws_snps,
      highlight.col = "darkred",
      highlight.cex = 0.32,
      highlight.text = hl_text,
      highlight.text.cex = 0.65,
      highlight.text.font = 3,
      amplify = FALSE,
      cex = 0.30,
      file.output = FALSE)


    if (has_zoom) {
      dt_cis <- df_manhattan[Chromosome == as.numeric(cis_info$chr)]
      if (nrow(dt_cis) > 0) {
        dt_cis <- data.table::copy(dt_cis)
        dt_cis[, logp   := -log10(P)]
        dt_cis[, pos_mb := Position / 1e6]
        cis_start_mb <- cis_info$start / 1e6
        cis_end_mb   <- cis_info$end   / 1e6
        x_range      <- c(cis_start_mb - 2, cis_end_mb + 2)
        dt_cis       <- dt_cis[pos_mb >= x_range[1] & pos_mb <= x_range[2]]
        y_max        <- if (nrow(dt_cis) > 0) max(dt_cis$logp) else 1
        gws_y        <- -log10(p_gws)
        ylim_upper   <- max(y_max * 1.08, gws_y * 1.05)
        graphics::par(fig = c(0.022, 0.98, 0.02, zoom_frac - 0.02), new = TRUE)
        graphics::par(mar = c(2.5, 2.9, 1.6, 0.5), mgp = c(1.6, 0.55, 0))
        graphics::plot(dt_cis$pos_mb, dt_cis$logp,
          type = "n",
          xlab = paste0("Chr", cis_info$chr, " position (Mb)"),
          ylab = expression(-log[10](italic(P))),
          xlim = x_range, ylim = c(0, ylim_upper),
          cex.axis = 0.85, cex.lab = 0.95,
          main = "")
        graphics::mtext("Cis-locus zoom, gene ±500 kb",
          side = 3, line = 0.3, cex = 0.95, font = 2)
        # Restore yellow shade for cis region
        graphics::rect(cis_start_mb, 0, cis_end_mb, ylim_upper,
          col = grDevices::adjustcolor("#E69F00", alpha.f = 0.12), border = NA)
        graphics::abline(h = gws_y, lty = 2, col = "grey40", lwd = 0.8)
        dt_ns  <- dt_cis[P >= p_gws]
        dt_sig <- dt_cis[P <  p_gws]
        if (nrow(dt_ns)  > 0) graphics::points(dt_ns$pos_mb,  dt_ns$logp,  pch = 20, cex = 0.25, col = "#8DA0CB")
        if (nrow(dt_sig) > 0) graphics::points(dt_sig$pos_mb, dt_sig$logp, pch = 20, cex = 0.5,  col = "darkred")
      }
    }
    grDevices::dev.off()
  }

  # --- QQ PNG (base R; CMplot opens its own device even with file.output=FALSE) ---
  if (!is.null(out_qq)) {
    p_all    <- df_manhattan$P
    n        <- length(p_all)
    p_sorted <- sort(p_all)
    obs      <- -log10(p_sorted)
    exp      <- -log10(stats::ppoints(n))
    # thin non-significant points for rendering speed
    is_sig   <- p_sorted < p_gws
    idx_ns   <- which(!is_sig)
    if (length(idx_ns) > 50000L) idx_ns <- sort(sample(idx_ns, 50000L))
    idx_plot <- sort(c(idx_ns, which(is_sig)))

    grDevices::png(out_qq, width = qq_size_in, height = qq_size_in, units = "in", res = dpi, bg = "white")
    graphics::par(mar = c(3.8, 3.8, 0.5, 0.5), mgp = c(2.2, 0.6, 0))
    # Use same xlim and ylim for square layout
    max_val <- max(exp[idx_plot], obs[idx_plot], na.rm = TRUE)
    graphics::plot(exp[idx_plot], obs[idx_plot],
      pch = 20, cex = 0.3, col = "#8DA0CB",
      xlab = expression(Expected ~ -log[10](italic(P))),
      ylab = expression(Observed ~ -log[10](italic(P))),
      xlim = c(0, max_val), ylim = c(0, max_val),
      cex.axis = 0.85, cex.lab = 1.0)
    idx_sig_plot <- which(is_sig[idx_plot])
    if (length(idx_sig_plot) > 0)
      graphics::points(exp[idx_plot][idx_sig_plot], obs[idx_plot][idx_sig_plot],
        pch = 20, cex = 0.55, col = "darkred")
    graphics::abline(0, 1, col = "grey40", lty = 2, lwd = 1.2)
    usr <- graphics::par("usr")
    graphics::text(usr[1] + 0.03*(usr[2]-usr[1]), usr[4] - 0.05*(usr[4]-usr[3]),
      bquote(lambda[GC] == .(sprintf("%.3f", lam))),
      adj = c(0, 1), cex = 1.1)
    grDevices::dev.off()
  }

  invisible(TRUE)
}

# PNG written successfully (Slurm kill after png() can leave 0-byte files; file.exists alone is not enough).
png_output_ok <- function(path) {
  if (!nzchar(path)) return(FALSE)
  if (!file.exists(path)) return(FALSE)
  sz <- file.info(path)$size
  length(sz) == 1L && !is.na(sz) && sz > 0L
}

#' Convenience: fread clean .tsv.gz, prep, plot; optional skip if both PNGs exist.
plot_manhattan_qq <- function(
    clean_path,
    out_manhattan,
    out_qq,
    skip_if_exists = TRUE,
    chr_levels = as.character(1:22),
    p_gws = 5e-8,
    p_suggestive = NULL,
    width_in = 18,
    height_in = 4.5,
    qq_size_in = 5,
    dpi = 300,
    cis_info = NULL,
    zoom_height_in = 3.5,
    label_snps = NULL) {
  for (p in c(out_manhattan, out_qq)) {
    if (!nzchar(p)) next
    if (!file.exists(p)) next
    sz <- file.info(p)$size
    if (length(sz) == 1L && !is.na(sz) && sz == 0L) unlink(p)
  }
  if (isTRUE(skip_if_exists) && png_output_ok(out_manhattan) && png_output_ok(out_qq)) {
    message("Skip plot exists: ", basename(out_manhattan))
    return(invisible(FALSE))
  }
  if (!file.exists(clean_path)) {
    message("Skip plot (no clean): ", clean_path)
    return(invisible(FALSE))
  }
  # prep_manhattan_dt uses four columns; the clean files carry ten (EA, NEA,
  # EA_FREQ, BETA, SE, N are dead weight here). Same header in bf1, bf2, finngen
  # and other, so the select is safe and halves the read.
  dt <- data.table::fread(clean_path, select = c("CHR", "POS", "SNP", "P"),
                          showProgress = FALSE)
  dfm <- prep_manhattan_dt(dt, chr_levels = chr_levels)
  if (nrow(dfm) == 0L) {
    message("Skip plot (no rows): ", clean_path)
    return(invisible(FALSE))
  }
  plot_gwas_manhattan_qq(
    dfm,
    out_manhattan = if (isTRUE(skip_if_exists) && png_output_ok(out_manhattan)) NULL else out_manhattan,
    out_qq        = if (isTRUE(skip_if_exists) && png_output_ok(out_qq))        NULL else out_qq,
    p_gws = p_gws, p_suggestive = p_suggestive,
    width_in = width_in, height_in = height_in, qq_size_in = qq_size_in, dpi = dpi,
    chr_levels = chr_levels,
    cis_info = cis_info, zoom_height_in = zoom_height_in,
    label_snps = label_snps)
}

#' Write a tab-separated table through a temp file and a rename, so a killed job
#' never leaves a half-written file that a skip-if-exists check would accept.
#' Gzipped when `path` ends in .gz.
write_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  data.table::fwrite(x, tmp, sep = "\t", compress = if (grepl("[.]gz$", path)) "gzip" else "none")
  if (!file.rename(tmp, path)) stop("could not write ", path)
  invisible(path)
}

#' Round harmonized GWAS columns to compact significant figures, and stop on a row with
#' no CHR, POS or SNP. Apply after harmonize_to_1kg() before fwrite(): every clean
#' script ends here, so this is the one check that no such row reaches a clean file.
round_gwas_clean <- function(d) {
  no_key <- is.na(d$CHR) | is.na(d$POS) | is.na(d$SNP) | !nzchar(d$SNP)
  if (any(no_key))
    stop(sum(no_key), " clean rows without CHR, POS or SNP; lift_to_hg19 and ",
         "harmonize_to_1kg should have dropped them")
  d[, BETA    := signif(as.numeric(BETA),    3)]
  d[, SE      := signif(as.numeric(SE),      4)]
  d[, EA_FREQ := signif(as.numeric(EA_FREQ), 3)]
  d[, P       := signif(as.numeric(P),       4)]
  n_before <- nrow(d)
  d <- d[!is.na(BETA) & !is.na(SE) & !is.na(P)]
  n_drop <- n_before - nrow(d)
  if (n_drop > 0) message("round_gwas_clean: dropped ", n_drop, " rows with NA in BETA/SE/P")
  d
}


#' BF1 protein panels that have GWAS output, with their assay counts, from the
#' "gwas (sample x protein)" column of info/bf_protein/qc_panel_summary.csv (written
#' by code/02_protein_qc/02 with the GWAS input). The BF1 clean and its QC
#' (gwas/submit_bf1_pqtl_clean.sh, 03_bf1_sumstats_qc.R) expect exactly these.
bf1_gwas_panels <- function(path = "info/bf_protein/qc_panel_summary.csv") {
  s <- data.table::fread(path)
  s <- s[cohort == "bf1"]
  n <- as.integer(sub(".* x ", "", s[["gwas (sample x protein)"]]))
  if (anyNA(n)) stop("unparsed gwas column in ", path)
  out <- data.table::data.table(panel = s$panel, expected_assays = n)[expected_assays > 0]
  if (!nrow(out)) stop("no BF1 GWAS panels in ", path)
  out
}
