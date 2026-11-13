# Locus-style figures: the whole toolkit, shared by every plotter that draws one.
#
# Three families use it -- the fine-mapping regional plot, the fine-mapping PIP
# plot (code/04_finemapping/plot_protein_loci.R) and the coloc two-panel plot
# (code/07_colocalization/plot_coloc_loci.R) -- and they now share the drawing code
# itself, not just a table of numbers.
#
# It used to be two implementations in two languages: the fine-mapping plots were
# Python/matplotlib and the coloc plots R/ggplot2, so a palette, an LD reader, a
# gene track, a legend and a layout existed twice. Keeping them consistent meant
# keeping code/lib/locus_plot_spec.json in sync by hand, and they still drifted
# (the "no LD" grey was #c7cbd0 on one side and #bdbdbd on the other). Python was
# the side that could move: it only needed .tsv.gz inputs and shelled out to plink,
# while the coloc plotter reads SuSiE fits from .rds, which R alone can open.
#
# Numbers still come from locus_plot_spec.json. It is now one language's constants
# file rather than a bridge between two, but it stays a data file: it is the one
# place to change a colour or a gap without reading code.
#
# Usage:
#   source(file.path(REPO, "code/lib/locus_plot.R"))
#   SP <- locus_spec()
#   ld_panel_init("~/db/shared_data/ref_geno/1kg_v3/EUR")
#   gene_bed_init("~/db/shared_data/coloc/info/ensembl87_grch37_genes.bed.gz")
#   g   <- gene_track_rows(chr, from, to)
#   lay <- layout_heights(SP$order$locus, gene_rows = max(g$row))
#   fig <- stack_plots(list(ld_legend_plot(...),
#                           assoc_panel(...), gene_panel(g, ...)),
#                      lay$drawn + lay$gap_after)
#   ggsave(f, fig, width = lay$width, height = lay$total, dpi = SP$render$dpi)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})

## ---------------------------------------------------------------- the spec

locus_spec <- function(path = Sys.getenv("BF_LOCUS_PLOT_SPEC", "")) {
  if (!nzchar(path)) {
    here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
              error = function(e) ".")
    path <- normalizePath(file.path(here, "locus_plot_spec.json"),
                          mustWork = FALSE)
    if (!file.exists(path))
      path <- path.expand("code/lib/locus_plot_spec.json")
  }
  if (!file.exists(path)) stop("Missing locus plot spec: ", path)
  fromJSON(path, simplifyVector = TRUE)
}

# ggplot's size= is roughly a diameter in mm; the spec stores points.
pt_to_mm <- function(pt) as.numeric(pt) * 25.4 / 72

ld_palette <- function(sp) setNames(c(sp$ld$colours, sp$ld$na_colour),
                                    c(sp$ld$bins, sp$ld$na_label))

# An ordered factor with fixed levels, so a window that uses three bins still
# shows all six keys.
ld_bin <- function(r2, sp) {
  lv <- c(sp$ld$bins, sp$ld$na_label)
  cut_lab <- cut(r2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                 labels = sp$ld$bins, right = FALSE)
  factor(ifelse(is.na(r2), sp$ld$na_label, as.character(cut_lab)), levels = lv)
}

base_theme <- function(sp) {
  theme_bw(base_size = sp$type$axis_title_pt) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = sp$type$panel_title_pt,
                                    colour = "#2c3e50", margin = margin(b = 3)),
          axis.text = element_text(size = sp$type$tick_pt),
          legend.position = "none")
}

## ---------------------------------------------------------------- LD reference
# The plink1 .bed is variant-major and the .bim is position-sorted, so a window is
# one contiguous block: a single seek + read. The .bim is 10M rows, so it is read
# once per job, not once per figure. The fine-mapping plots used to shell out to
# plink --r2 once per locus, i.e. 25,837 process launches per full render.

.LD <- new.env(parent = emptyenv())

ld_panel_init <- function(bfile, quiet = FALSE) {
  .LD$bfile <- path.expand(bfile)
  .LD$bim <- fread(paste0(.LD$bfile, ".bim"), header = FALSE,
                   col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))
  .LD$bim[, CHR := as.character(CHR)]
  .LD$nind <- nrow(fread(paste0(.LD$bfile, ".fam"), header = FALSE, select = 1L))
  if (!quiet) message(format(Sys.time(), "%H:%M:%S "), "LD panel: ",
                      nrow(.LD$bim), " variants x ", .LD$nind, " individuals")
  invisible(TRUE)
}

ld_r2 <- function(chr, from, to, ref_snp, ref_pos = NA_real_) {
  if (is.null(.LD$bim)) stop("call ld_panel_init() first")
  idx <- .LD$bim[, which(CHR == as.character(chr) & POS >= from & POS <= to)]
  if (!length(idx)) return(NULL)
  sub <- .LD$bim[idx]
  ref_i <- match(ref_snp, sub$SNP)
  if (is.na(ref_i)) {
    if (is.na(ref_pos)) return(NULL)
    ref_i <- which.min(abs(sub$POS - ref_pos))
    # Said out loud: every r2 in the figure is now to a stand-in variant.
    message("  ld_r2: ", ref_snp, " not in the LD panel; r2 is to the nearest panel variant ",
            sub$SNP[ref_i], " (", abs(sub$POS[ref_i] - ref_pos), " bp away)")
  }

  nb <- ceiling(.LD$nind / 4)
  con <- file(paste0(.LD$bfile, ".bed"), "rb")
  on.exit(close(con))
  if (!identical(as.integer(readBin(con, "raw", 3L)), c(108L, 27L, 1L)))
    stop("not a variant-major plink1 .bed")
  seek(con, 3 + (min(idx) - 1) * nb)
  raw  <- readBin(con, "raw", n = length(idx) * nb)
  bits <- as.integer(rawToBits(raw))                      # LSB first, 8 bits per byte
  code <- matrix(bits[c(TRUE, FALSE)] + 2L * bits[c(FALSE, TRUE)], nrow = 4 * nb)
  g <- code[seq_len(.LD$nind), , drop = FALSE]            # variants in columns
  g[] <- c(2L, NA, 1L, 0L)[g + 1L]                        # 00=2 01=NA 10=1 11=0
  # Correlation with the reference column as one mat-vec: 1kg has few missing
  # calls, so mean-impute rather than looping cor(use = "pairwise") per variant.
  x  <- g[, ref_i]
  ok <- !is.na(x)
  gg <- g[ok, , drop = FALSE]; xx <- as.numeric(x[ok])
  na <- which(is.na(gg), arr.ind = TRUE)
  if (nrow(na)) gg[na] <- colMeans(gg, na.rm = TRUE)[na[, 2L]]
  gg <- sweep(gg, 2L, colMeans(gg))
  xc <- xx - mean(xx)
  den <- sqrt(colSums(gg^2) * sum(xc^2))
  r2  <- ifelse(den > 0, as.vector(crossprod(gg, xc)) / den, NA_real_)^2
  data.table(POS = sub$POS, SNP = sub$SNP, r2 = r2)
}

## ---------------------------------------------------------------- locus window
# One discovery locus's variants with r2 to the clump lead: what the regional PNG
# draws and what the app's interactive plot reads, built in one place so the two
# cannot disagree. `ss` is the trait's clean sumstats, `sig` one row of
# signals.tsv.gz. The window is signals.tsv.gz's own -- exactly the window SuSiE
# ran in, never recomputed from a credible set's position, which used to push CS
# variants at the far edge outside the drawn range.
locus_window <- function(ss, sig, sp, p_min_vis = 1e-100) {
  if (is.null(.LD$bim)) stop("call ld_panel_init() first")
  from <- as.integer(sig$window_start); to <- as.integer(sig$window_end)
  d <- ss[CHR == as.character(sig$chr) & POS >= from & POS <= to]
  if (nrow(d) < 10L) { message("  ", sig$signal_id, ": too few variants, skipping"); return(NULL) }
  ld <- ld_r2(sig$chr, from, to, sig$clump_lead_snp, sig$pos)
  n_in <- nrow(d)
  d <- if (is.null(ld)) copy(d)[, r2 := NA_real_] else
    merge(d, ld[, .(POS, r2)], by = "POS", all.x = TRUE)
  # merge() on POS multiplies rows if the LD panel has two variants at one position
  # (none on the 1kg_v3 EUR panel today); a longer table would draw such SNPs twice.
  if (nrow(d) != n_in)
    stop(sig$signal_id, ": window has ", nrow(d), " rows from ", n_in, " variants")
  d[SNP == sig$clump_lead_snp, r2 := 1]
  d <- d[P >= 0 & !is.na(POS)]
  if (!nrow(d)) return(NULL)
  d[, mlog10p := -log10(pmax(P, p_min_vis))]
  d[, bin := ld_bin(r2, sp)]
  list(signal_id = as.character(sig$signal_id), chr = as.character(sig$chr),
       from = from, to = to, lead_pos = as.integer(sig$pos), d = d,
       cognate = sub("\\..*$", "", as.character(sig$ensembl_gene_id)))
}

# r2 to 3 decimals for the exported tables, rounded DOWN: the LD bins break at
# 0.2/0.4/0.6/0.8 with the lower edge closed, so round() moved 0.7996 up into
# >= 0.8 and the app drew a point in a different colour from the PNG. floor()
# never crosses a bin edge. The 1e-9 absorbs 0.6 * 1000 = 599.99999...
r2_floor3 <- function(r2) floor(r2 * 1000 + 1e-9) / 1000

# An exported file counts as done only if no input is newer than it. Existence
# alone kept old windows and PIPs after fine-mapping or signals.tsv.gz was rerun.
# The submit scripts' bash pre-filters apply the same rule (`[[ in -nt out ]]`).
# A change to the LD panel, the chain or the code is not an input: use FORCE=1.
up_to_date <- function(out, inputs) {
  file.exists(out) && isTRUE(file.size(out) > 0) &&
    all(file.exists(inputs)) && all(file.mtime(inputs) <= file.mtime(out))
}

## ---------------------------------------------------------------- liftover
# GRCh37 -> GRCh38 for positions on chromosomes `chr` (no "chr" prefix). NA
# unless the point lifts to exactly one place on the same chromosome -- the rule
# of code/06_external_data/build_pos_b38.R. The chain is read once per process; rtracklayer
# is loaded only here, so plotters that never lift do not need it.

.LIFT <- new.env(parent = emptyenv())

lift_b38 <- function(chr, pos, chain_gz = Sys.getenv("CHAIN_HG19_TO_HG38",
                       "~/db/shared_data/ref_geno/liftover/hg19ToHg38.over.chain.gz")) {
  if (is.null(.LIFT$chain)) {
    suppressPackageStartupMessages(library(rtracklayer))
    chain_gz <- path.expand(chain_gz)
    if (!file.exists(chain_gz)) stop("Missing chain: ", chain_gz)
    txt <- tempfile(fileext = ".chain")          # import.chain reads plain text only
    R.utils::gunzip(chain_gz, txt, remove = FALSE)
    .LIFT$chain <- import.chain(txt)
  }
  out <- rep(NA_integer_, length(pos))
  ok <- !is.na(pos)
  if (!any(ok)) return(out)
  hit <- liftOver(GRanges(paste0("chr", chr[ok]), IRanges(pos[ok], width = 1L)), .LIFT$chain)
  n <- lengths(hit)
  one <- unlist(hit[n == 1L])
  same <- as.character(seqnames(one)) == paste0("chr", chr[ok][n == 1L])
  out[which(ok)[which(n == 1L)[same]]] <- start(one)[same]
  out
}

## ---------------------------------------------------------------- gene track

.GENES <- new.env(parent = emptyenv())

gene_bed_init <- function(bed) {
  bed <- path.expand(bed)
  # The build is part of the cache filename so a GRCh38 BED cannot be picked up
  # silently: it was, once, and every gene label sat ~1.9 Mb off its true
  # position while the track named the wrong genes.
  if (!grepl("grch37", basename(bed), ignore.case = TRUE))
    stop("gene BED must be a GRCh37 cache, got: ", basename(bed))
  if (!file.exists(bed)) stop("Missing gene BED: ", bed)
  g <- fread(cmd = paste("zcat", shQuote(bed)), header = FALSE,
             col.names = c("CHR", "START", "END", "STRAND", "ENSG", "SYMBOL", "TYPE"))
  g <- g[TYPE == "protein_coding" & !grepl("^ENSG", SYMBOL)]
  g[, CHR := as.character(CHR)]
  .GENES$g <- g
  invisible(TRUE)
}

# Genes overlapping the window, packed greedily into rows so labels do not collide.
gene_track_rows <- function(chr, from, to, sp = locus_spec()) {
  if (is.null(.GENES$g)) stop("call gene_bed_init() first")
  g <- .GENES$g[CHR == as.character(chr) & END >= from & START <= to]
  if (!nrow(g)) return(g)
  g <- copy(g); setorder(g, START)
  g[, `:=`(x0 = pmax(START, from), x1 = pmin(END, to))]   # clamp: a gene running
  mid  <- (g$x0 + g$x1) / 2                               # past the edge still draws
  edge <- 0.08 * (to - from)                              # names this close to the
  g[, lab_x := pmin(pmax(mid, from), to)]                 # edge are anchored to it
  g[, lab_h := fifelse(mid > to - edge, 1, fifelse(mid < from + edge, 0, 0.5))]
  row_end <- numeric(0); g[, row := 0L]
  pad <- sp$spacing_in$gene_label_pad_frac * (to - from)  # room for the italic label
  for (i in seq_len(nrow(g))) {
    k <- which(row_end < g$START[i] - pad)
    k <- if (length(k)) k[1] else length(row_end) + 1L
    g$row[i] <- k; row_end[k] <- g$END[i]
  }
  g
}

# `cognate` is the ENSG (or ENSGs) the figure is about. Matched by ENSG, never by
# symbol -- 144 Olink assays carry a legacy symbol that no longer matches GENCODE.
# One chip style even when there are two: the protein's gene and the eQTL's gene
# are the same in 10,312 of 10,331 colocalizations, and in the 19 exceptions they
# are neighbours, so two chips read as "two genes in play".
gene_panel <- function(g, chr, from, to, sp, cognate = character()) {
  if (!nrow(g)) return(NULL)
  g <- copy(g)
  g[, is_cognate := sub("\\..*$", "", as.character(ENSG)) %in% cognate]
  ggplot(g) +
    geom_segment(aes(x = x0 / 1e6, xend = x1 / 1e6, y = row, yend = row),
                 linewidth = sp$marks$gene_line_width_pt * .35,
                 colour = sp$marks$gene_line_colour, lineend = "butt") +
    geom_label(data = g[is_cognate == TRUE],
               aes(x = lab_x / 1e6, y = row + .40, label = SYMBOL, hjust = lab_h),
               size = sp$type$gene_label_pt / .pt, fontface = "italic",
               colour = "grey20", fill = sp$cognate$chip_fill,
               linewidth = sp$cognate$chip_edge_width_pt * .35,
               label.r = unit(1.2, "pt"), label.padding = unit(1.1, "pt")) +
    geom_text(data = g[is_cognate == FALSE],
              aes(x = lab_x / 1e6, y = row + .40, label = SYMBOL, hjust = lab_h),
              size = sp$type$gene_label_pt / .pt, fontface = "italic", colour = "grey20") +
    scale_x_continuous(limits = c(from, to) / 1e6, expand = c(0.01, 0)) +
    scale_y_reverse(expand = expansion(add = .7)) +
    labs(x = sub("\\{chr\\}", chr, sp$labels$xlabel)) +
    theme_void(base_size = sp$type$axis_title_pt) +
    theme(axis.title.x = element_text(size = sp$type$axis_title_pt, margin = margin(t = 4)),
          axis.text.x = element_text(size = sp$type$tick_pt),
          axis.ticks.x = element_line(linewidth = .3),
          axis.line.x = element_line(linewidth = .3),
          plot.background = element_rect(fill = "white", colour = NA))
}

## ---------------------------------------------------------------- data panels

# -log10(P) over a window, coloured by LD to whatever variant the caller anchored
# on. `d` needs POS, mlog10p and bin; `lead` is a POS.
assoc_panel <- function(d, lead, title, from, to, sp, cs_pos = integer(),
                        label_lead = TRUE, x_axis = FALSE) {
  d <- d[order(!is.na(r2), r2)]                           # strong LD drawn last
  lab <- d[POS == lead][1]
  p <- ggplot(d, aes(POS / 1e6, mlog10p)) +
    geom_hline(yintercept = -log10(sp$marks$sig_line), linetype = 3,
               colour = "grey55", linewidth = .3) +
    geom_point(aes(fill = bin), shape = 21, colour = "grey30", stroke = .1,
               size = pt_to_mm(sp$marks$point_diameter_pt)) +
    scale_fill_manual(values = ld_palette(sp), drop = FALSE, guide = "none") +
    scale_x_continuous(limits = c(from, to) / 1e6, expand = c(0.01, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = title, y = expression(-log[10](italic(P)))) + base_theme(sp)
  p <- p + if (x_axis) labs(x = sub("\\{chr\\}", d$CHR[1], sp$labels$xlabel))
           else theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                      axis.ticks.x = element_blank())
  if (length(cs_pos)) p <- p +
    geom_point(data = d[POS %in% cs_pos], shape = 21,
               size = pt_to_mm(sp$marks$point_diameter_pt_sparse), fill = NA,
               colour = sp$marks$cs_ring_colour, stroke = sp$marks$cs_ring_width_pt)
  if (nrow(lab) && !is.na(lab$POS)) {
    p <- p + geom_point(data = lab, shape = 23,
                        size = pt_to_mm(sp$marks$top_variant_diameter_pt),
                        fill = sp$marks$top_variant_colour, colour = "black", stroke = .2)
    if (label_lead) p <- p +
      geom_repel_label(lab, sp, max(d$mlog10p))
  }
  p
}

geom_repel_label <- function(lab, sp, ymax) {
  geom_text_repel(data = lab, aes(label = SNP), size = sp$type$tick_pt / .pt,
                  nudge_y = .06 * ymax, min.segment.length = 0,
                  box.padding = .45, seed = 1)
}

# PIP over the same window. Credible-set variants only: the association-only rows
# merge in at PIP 0 and would draw an uninformative carpet.
pip_panel <- function(d, top, title, from, to, sp) {
  d <- d[order(!is.na(r2), r2)]
  lab <- d[POS == top][1]
  ggplot(d, aes(POS / 1e6, pip)) +
    geom_linerange(aes(ymin = 0, ymax = pip), colour = "#a9afb7", linewidth = .25) +
    geom_point(aes(fill = bin), shape = 21, colour = "grey30", stroke = .1,
               size = pt_to_mm(sp$marks$point_diameter_pt_sparse)) +
    scale_fill_manual(values = ld_palette(sp), drop = FALSE, guide = "none") +
    scale_x_continuous(limits = c(from, to) / 1e6, expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(-0.055, 1.05), breaks = c(0, 0.5, 1)) +
    geom_point(data = lab, shape = 23, size = pt_to_mm(sp$marks$top_variant_diameter_pt),
               fill = sp$marks$top_variant_colour, colour = "black", stroke = .2) +
    labs(title = title, y = "PIP", x = sub("\\{chr\\}", d$CHR[1], sp$labels$xlabel)) +
    base_theme(sp)
}

## ---------------------------------------------------------------- furniture

# The LD key as its own strip, so its height is exactly spacing_in$legend_row and
# it sits the same distance above the data in every family. A guide collected by
# ggplot has no height we can pin, and pinning it is the entire point.
# `stat` rides at the far right of the same row. It has no strip of its own: the
# bold header row that used to carry it was removed because the page states the
# same facts above the image, but a statistic is not a restatement -- it is the
# result, and it should travel with an exported PNG.
ld_legend_plot <- function(sp, label, extra_key = NULL, stat = NULL) {
  keys <- c(rev(sp$ld$bins), sp$ld$na_label)
  cols <- c(rev(sp$ld$colours), sp$ld$na_colour)
  ok <- function(x) length(x) == 1L && !is.na(x) && nzchar(x)
  has_key <- ok(extra_key); has_stat <- ok(stat)
  # The bins stop short of whatever follows them on the right, and the label on
  # their left shrinks the run further, so the row is laid out from both ends:
  # six keys crammed into half the width overlap their own dots.
  left  <- 0.235
  right <- if (has_key && has_stat) 0.66 else if (has_key) 0.78 else if (has_stat) 0.84 else 0.96
  x <- seq(left + 0.05, right, length.out = length(keys))
  d <- data.frame(x = x, y = 0.5, key = keys, col = cols, stringsAsFactors = FALSE)
  p <- ggplot(d) +
    geom_point(aes(x, y), colour = d$col, size = pt_to_mm(sp$marks$point_diameter_pt_sparse)) +
    geom_text(aes(x + 0.012, y, label = key), hjust = 0,
              size = sp$type$legend_pt / .pt, colour = "#111111") +
    annotate("text", x = left, y = 0.5, label = label, hjust = 1,
             size = sp$type$legend_pt / .pt, colour = "#111111") +
    scale_x_continuous(limits = c(0, 1.06), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))
  if (has_key)
    p <- p + annotate("text", x = if (has_stat) 0.945 else 1.04, y = 0.5,
                      label = extra_key, hjust = 1,
                      size = sp$type$legend_pt / .pt, colour = sp$marks$top_variant_colour)
  if (has_stat)
    p <- p + annotate("text", x = 1.06, y = 0.5, label = stat, hjust = 1,
                      size = sp$type$legend_pt / .pt, colour = "#5a6b7d")
  p
}

## ---------------------------------------------------------------- layout

# cowplot, not patchwork: patchwork's absolute unit heights are PANEL heights, so
# a theme_void() strip (title, LD key, gene track) collapses to nothing and the
# bottom axis is pushed off the canvas. plot_grid's rel_heights apply to the whole
# sub-plot, and align="v"/axis="lr" still lines the panels up.
# Note also that Reduce(`/`, parts) builds a left-nested binary tree, so
# plot_layout(heights=) would reach only the outermost two children.
stack_plots <- function(parts, heights) {
  plot_grid(plotlist = parts, ncol = 1, rel_heights = heights,
            align = "v", axis = "lr") +
    theme(plot.background = element_rect(fill = "white", colour = NA))
}

# Inches per stacked element, in the spec's order. Every gap is an absolute
# distance: sizing by fraction of figure height is what once put the legend
# 0.63 in below the axis on a short figure and 0.37 in below on a tall one.
# `drawn` and `gap_after` are what stack_plots() wants; gaps ride as bottom
# margins on the element above.
layout_heights <- function(order, gene_rows = 0, sp = locus_spec()) {
  s <- sp$spacing_in
  gene_h <- max(s$gene_min, gene_rows * s$gene_row)
  nm <- c(); hh <- c()
  add <- function(name, h) { nm <<- c(nm, name); hh <<- c(hh, h) }
  for (seg in order) {
    if (seg == "title") add("title", s$title_row)
    else if (seg == "legend") { add("legend", s$legend_row); add("_gap", s$gap_legend_plot) }
    else if (seg == "gene") { add("_gap", s$gap_plot_gene); add("gene", gene_h) }
    else if (seg == "xlabel") add("_gap", s$gap_plot_xlabel)
    else { add(seg, s[[seg]]); add("_gap_panel", s$gap_panel_panel) }
  }
  while (length(nm) && nm[length(nm)] == "_gap_panel") { nm <- nm[-length(nm)]; hh <- hh[-length(hh)] }
  drop <- which(nm == "_gap_panel" & c(nm[-1], "") == "_gap")
  if (length(drop)) { nm <- nm[-drop]; hh <- hh[-drop] }
  keep <- !grepl("^_gap", nm)
  gap_after <- vapply(seq_along(nm), function(i) {
    if (!keep[i]) return(0)
    j <- i + 1L
    if (j <= length(nm) && grepl("^_gap", nm[j])) hh[j] else 0
  }, numeric(1))
  list(names = nm, heights = hh, drawn = hh[keep], gap_after = gap_after[keep],
       total = s$pad_top + sum(hh) + s$pad_bottom, width = s$width)
}
