#!/usr/bin/env Rscript
# Overview schematic: a stylised regional pQTL plot ("identify genetic signals for
# protein levels") for the study-overview figure.
#
# The data are SIMULATED -- this is an illustration, not a result. Background SNPs
# are drawn as noise in LD blocks coloured blue/grey; one signal peaks at a lead SNP,
# its neighbours coloured by a simulated r2 to the lead (yellow -> orange -> red).
#
# A second panel, overview_coloc.svg, shows colocalization: a protein pQTL and a gene
# eQTL track peaking at the same (purple) variant. A third, overview_cistrans.svg,
# shows a cis and a trans pQTL on one protein track. All share one
# theme and palette so they sit together in the schematic.
#
# Output (in output/overview_plot/):
#   overview_locus.svg  (regional pQTL plot)
#   overview_coloc.svg  (pQTL / eQTL colocalization)
#   overview_cistrans.svg (cis- vs trans-pQTL, one coloc-style pQTL track)
#   overview_cistrans_coloc.svg (cis/trans beside an eQTL / pQTL / GWAS coloc stack,
#     one image so both halves share a font size)
# Usage:
#   Rscript code/09_statistical_analysis/overview_plot.R [out_dir]

suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

DEFAULT_OUT <- "output/overview_plot"
args <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- if (length(args) >= 1) args[1] else DEFAULT_OUT
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(7)

# ---- parameters --------------------------------------------------------------
LEAD_POS  <- 45.55      # Mb
LEAD_LOGP <- 14
LEAD_ID   <- "rs123456"
WIN       <- c(44.9, 46.45)
SHADE     <- c(LEAD_POS - 0.08, LEAD_POS + 0.08)
N_BG      <- 1800
N_PEAK    <- 90

# ---- background: noise in LD blocks, blue/grey ------------------------------
bg <- data.frame(pos = sort(runif(N_BG, WIN[1], WIN[2])))
bg$block <- floor((bg$pos - WIN[1]) / 0.06)
block_col <- c("#1f4e9c", "#3a78d0", "#7fb0ea", "#9a9a9a", "#bdbdbd", "#2b62b8")
bg$col <- block_col[(bg$block %% length(block_col)) + 1]
# a few sub-threshold bumps that taper away from their centre, so the background
# is not flat but no second signal competes with the lead
bump_at <- c(45.02, 45.22, 45.43, 45.70, 45.98, 46.20, 46.35)
bump_h  <- c(3, 5, 4, 5.5, 6, 3.5, 6)
near <- sapply(bg$pos, function(x) {
  d <- abs(x - bump_at); i <- which.min(d); bump_h[i] * exp(-d[i] / 0.015)
})
bg$logp <- 0.4 + -log10(runif(N_BG)) * 0.6 + near * runif(N_BG, 0.2, 1)

# ---- signal: SNPs in LD with the lead ---------------------------------------
pk <- data.frame(pos = LEAD_POS + rnorm(N_PEAK, 0, 0.03))
pk$r2 <- pmin(1, exp(-abs(pk$pos - LEAD_POS) / 0.03) * runif(N_PEAK, 0.4, 1.05))
pk$logp <- pmax(0.5, LEAD_LOGP * pk$r2 + rnorm(N_PEAK, 0, 0.6))
pk <- pk[order(pk$r2), ]
lead <- data.frame(pos = LEAD_POS, logp = LEAD_LOGP)

# ---- shared style -----------------------------------------------------------
BASE_SIZE <- 14
NAVY      <- "#1f3a7a"
PURPLE    <- "#7b2cbf"
theme_overview <- theme_classic(base_size = BASE_SIZE) +
  theme(axis.line = element_line(colour = "grey20"),
        axis.text = element_text(colour = "grey20"))
ylab_logp <- expression(-log[10](italic(P)))

# ---- figure 1: regional pQTL plot --------------------------------------------
p <- ggplot() +
  annotate("rect", xmin = SHADE[1], xmax = SHADE[2], ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  geom_point(data = bg, aes(pos, logp, colour = I(col)), size = 1.1, alpha = 0.85) +
  geom_point(data = pk, aes(pos, logp, colour = r2), size = 1.8) +
  scale_colour_gradientn(colours = c("#f7d34a", "#f59a23", "#e8472b", "#d7191c"),
                         limits = c(0, 1), guide = "none") +
  geom_point(data = lead, aes(pos, logp), colour = "#e3261d", size = 4.5) +
  annotate("text", x = LEAD_POS + 0.05, y = LEAD_LOGP + 1, label = LEAD_ID,
           hjust = 0, colour = NAVY, fontface = "bold", size = 5.5) +
  scale_x_continuous(breaks = c(45, 45.5, 46), labels = c("45.0", "45.5", "46.0"),
                     limits = WIN, expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 16), breaks = c(0, 5, 10, 15),
                     expand = c(0, 0)) +
  labs(x = "Genomic position (Mb)", y = ylab_logp) +
  theme_overview

# ---- shared QTL track (figures 2 and 3) --------------------------------------
# Grey background SNPs, one or more peaks, a dashed 5e-8 line; peak SNPs above the
# line take the track colour, those below are grey; each lead is a purple diamond.
# Protein tracks are blue, gene green -- the database portal's entity colours
# (code/proteogenomic_database_app/www/bf_home.css), one flat colour per track.
N_BG2   <- 220
N_PK2   <- 40
TOP2    <- 8
YMAX2   <- 10.5
SIG2    <- 3        # where the 5e-8 line sits on this unlabelled axis
DIAMOND <- 5
BLUE_PAL  <- "#3d6fa5"
GREEN_PAL <- "#3f8f6b"
AMBER_PAL <- "#b5761f"   # disease

# evenly spaced around the lead (plus jitter), so every draw has SNPs close to it
track_peak <- function(lead, top, cols) {
  pk <- data.frame(pos = lead + seq(-0.15, 0.15, length.out = N_PK2) +
                     rnorm(N_PK2, 0, 0.004))
  pk$r2 <- exp(-abs(pk$pos - lead) / 0.06) * runif(N_PK2, 0.6, 1)
  pk$logp <- pmax(0.5, (top - 0.8) * pk$r2 + rnorm(N_PK2, 0, 0.3))
  pk$col <- colorRampPalette(cols)(100)[ceiling(pmax(0.01, pk$r2) * 100)]
  pk[order(pk$r2), ]
}

# `under`: layers drawn beneath the points (the coloc line, the cis band); `win` is
# the x range, and the background SNP count scales with its width so a narrower
# track keeps the same density
qtl_track <- function(seed, leads, tops, cols, under = NULL, win = WIN) {
  set.seed(seed)
  n_bg <- round(N_BG2 * diff(win) / diff(WIN))
  bg <- data.frame(pos = runif(n_bg, win[1], win[2]))
  bg$logp <- 0.15 + rexp(n_bg, 1 / 0.45)
  pk <- do.call(rbind, lapply(seq_along(leads), function(i)
    track_peak(leads[i], tops[i], cols)))
  sig <- pk$logp >= SIG2
  ggplot() + under +
    geom_point(data = rbind(bg, pk[!sig, c("pos", "logp")]), aes(pos, logp),
               colour = "grey65", size = 1.6, alpha = 0.85) +
    geom_hline(yintercept = SIG2, linetype = "dashed", colour = "grey55") +
    geom_point(data = pk[sig, ], aes(pos, logp, colour = I(col)), size = 2.2) +
    annotate("point", x = leads, y = tops, shape = 18, size = DIAMOND,
             colour = PURPLE) +
    scale_x_continuous(limits = win, expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, YMAX2), expand = c(0, 0)) +
    labs(x = "Genomic position", y = ylab_logp) +
    theme_overview +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_text(margin = margin(t = 8)))
}

# ---- figure 2: pQTL / eQTL colocalization ------------------------------------
# One simulated causal variant at LEAD_POS drives both traits; each track gets its
# own noise so the two peaks look alike but not identical.
# A narrower window centred on the lead keeps the peak in the middle.
COLOC_WIN  <- LEAD_POS + c(-0.3, 0.3)
coloc_line <- geom_vline(xintercept = LEAD_POS, linetype = "dashed", colour = "grey55")
p_pqtl <- qtl_track(11, LEAD_POS, TOP2, BLUE_PAL, coloc_line, win = COLOC_WIN)
p_eqtl <- qtl_track(12, LEAD_POS, TOP2, GREEN_PAL, coloc_line, win = COLOC_WIN)
p2 <- p_pqtl | p_eqtl

# ---- figure 3: cis vs trans --------------------------------------------------
# One protein track with a cis peak inside a grey band marking the cis region and
# a trans peak outside it.
CIS_POS   <- 45.30
TRANS_POS <- 46.05
CIS_WIN   <- 0.22                  # half-width of the cis region, grey band
CT_WIN    <- c(44.98, 46.32)       # x range, trimmed of empty background at both ends
# The band's label matches the pipeline's cis definition (code/02_protein_qc/
# 03_annot_ensg.R: gene start/end -/+ 500 kb). It sits beside the band with a leader
# line, keeping the space above the peaks for the cis-/trans-pQTL labels.
NOTE_Y <- 5.8
cis_band <- list(
  annotate("rect", xmin = CIS_POS - CIS_WIN, xmax = CIS_POS + CIS_WIN,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6),
  annotate("segment", x = CIS_POS + CIS_WIN - 0.03, xend = CIS_POS + CIS_WIN + 0.04,
           y = NOTE_Y, yend = NOTE_Y, colour = "grey35", linewidth = 0.4),
  annotate("text", x = CIS_POS + CIS_WIN + 0.05, y = NOTE_Y, hjust = 0,
           label = "encoding gene\n±500 kb", colour = "grey35", size = 4.6,
           lineheight = 0.9))
p3 <- qtl_track(13, c(CIS_POS, TRANS_POS), c(TOP2, 7), BLUE_PAL, cis_band,
                win = CT_WIN) +
  annotate("text", x = c(CIS_POS, TRANS_POS), y = c(TOP2, 7) + 1.2,
           label = c("bolditalic(\"cis\")*bold(\"-pQTL\")", "bolditalic(\"trans\")*bold(\"-pQTL\")"),
           parse = TRUE, colour = BLUE_PAL, size = 4.6)

# ---- figure 4: cis/trans beside a coloc stack ----------------------------------
# eQTL -> pQTL -> GWAS share one lead; each track names itself top-left, only the
# middle one carries the y title and only the bottom one the x title.
stack_track <- function(seed, col, lab, pos = c("top", "mid", "bottom")) {
  pos <- match.arg(pos)
  qtl_track(seed, LEAD_POS, TOP2, col, coloc_line, win = COLOC_WIN) +
    annotate("text", x = COLOC_WIN[1], y = Inf, label = lab, colour = col,
             fontface = "bold", size = 4.6, hjust = 0, vjust = 1.3) +
    theme(axis.title.y = if (pos == "mid") element_text() else element_blank(),
          axis.title.x = if (pos == "bottom") element_text(margin = margin(t = 8))
                         else element_blank(),
          plot.margin = margin(2, 5.5, 2, 5.5))
}
p_stack <- stack_track(12, GREEN_PAL, "eQTL", "top") /
  stack_track(11, BLUE_PAL, "pQTL", "mid") /
  stack_track(26, AMBER_PAL, "GWAS", "bottom")
p4 <- (p3 | plot_spacer() | p_stack) + plot_layout(widths = c(4.2, 0.6, 3.2))

# ---- write -------------------------------------------------------------------
outs <- list(overview_locus.svg = list(p, 5.5, 3.2),
             overview_coloc.svg = list(p2, 5, 2.5),
             overview_cistrans.svg = list(p3, 4.2, 3.8),
             overview_cistrans_coloc.svg = list(p4, 8.2, 3.8))
for (nm in names(outs)) {
  f <- file.path(OUT_DIR, nm)
  ggsave(f, outs[[nm]][[1]], width = outs[[nm]][[2]], height = outs[[nm]][[3]])
  if (!file.exists(f)) stop("not written: ", f)
}
message("overview_plot.R: wrote ", paste(names(outs), collapse = ", "), " to ", OUT_DIR)
