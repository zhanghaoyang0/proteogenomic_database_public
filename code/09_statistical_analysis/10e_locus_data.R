#!/usr/bin/env Rscript
# Data for the three locus panels of the paired-fluid figure: a brain eQTL,
# the CSF pQTL and the plasma pQTL of the same gene, as -log10 P coloured by LD with
# the lead of the fluid that carries the signal (CSF for GRN and ROR1, plasma for THY1).
#
#   GRN  - MetaBrain cortex       shared: eQTL, CSF and plasma peak at rs5848
#   ROR1 - MetaBrain cerebellum   CSF only: eQTL and CSF peak at rs3861946 (EAF 0.56),
#                                 plasma P > 1e-4 across the window
#   THY1 - MetaBrain cortex       plasma only: eQTL and plasma peak at rs10736497,
#                                 CSF P > 1e-3 across the window
#
# ROR1 was the only CSF-only brain colocalization (head-to-head proteins) with a CSF
# peak at P < 1e-8, a brain eQTL at P < 1e-6 and no plasma variant at P < 1e-4 within
# 50 kb (scan of 174 gene-context pairs, 2026-10-27). The window is centred on the lead
# fluid's peak nearest the gene (the eQTL file spans the gene's cis window, so its midpoint is
# the gene); the coloc table's pos column is not the lead and must not be used.
#
# Analysis only; reads the clean GWAS from the secure tree so 10d never touches it.
# Outputs (output/):
#   10e_locus_tracks.csv.gz  gene label lead_fluid track chr pos p r2 lead
# Usage (repo root):
#   Rscript code/09_statistical_analysis/10e_locus_data.R

suppressPackageStartupMessages(library(data.table))
source("code/lib/locus_plot.R")         # ld_panel_init, ld_r2

OUT_DIR <- "output"
source("code/09_statistical_analysis/lib/common.R")   # GWAS_BASE, LD_REF
GWAS <- GWAS_BASE
LD_BFILE <- LD_REF
HALF <- 150e3

CASES <- data.table(
  gene = c("GRN", "ROR1", "THY1"),
  label = c("shared", "CSF only", "plasma only"),
  ensg = c("ENSG00000030582", "ENSG00000185483", "ENSG00000154096"),
  resource = "metabrain",
  tissue = c("cortex", "cerebellum", "cortex"),
  eqtl_track = c("Cortex eQTL", "Cerebellum eQTL", "Cortex eQTL"),
  lead_fluid = c("CSF", "CSF", "Plasma")
)

kept <- fread(file.path(OUT_DIR, "00_kept_assays.csv"))[platform != "SomaLogic" & gene %in% CASES$gene]
stopifnot(kept[, .N, by = .(gene, fluid)][, all(N == 1) && .N == 2 * nrow(CASES)])
assay <- function(g, f) file.path(GWAS, "bf2/clean/protein",
                                  paste0("bf2__", kept[gene == g & fluid == f, protein_assay], ".tsv.gz"))

read_win <- function(f, chr, from, to) {
  if (!file.exists(f)) stop("missing ", f)
  d <- fread(f, select = c("CHR", "POS", "P"))[CHR == chr & POS >= from & POS <= to]
  if (!nrow(d)) stop("no variants in window: ", f)
  setorder(d, POS, P)
  unique(d, by = "POS")   # one row per position, the most significant
}

ld_panel_init(LD_BFILE)
tracks <- list()
for (i in seq_len(nrow(CASES))) {
  cs <- CASES[i]
  ef <- file.path(GWAS, cs$resource, "clean", cs$tissue, paste0(cs$ensg, ".tsv.gz"))
  eq <- fread(ef, select = c("CHR", "POS"))
  chr <- eq$CHR[1]; gmid <- mean(range(eq$POS))
  near <- read_win(assay(cs$gene, cs$lead_fluid), chr, gmid - 250e3, gmid + 250e3)
  lead <- near[which.min(P), POS]
  from <- lead - HALF; to <- lead + HALF
  d <- rbind(read_win(ef, chr, from, to)[, track := cs$eqtl_track],
             read_win(assay(cs$gene, "CSF"), chr, from, to)[, track := "CSF pQTL"],
             read_win(assay(cs$gene, "Plasma"), chr, from, to)[, track := "Plasma pQTL"])
  ref <- .LD$bim[CHR == as.character(chr) & POS == lead, SNP[1]]
  ld <- ld_r2(chr, from, to, if (length(ref)) ref else "", ref_pos = lead)
  if (is.null(ld)) stop("no LD for ", cs$gene)
  d[, r2 := ld$r2[match(POS, ld$POS)]]
  d[POS == lead, r2 := 1]
  d[, `:=`(gene = cs$gene, label = cs$label, lead_fluid = cs$lead_fluid,
           lead = POS == lead)]
  top <- d[, .(pos = POS[which.min(P)], p = min(P)), by = track]
  message(cs$gene, ": ", chr, ":", from, "-", to, "; lead ", lead, " (", ref, "); tops ",
          top[, paste0(track, " ", pos, " ", signif(p, 2), collapse = " | ")])
  tracks[[i]] <- d
}
tracks <- rbindlist(tracks)
setnames(tracks, c("CHR", "POS", "P"), c("chr", "pos", "p"))
setcolorder(tracks, c("gene", "label", "lead_fluid", "track", "chr", "pos", "p", "r2", "lead"))

# The panels claim these patterns; stop if the inputs no longer show them.
pk <- tracks[, .(top = pos[which.min(p)], p = min(p), lead = pos[lead][1]), by = .(gene, track)]
stopifnot(pk[gene == "GRN", all(top == lead & p < 5e-8)],
          pk[gene == "ROR1" & track != "Plasma pQTL", all(top == lead & p < 5e-8)],
          pk[gene == "ROR1" & track == "Plasma pQTL", p > 1e-4],
          pk[gene == "THY1" & track != "CSF pQTL", all(top == lead & p < 5e-8)],
          pk[gene == "THY1" & track == "CSF pQTL", p > 1e-3])
tracks[, `:=`(p = signif(p, 4), r2 = signif(r2, 3))]
fwrite(tracks, file.path(OUT_DIR, "10e_locus_tracks.csv.gz"))
if (!file.exists(file.path(OUT_DIR, "10e_locus_tracks.csv.gz"))) stop("not written")
message("10e_locus_data done: ", nrow(tracks), " rows")
