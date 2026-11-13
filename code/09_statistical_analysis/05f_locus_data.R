#!/usr/bin/env Rscript
# Data for the two locus panels of the redesigned protein-disease figure: one
# fluid-specific and one shared protein-disease link, drawn as stacked regional tracks
# (disease GWAS, CSF pQTL, plasma pQTL), each -log10 P coloured by LD with the disease
# lead and a SuSiE PIP strip below (layout after the GWAS / eQTL locus panels the
# author pointed to, 2026-10-27).
#
#   CD58 - multiple sclerosis   CSF only: no plasma cis instrument. The plasma GWAS
#                               still peaks at the same variants (P ~ 3e-5, r2 > 0.8
#                               with the MS lead): an attenuated, sub-threshold
#                               effect, not an absent one (found 2026-10-27)
#   GRN  - Alzheimer disease    joint support in both fluids, peaks aligned
#
# The plasma track is drawn even when it has no fine-mapped signal, so the reader sees
# what the plasma GWAS does at the locus. PIP is NA where a track has no SuSiE fit.
#
# Analysis only; reads the clean GWAS, the SuSiE fits and the coloc table from the
# secure tree, so 05g (the plot) never touches it. Windows, LD and genes come from
# code/lib/locus_plot.R, the helpers the website's locus plots use.
#
# Outputs (output/):
#   05f_locus_tracks.csv.gz  gene disease track chr pos snp p pip in_cs r2
#   05f_locus_genes.csv      gene chr start end strand symbol (protein-coding genes
#                            in each window)
#   05f_locus_meta.csv       one row per case x fluid with a coloc result: assay,
#                            PP4, protein and disease leads
# Usage (repo root):
#   Rscript code/09_statistical_analysis/05f_locus_data.R

suppressPackageStartupMessages(library(data.table))

source("code/07_colocalization/functions_coloc.R")  # protein_susie_path, disease_susie_path,
                                                 # protein_src, coloc_reported
source("code/lib/locus_plot.R")         # ld_panel_init, ld_r2, gene_bed_init

OUT_DIR <- "output"
LD_BFILE <- LD_REF                                            # code/lib/common.R,
GENE_BED <- file.path(COLOC_ROOT, "info/ensembl87_grch37_genes.bed.gz")   # via functions_coloc.R
COLOC_FILE <- file.path(COLOC_ROOT, "protein_x_disease/result/coloc_protein_x_disease.tsv")

CASES <- data.table(
  gene = c("CD58", "GRN"),
  outcome_id = c("G6_MS_meta_fin_ukbb_mvp", "finngen_R12_G6_ALZHEIMER_meta_fin_ukbb"),
  disease = c("Multiple sclerosis", "Alzheimer disease")
)

clean_path <- function(trait) path.expand(file.path(
  GWAS_BASE,
  if (startsWith(trait, "bf2__")) "bf2/clean/protein" else "finngen/clean",
  paste0(trait, ".tsv.gz")))

# ---- assays and coloc rows ------------------------------------------------------
kept <- fread(file.path(OUT_DIR, "00_kept_assays.csv"))[
  platform != "SomaLogic" & gene %in% CASES$gene]
kept[, exposure_id := paste0("bf2__", protein_assay)]
stopifnot(kept[, .N, by = .(gene, fluid)][, all(N == 1) && .N == 2 * nrow(CASES)])

cod <- fread(path.expand(COLOC_FILE), select = c("exposure_id", "outcome_id", "status",
             "PP4", "lead_snp", "o_lead", "chr", "pos", "primary"))
cod <- merge(cod, kept[, .(exposure_id, gene, fluid)], by = "exposure_id")
cod <- merge(cod, CASES, by = c("gene", "outcome_id"))
cod <- cod[coloc_reported(status, primary, PP4)]
setorder(cod, gene, fluid, -PP4)
best <- cod[, .SD[1], by = .(gene, fluid)]
# CD58 colocalizes in CSF only, GRN in both: anything else means the inputs moved
stopifnot(setequal(best[, paste(gene, fluid)], c("CD58 CSF", "GRN CSF", "GRN Plasma")),
          best[, all(PP4 > PP4_STRONG)])   # code/lib/common.R

# The coloc result behind each panel; 05g reads disease and disease_lead. The MR
# columns it once carried had no reader (code review 2026-10-30, R17).
meta <- merge(best[, .(gene, fluid, outcome_id, exposure_id, pp4 = PP4,
                       protein_lead = lead_snp, disease_lead = o_lead, chr)],
              CASES[, .(gene, disease)], by = "gene")

# ---- tracks ----------------------------------------------------------------------
ld_panel_init(LD_BFILE)
gene_bed_init(GENE_BED)

# PIP and credible-set membership from a SuSiE fit, by position (fits carry rsids,
# the protein sumstats chr:pos:alleles ids, so position is the shared key).
fit_pip <- function(path, chr) {
  if (!file.exists(path)) stop("Missing SuSiE fit: ", path)
  fit <- readRDS(path)
  cs_snp <- fit$snp[unique(unlist(fit$sets$cs, use.names = FALSE))]
  pos <- .LD$bim[CHR == as.character(chr), .(snp = SNP, pos = POS)]
  d <- merge(data.table(snp = fit$snp, pip = fit$pip), pos, by = "snp")
  stopifnot(nrow(d) > 0.9 * length(fit$snp))
  d[, .(pos, pip, in_cs = snp %in% cs_snp)]
}

track_rows <- function(trait, chr, from, to, label, fit = NULL) {
  ss <- fread(clean_path(trait), select = c("CHR", "POS", "SNP", "P"))[
    CHR == chr & POS >= from & POS <= to]
  if (!nrow(ss)) stop("no variants in window for ", trait)
  # one row per position, the most significant, so the PIP and r2 joins cannot double
  setorder(ss, POS, P)
  d <- unique(ss, by = "POS")[, .(pos = POS, snp = SNP, p = P)]
  if (is.null(fit)) {
    d[, `:=`(pip = NA_real_, in_cs = FALSE)]
  } else {
    d <- merge(d, fit[!duplicated(pos)], by = "pos", all.x = TRUE)
    d[is.na(in_cs), in_cs := FALSE]
  }
  d[, track := label]
  d
}

tracks <- list(); genes <- list()
for (i in seq_len(nrow(CASES))) {
  cs <- CASES[i]
  csf <- best[gene == cs$gene & fluid == "CSF"]
  chr <- as.character(csf$chr)
  # the window is the CSF SuSiE fit's, so every track covers what fine-mapping saw
  pfit <- fit_pip(protein_susie_path(protein_src(csf$exposure_id), csf$exposure_id,
                                     csf$lead_snp), chr)
  from <- min(pfit$pos); to <- max(pfit$pos)
  dfit <- fit_pip(disease_susie_path("finngen", cs$outcome_id, csf$o_lead), chr)

  pl <- best[gene == cs$gene & fluid == "Plasma"]
  plfit <- if (nrow(pl)) fit_pip(protein_susie_path(protein_src(pl$exposure_id),
                                                    pl$exposure_id, pl$lead_snp), chr)
  pl_assay <- kept[gene == cs$gene & fluid == "Plasma", exposure_id]

  d <- rbind(
    track_rows(cs$outcome_id, chr, from, to, "Disease", dfit),
    track_rows(csf$exposure_id, chr, from, to, "CSF", pfit),
    track_rows(pl_assay, chr, from, to, "Plasma", plfit))

  # r2 to the disease lead, the variant every track is read against
  lead_pos <- fread(clean_path(cs$outcome_id), select = c("SNP", "POS"))[
    SNP == csf$o_lead, POS[1]]
  ld <- ld_r2(chr, from, to, csf$o_lead, ref_pos = lead_pos)
  if (is.null(ld)) stop("no LD for ", cs$gene)
  d[, r2 := ld$r2[match(pos, ld$POS)]]
  d[pos == lead_pos, r2 := 1]
  d[, `:=`(gene = cs$gene, disease = cs$disease, chr = chr)]
  tracks[[i]] <- d

  g <- .GENES$g[CHR == chr & END >= from & START <= to,
                .(chr = CHR, start = START, end = END, strand = STRAND, symbol = SYMBOL)]
  genes[[i]] <- g[, gene := cs$gene]
  message(cs$gene, ": ", d[, .N, by = track][, paste(track, N, collapse = ", ")],
          "; window ", chr, ":", from, "-", to, "; ", nrow(g), " genes")
}
tracks <- rbindlist(tracks)
setcolorder(tracks, c("gene", "disease", "track", "chr", "pos", "snp", "p", "pip", "in_cs", "r2"))
tracks[, `:=`(p = signif(p, 4), pip = signif(pip, 4), r2 = signif(r2, 3))]
# the panels promise CS lines in the disease and CSF tracks of both cases
stopifnot(tracks[track %in% c("Disease", "CSF"), any(in_cs), by = .(gene, track)][, all(V1)])

fwrite(tracks, file.path(OUT_DIR, "05f_locus_tracks.csv.gz"))
fwrite(rbindlist(genes), file.path(OUT_DIR, "05f_locus_genes.csv"))
fwrite(meta, file.path(OUT_DIR, "05f_locus_meta.csv"))
for (f in c("05f_locus_tracks.csv.gz", "05f_locus_genes.csv", "05f_locus_meta.csv"))
  if (!file.exists(file.path(OUT_DIR, f))) stop("not written: ", f)
message("05f_locus_data done")
