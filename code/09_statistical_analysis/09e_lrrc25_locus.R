#!/usr/bin/env Rscript
# Extract complete LRRC25 lead +/-150 kb summary-statistic windows and public-reference LD.
# Run from the repository root on a compute node; no new fine-mapping or coloc is performed.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
root <- normalizePath(getwd())
project <- root
dest <- file.path(project, "output")
source(file.path(root, "code/lib/locus_plot.R"))
lead <- fread(file.path(project, "output/09d_novel_evidence.csv"))[gene == "LRRC25"]
stopifnot(nrow(lead) == 1L, lead$snp == "rs1985157")
chr <- lead$chr; pos <- lead$pos; from <- pos - 150000L; to <- pos + 150000L
gw <- path.expand("~/db/shared_data/gwas")
files <- c("AD GWAS" = file.path(gw, "other/clean/ad_uffelmann_2026.tsv.gz"),
           "CSF pQTL" = file.path(gw, "bf2/clean/protein/bf2__LRRC25__Oncology_o2csf.tsv.gz"),
           "Plasma pQTL" = file.path(gw, "bf2/clean/protein/bf2__LRRC25__ExploreHT_o7plasma.tsv.gz"))
tracks <- rbindlist(lapply(names(files), function(track_name) {
  f <- files[[track_name]]
  stopifnot(file.exists(f))
  # Stream the complete file; retain every variant in this declared regional window.
  filt <- sprintf('NR==1 || ($1==%d && $2>=%d && $2<=%d)', chr, from, to)
  cmd <- paste("gzip -dc", shQuote(f), "| awk", shQuote(filt))
  x <- fread(cmd = cmd, select = c("CHR", "POS", "SNP", "EA", "NEA", "BETA", "SE", "P"))
  stopifnot(nrow(x) > 10L, all(is.finite(x$P)), all(x$P > 0 & x$P <= 1),
            all(x$CHR == chr), any(x$SNP == lead$snp))
  x[, track := track_name]
  message(track_name, ": ", nrow(x), " variants retained in the +/-150 kb window")
  x
}))
ld_panel_init("~/db/shared_data/ref_geno/1kg_v3/EUR")
stopifnot(.LD$bim[CHR == as.character(chr) & SNP == lead$snp & POS == pos, .N] == 1L)
ld <- ld_r2(chr, from, to, lead$snp, ref_pos = pos)
stopifnot(!is.null(ld), ld[SNP == lead$snp, abs(r2-1) < 1e-8])
tracks[, r2 := ld$r2[match(SNP, ld$SNP)]]
tracks[, lead_variant := SNP == lead$snp]
fwrite(tracks, file.path(dest, "09e_lrrc25_tracks.csv.gz"))
co <- fread(path.expand("~/db/shared_data/coloc/protein_x_disease/result/coloc_protein_x_disease.tsv"))
ev <- co[outcome_id == "ad_uffelmann_2026" & grepl("LRRC25__", exposure_id, fixed = TRUE)]
stopifnot(nrow(ev) > 0)
fwrite(ev, file.path(dest, "09e_lrrc25_coloc.csv"))
manifest <- data.table(track = names(files), input_file = basename(files),
                       chr = chr, start = from, end = to,
                       lead_snp = lead$snp,
                       input_md5 = unname(tools::md5sum(files)))
fwrite(manifest, file.path(dest, "09e_lrrc25_provenance.csv"))
stopifnot(file.info(file.path(dest, "09e_lrrc25_tracks.csv.gz"))$size > 0)
message("LRRC25 regional data ready")
