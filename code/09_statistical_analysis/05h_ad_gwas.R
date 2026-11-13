# Colocalization of the analysed CSF and plasma proteins with two Alzheimer's disease
# GWAS outside the FinnGen disease set: Bellenguez et al. 2022 and the larger PGC-ALZ3
# meta-analysis (Uffelmann et al. 2026, EUR case-control; ref/nature_ad_gwas_2026.md).
# Same rule as 05a: a reported component (coloc_reported) with PP4 > COLOC_PP4, on the
# Olink/NuLISA assays of the deduped analysed set (load_assay_annot). Read by Results 2
# (LRRC25) and Results 5; Fig. 6 and every 05a count stay FinnGen-only.
#
# Server (reads the secure coloc table). Writes:
#   output/05h_ad_coloc.csv   one row per outcome x gene with a strong colocalization:
#                             fluids, max PP4, measured in both fluids, and the PGC-ALZ3
#                             novel locus (Table 1 of the paper) the disease lead falls in
#
# Usage (repo root): Rscript code/09_statistical_analysis/05h_ad_gwas.R

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/common.R")             # coloc_reported, COLOC_PP4
source("code/09_statistical_analysis/lib/eqtl_coloc_scope.R")   # load_assay_annot

OUT     <- "output/05h_ad_coloc.csv"
AD      <- c("ad_bellenguez_2022", "ad_uffelmann_2026")
COLOC   <- file.path(COLOC_ROOT, "protein_x_disease/result/coloc_protein_x_disease.tsv")
AD_FM   <- file.path(GWAS_BASE, "other/finemapping/ad_uffelmann_2026.finemapped.tsv.gz")
PAPER   <- "ref/nature_ad_gwas_2026.md"
LOCUS_PAD <- 1e5   # a disease lead within 100 kb of a novel locus's bounds falls in it

annot <- load_assay_annot()
both <- annot[, uniqueN(fluid), by = gene][V1 == 2, gene]

co <- fread(COLOC)[outcome_id %in% AD]
if (!setequal(unique(co$outcome_id), AD)) stop("coloc table lacks ", setdiff(AD, co$outcome_id))
n0 <- nrow(co)
co <- merge(co, annot[, .(exposure_id, fluid, gene)], by = "exposure_id")   # inner: drops SomaLogic
message(sprintf("AD coloc rows: %d, %d on analysed assays", n0, nrow(co)))
s <- co[(coloc_reported(status, primary, PP4) & PP4 > COLOC_PP4) %in% TRUE]

# Disease-lead positions (hg19) from the PGC-ALZ3 fine-map; APOE-region leads are absent
# from it (no position, so never matched to a novel locus, none of which is on chr19 q13.32).
fm <- unique(fread(AD_FM, select = c("snp", "chr", "pos")))[, .(o_lead = snp, o_chr = chr, o_pos = pos)]
if (anyDuplicated(fm$o_lead)) stop("a fine-mapped variant has two positions")
n1 <- nrow(s)
s <- merge(s, fm, by = "o_lead", all.x = TRUE)
if (nrow(s) != n1) stop("lead-position merge changed the row count")
# A PGC-ALZ3 lead without a position can never match a novel locus; allowed only in the
# APOE region, judged by the coloc window (code review 2026-10-30, R14).
no_pos <- s[outcome_id == "ad_uffelmann_2026" & is.na(o_pos)]
if (nrow(no_pos) && !all(in_apoe_region(no_pos$chr, no_pos$pos)))
  stop("PGC-ALZ3 leads without a fine-map position outside APOE: ",
       paste(unique(no_pos[!in_apoe_region(chr, pos), paste(gene, o_lead)]), collapse = ", "))

# Novel loci: Table 1 of the paper (locus, predicted effector, chromosome, start, end).
md <- readLines(PAPER)
rows <- grep("^\\|[0-9]+[a-z]?\\|", md, value = TRUE)
tb <- rbindlist(lapply(strsplit(rows, "|", fixed = TRUE), function(x)
  data.table(locus = x[2], effector = gsub("_", "", x[3]), chr = x[4], start = x[5], end = x[6])))
tb <- tb[grepl("^[0-9]+$", chr)]
tb[, `:=`(chr = as.integer(chr), start = as.numeric(gsub("[^0-9]", "", start)),
          end = as.numeric(gsub("[^0-9]", "", end)))]
if (nrow(tb) != 48) stop("expected the 48 novel loci of Table 1, parsed ", nrow(tb))
s[, novel_locus := NA_character_]
for (i in seq_len(nrow(tb)))
  s[outcome_id == "ad_uffelmann_2026" & o_chr == tb$chr[i] &
    o_pos >= tb$start[i] - LOCUS_PAD & o_pos <= tb$end[i] + LOCUS_PAD,
    novel_locus := paste0(tb$locus[i], " (", tb$effector[i], ")")]

res <- s[, .(fluids = paste(sort(unique(fluid)), collapse = "+"), max_pp4 = max(PP4),
             novel_locus = paste(sort(unique(na.omit(novel_locus))), collapse = "; ")),
         by = .(outcome_id, gene)]
res[, in_both := gene %in% both]
setorder(res, outcome_id, fluids, gene)
print(res)
fwrite(res, OUT)
if (!file.exists(OUT)) stop("not written: ", OUT)
message("05h_ad_gwas done -> ", OUT)
