#!/usr/bin/env Rscript
# Pathway enrichment of each trans-pQTL hotspot's target proteins, per fluid.
#
# For every hotspot that clears HOTSPOT_MIN in a fluid (01b), its target genes in that
# fluid are tested against every gene set by a one-sided Fisher test. The background is
# the genes measured in that fluid (00_kept_assays.csv), not the genome: the panels are
# already rich in secreted, glycosylated proteins, and a genome background would call
# that panel composition a hotspot property.
#
# Gene sets: MSigDB 2024.1 Hallmark, Reactome and GO BP/CC/MF (symbols), downloaded to
# GENESET_DIR (URL below), and UniProt keywords (info/bf_protein/uniprot/features.csv).
# Each set is restricted to the background; sets with fewer than MIN_SET or more than
# MAX_SET measured genes are not tested. BH-FDR runs within each hotspot x fluid over
# all collections together.
#
# Output: output/01b_hotspot_enrichment.csv, one row per hotspot x fluid x tested set.
suppressPackageStartupMessages(library(data.table))
OUT_DIR <- "output"
source("code/09_statistical_analysis/lib/common.R")

# https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/<name>.v2024.1.Hs.symbols.gmt
GENESET_DIR <- path.expand("~/db/shared_data/ref_geno/geneset")
GMT <- c(Hallmark = "h.all", Reactome = "c2.cp.reactome", `GO BP` = "c5.go.bp",
         `GO CC` = "c5.go.cc", `GO MF` = "c5.go.mf")
MIN_SET <- 5L
MAX_SET <- 500L

read_gmt <- function(f) {
  l <- strsplit(readLines(f), "\t", fixed = TRUE)
  rbindlist(lapply(l, function(x) data.table(set = x[1L], gene = x[-(1:2)])))
}
gmt_files <- file.path(GENESET_DIR, paste0(GMT, ".v2024.1.Hs.symbols.gmt"))
if (!all(file.exists(gmt_files))) stop("missing gene set files in ", GENESET_DIR)
gs <- rbindlist(Map(function(f, coll) read_gmt(f)[, collection := coll], gmt_files, names(GMT)))

# UniProt keywords, mapped to gene through the assay annotation.
feat <- fread("info/bf_protein/uniprot/features.csv",
              select = c("uniprot_id", "keywords"))
annot <- unique(fread("info/bf_protein/assay_annot.csv",
                      select = c("gene", "uniprot_id"))[uniprot_id != ""])
kw <- merge(annot, feat, by = "uniprot_id")[keywords != ""]
kw <- kw[, .(set = trimws(unlist(strsplit(keywords, ";", fixed = TRUE)))), by = gene]
# Bookkeeping keywords say how a protein was studied, not what it does.
kw <- kw[!set %in% c("Reference proteome", "Proteomics identification", "3D-structure",
                     "Alternative splicing", "Direct protein sequencing", "Polymorphism")]
gs <- rbind(gs, unique(kw)[, .(set = paste0("KW_", toupper(gsub("[^A-Za-z0-9]+", "_", set))),
                               gene, collection = "UniProt keyword")])
gs <- unique(gs)

kept <- fread(file.path(OUT_DIR, "00_kept_assays.csv"), select = c("fluid", "gene"))
hs <- fread(file.path(OUT_DIR, "01b_hotspot_summary.csv"))
stopifnot(c("target_genes_csf", "target_genes_plasma") %in% names(hs))

hs_long <- rbind(
  hs[is_hotspot_csf == TRUE, .(regulator, lead_snp, fluid = "CSF", targets = target_genes_csf)],
  hs[is_hotspot_plasma == TRUE, .(regulator, lead_snp, fluid = "Plasma", targets = target_genes_plasma)])

res <- rbindlist(lapply(seq_len(nrow(hs_long)), function(i) {
  h <- hs_long[i]
  bg <- unique(kept[fluid == h$fluid, gene])
  tg <- intersect(strsplit(h$targets, ";", fixed = TRUE)[[1L]], bg)
  stopifnot(length(tg) >= HOTSPOT_MIN)
  g <- gs[gene %chin% bg]
  g <- g[, .(n_set = .N, n_hit = sum(gene %chin% tg),
             hits = paste(sort(gene[gene %chin% tg]), collapse = ";")), by = .(collection, set)]
  g <- g[n_set >= MIN_SET & n_set <= MAX_SET]
  N <- length(bg); k <- length(tg)
  # one-sided Fisher = hypergeometric upper tail
  g[, p := phyper(n_hit - 1L, n_set, N - n_set, k, lower.tail = FALSE)]
  g[, `:=`(fold = (n_hit / k) / (n_set / N), fdr = p.adjust(p, "BH"))]
  g[, `:=`(regulator = h$regulator, lead_snp = h$lead_snp, fluid = h$fluid,
           n_targets = k, n_background = N)]
}))
setcolorder(res, c("regulator", "lead_snp", "fluid", "n_targets", "n_background",
                   "collection", "set", "n_set", "n_hit", "fold", "p", "fdr", "hits"))
setorder(res, regulator, fluid, p)

out <- file.path(OUT_DIR, "01b_hotspot_enrichment.csv")
fwrite(res, out)
stopifnot(file.exists(out))
message("Top set per hotspot x fluid (FDR < 0.05 count in n_sig):")
print(res[, .(n_targets = n_targets[1], n_sig = sum(fdr < 0.05), top = set[1],
              hit = paste0(n_hit[1], "/", n_set[1]), fold = round(fold[1], 1),
              fdr = signif(fdr[1], 2)), by = .(regulator, fluid)])
message("01b_hotspot_enrichment.R done: ", nrow(res), " tests, ", out)
