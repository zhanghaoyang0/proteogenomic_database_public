#!/usr/bin/env Rscript
# Functional annotation of the fine-mapped pQTL signals: what the credible-set variants
# are (VEP consequence) and where they sit (Roadmap chromatin states), compared between
# signals shared by CSF and plasma and signals seen in one fluid only.
#
# Units are the reported distinct signals (lib/cs_signals.R). Every member variant
# counts by its PIP, so a signal contributes one unit of probability however large its
# credible set. Fluid specificity is 03's credible-set call (03_fluid_specificity_
# credible_sets.csv), defined only for proteins assayed in both fluids.
#
# Consequence: VEP (code/03_pqtl_mapping, one picked consequence per variant), joined
# by rsid; each variant takes the most severe class in CSQ_CLASS order.
#
# Chromatin: Roadmap Epigenomics core 15-state ChromHMM, 127 epigenomes (hg19, as the
# fine-mapping); URLs in ROADMAP_DIR/url.txt. Active regulatory = TssA, TssAFlnk, EnhG,
# Enh. For one signal and one epigenome, observed = PIP mass of its members in active
# states; expected = the share of the +/-500 kb around its top-PIP variant covered by
# active states (the PIP mass a signal would put there by chance). Enrichment of a set
# of signals = sum(observed) / sum(expected); 95% CI and the CSF-specific versus
# plasma-specific contrast by bootstrapping signals (N_BOOT). The MHC is left out.
#
# Outputs (output/): 03b_consequence.csv, 03b_pav_test.csv,
# 03b_chromatin_enrichment.csv, 03b_chromatin_contrast.csv
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})
OUT_DIR <- "output"
source("code/09_statistical_analysis/lib/common.R")
source("code/09_statistical_analysis/lib/cs_signals.R")

FM_DIR <- file.path(GWAS_BASE, "bf2/finemapping/protein")
VEP_TSV <- file.path(GWAS_BASE, "bf2/annotation/protein/variant_annotation.tsv.gz")
ROADMAP_DIR <- path.expand("~/db/shared_data/ref_geno/roadmap_chromhmm15")
ACTIVE <- c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh")
HALF_WIN <- 5e5
MHC <- fread("code/lib/mhc_hg19.range", header = FALSE,   # hg19, the repo-wide range
             col.names = c("chr", "start", "end", "label"))
N_BOOT <- 1000L
set.seed(1)

# VEP consequence classes: CSQ_CLASS and csq_class() in lib/common.R.

# ---- signals and their members ------------------------------------------------
kept <- fread(file.path(OUT_DIR, "00_kept_assays.csv"))
cs <- rbindlist(Map(read_finemap_cs, kept$cohort, kept$protein_assay,
                    MoreArgs = list(fm_dir = FM_DIR,
                                    select = c("cs_id", "snp", "chr", "pos", "pip", "is_cis"))))
check_finemap_cs(cs, kept, FM_DIR)
sig <- distinct_signals(unique(cs[, .(cohort, protein_assay, cs_id)]))
cs <- merge(cs, sig[, .(cohort, protein_assay, cs_id)], by = c("cohort", "protein_assay", "cs_id"))
cs <- merge(cs, kept[, .(cohort, protein_assay, fluid, ensembl_gene_id)],
            by = c("cohort", "protein_assay"))
cs[, `:=`(chr = as.character(chr), is_cis = as.logical(is_cis))]
# one row per member: overlapping windows can list a variant twice within a set
cs <- cs[, .(pip = max(pip), chr = chr[1L], pos = pos[1L], is_cis = is_cis[1L]),
         by = .(fluid, ensembl_gene_id, protein_assay, cs_id, snp)]
cs[, cs_key := paste(fluid, ensembl_gene_id, cs_id, sep = "|")]

spec <- fread(file.path(OUT_DIR, "03_fluid_specificity_credible_sets.csv"),
              select = c("cs_key", "shared"))
cs <- merge(cs, spec, by = "cs_key", all.x = TRUE)
cs[, class := fcase(is.na(shared), NA_character_, shared, "Shared",
                    fluid == "CSF", "CSF-specific", default = "Plasma-specific")]
cs[, type := fifelse(is_cis, "cis", "trans")]
n_sig <- uniqueN(cs$cs_key)
message("signals: ", n_sig, "; members: ", nrow(cs))

# ---- VEP consequence ----------------------------------------------------------
vep <- unique(fread(VEP_TSV, select = c("SNP", "Consequence")), by = "SNP")
cs <- merge(cs, vep[, .(snp = SNP, csq = csq_class(Consequence))], by = "snp", all.x = TRUE)
miss <- cs[is.na(csq), sum(pip)] / cs[, sum(pip)]
message(sprintf("PIP mass without a VEP record: %.2f%%", 100 * miss))
if (miss > 0.01) stop("more than 1% of PIP mass has no VEP record; rerun code/03_pqtl_mapping")
cs[is.na(csq), csq := "Unannotated"]

# Each signal's PIP is normalised to 1, so sets with total PIP < 1 do not weigh less.
cs[, w := pip / sum(pip), by = cs_key]
csq_tab <- function(d, by) {
  o <- d[, .(pip = sum(w)), by = c(by, "csq")]
  o[, `:=`(share = pip / sum(pip), n_signals = round(sum(pip))), by = by]
  o[]
}
csq_out <- rbind(
  csq_tab(cs, c("fluid", "type"))[, set := "all signals"],
  csq_tab(cs[!is.na(class)], c("class", "type"))[, set := "proteins assayed in both fluids"],
  fill = TRUE)
fwrite(csq_out, file.path(OUT_DIR, "03b_consequence.csv"))

# Per signal: PIP mass on protein-altering variants. A coding change alters the protein
# wherever it is measured, so it is the class expected to be shared by both fluids;
# shared versus each single-fluid class by Wilcoxon rank-sum.
per_sig <- cs[, .(pav_pip = sum(w[csq == "Protein-altering"])),
              by = .(cs_key, fluid, type, class)]
pav_test <- rbindlist(lapply(c("cis", "trans"), function(t) rbindlist(lapply(
  c("CSF-specific", "Plasma-specific"), function(k) {
    a <- per_sig[type == t & class == "Shared", pav_pip]
    b <- per_sig[type == t & class == k, pav_pip]
    data.table(type = t, versus = k, n_shared = length(a), n_other = length(b),
               mean_shared = mean(a), mean_other = mean(b),
               p = wilcox.test(a, b)$p.value)
  }))))
fwrite(pav_test, file.path(OUT_DIR, "03b_pav_test.csv"))
print(pav_test, digits = 3)

# ---- Roadmap chromatin --------------------------------------------------------
meta <- fread(file.path(ROADMAP_DIR, "EID_metadata.tab"))
cs <- cs[!(chr == MHC$chr & pos >= MHC$start & pos <= MHC$end)]
top <- cs[order(-pip), .SD[1L], by = cs_key][, .(cs_key, chr, pos)]
win <- GRanges(paste0("chr", top$chr), IRanges(pmax(1, top$pos - HALF_WIN), top$pos + HALF_WIN))
mem <- GRanges(paste0("chr", cs$chr), IRanges(cs$pos, width = 1L))
sig_idx <- match(cs$cs_key, top$cs_key)

obs <- exp <- matrix(0, nrow(top), nrow(meta), dimnames = list(top$cs_key, meta$EID))
for (j in seq_len(nrow(meta))) {
  f <- file.path(ROADMAP_DIR, paste0(meta$EID[j], "_15_coreMarks_mnemonics.bed.gz"))
  b <- fread(f, header = FALSE, col.names = c("chr", "start", "end", "state"))[state %in% ACTIVE]
  act <- reduce(GRanges(b$chr, IRanges(b$start + 1L, b$end)))   # BED is 0-based
  obs[, j] <- rowsum(cs$w * overlapsAny(mem, act), sig_idx, reorder = TRUE)[, 1L]
  h <- findOverlaps(win, act)
  cov <- width(pintersect(win[queryHits(h)], act[subjectHits(h)]))
  exp[, j] <- tabulate(rep(queryHits(h), cov), nbins = length(win)) / width(win)
}
stopifnot(all(exp > 0 | obs == 0))

lab <- unique(cs[, .(cs_key, fluid, type, class)])[match(top$cs_key, cs_key)]
groups <- rbind(
  lab[, .(cs_key, set = paste(fluid, type, sep = " "))],
  lab[!is.na(class), .(cs_key, set = paste(class, type, sep = " "))])
enr <- rbindlist(lapply(split(groups$cs_key, groups$set), function(k) {
  i <- match(k, top$cs_key)
  o <- colSums(obs[i, , drop = FALSE]); e <- colSums(exp[i, , drop = FALSE])
  bt <- replicate(N_BOOT, {
    s <- sample(i, replace = TRUE)
    log(colSums(obs[s, , drop = FALSE]) / colSums(exp[s, , drop = FALSE]))
  })
  data.table(EID = meta$EID, n_signals = length(i), observed = o, expected = e,
             enrichment = o / e, lo = exp(apply(bt, 1, quantile, 0.025)),
             hi = exp(apply(bt, 1, quantile, 0.975)))
}), idcol = "set")
enr <- merge(enr, meta[, .(EID, GROUP, STD_NAME, ANATOMY)], by = "EID")
fwrite(enr[order(set, EID)], file.path(OUT_DIR, "03b_chromatin_enrichment.csv"))

# CSF-specific versus plasma-specific, cis: log enrichment ratio per epigenome, with a
# two-sided bootstrap P (share of resamples on the other side of 0, doubled).
contrast <- function(type_) {
  a <- match(lab[class == "CSF-specific" & type == type_, cs_key], top$cs_key)
  b <- match(lab[class == "Plasma-specific" & type == type_, cs_key], top$cs_key)
  lr <- function(ia, ib) log(colSums(obs[ia, , drop = FALSE]) / colSums(exp[ia, , drop = FALSE])) -
                         log(colSums(obs[ib, , drop = FALSE]) / colSums(exp[ib, , drop = FALSE]))
  est <- lr(a, b)
  bt <- replicate(N_BOOT, lr(sample(a, replace = TRUE), sample(b, replace = TRUE)))
  data.table(EID = meta$EID, type = type_, n_csf = length(a), n_plasma = length(b),
             log_ratio = est, lo = apply(bt, 1, quantile, 0.025),
             hi = apply(bt, 1, quantile, 0.975),
             p = pmin(1, 2 * pmin(rowMeans(bt <= 0), rowMeans(bt >= 0))))
}
con <- rbind(contrast("cis"), contrast("trans"))
con[, fdr := p.adjust(p, "BH"), by = type]
con <- merge(con, meta[, .(EID, GROUP, STD_NAME)], by = "EID")
fwrite(con[order(type, -log_ratio)], file.path(OUT_DIR, "03b_chromatin_contrast.csv"))

message("Consequence, cis signals of proteins assayed in both fluids (share of PIP):")
print(dcast(csq_out[set != "all signals" & type == "cis"], csq ~ class, value.var = "share"),
      digits = 2)
message("Mean enrichment by Roadmap group, cis:")
print(dcast(enr[grepl("cis$", set), .(enr = round(mean(enrichment), 2)), by = .(set, GROUP)],
            GROUP ~ set, value.var = "enr"))
message("CSF- vs plasma-specific cis, mean log ratio by group:")
print(con[type == "cis", .(log_ratio = round(mean(log_ratio), 2), n_fdr05 = sum(fdr < 0.05)),
          by = GROUP][order(-log_ratio)])
message("03b_variant_annotation.R done: ", n_sig, " signals")
