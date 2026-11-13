#!/usr/bin/env Rscript
# Cell-type-of-origin: within brain, WHICH cell type carries a protein's cis
# regulation? 04_protein_tissue.R settles tissue (does CSF map to brain, plasma to
# periphery?) using bulk resources; this is the next question down, and it is the one
# bulk tissue cannot answer at all.
#
# Resource: SingleBrain single-nucleus eQTL, 7 major brain cell types (Nat Genet 2026,
# doi:10.1038/s41588-026-02541-x). MetaBrain bulk brain is loaded alongside it purely
# as the comparison arm -- what does single-nucleus recover that bulk misses?
#
# Scope is identical to 04 by construction (lib/eqtl_coloc_scope.R): same assays, same
# cis + cognate filter, same PP4 cutoff. Numbers here and in 04 are comparable only
# because that code is shared rather than copied.
#
# Run from the repo root:
#   Rscript code/09_statistical_analysis/04b_protein_celltype.R
suppressPackageStartupMessages({
  library(data.table)
})

source("code/09_statistical_analysis/lib/common.R")   # assay_platform, assay_fluid, CHR_LEN_HG19
source("code/09_statistical_analysis/lib/eqtl_coloc_scope.R")

OUT <- "output"
out_summary   <- file.path(OUT, "04b_celltype_summary.csv")
out_spec      <- file.path(OUT, "04b_celltype_specificity.csv")
out_vs_bulk   <- file.path(OUT, "04b_celltype_vs_bulk.csv")
out_evidence  <- file.path(OUT, "04b_celltype_pair_evidence.csv.gz")

annot <- load_assay_annot()
cb <- rbindlist(list(
  read_coloc(COLOC_RESULT("sb"), "SingleBrain"),
  read_coloc(COLOC_RESULT("mb"), "MetaBrain")
), use.names = TRUE)
cb <- apply_cis_cognate_scope(cb, annot)

# Denominator. status "ok" and "no_sig_coloc" are the tests that actually ran, i.e.
# both sides had a credible set. "no_region" (no protein cis CS) and "no_sb_cs" (the
# eQTL side had none) were never tested, so counting them would report "not measured"
# as "no colocalization" -- the single easiest way to misstate this analysis.
tested <- cb[resource == "SingleBrain" & status %chin% c("ok", "no_sig_coloc"),
             unique(exposure_id)]
sig <- cb[coloc_retained == TRUE & PP4 > COLOC_PP4]
sb  <- sig[resource == "SingleBrain"]
message(sprintf("SingleBrain: %d assays tested, %d with a colocalization (%.1f%%)",
                length(tested), uniqueN(sb$exposure_id),
                100 * uniqueN(sb$exposure_id) / length(tested)))

# ---- 1. Per cell type ------------------------------------------------------------
summary_dt <- sb[, .(
  coloc_events = .N,
  assays       = uniqueN(exposure_id),
  genes        = uniqueN(GENE),
  PP4_median   = round(median(PP4), 3)
), by = .(cell_type = tissue)]
summary_dt[, cell_type_label := CELLTYPE_LAB[cell_type]]
# Split by fluid too: the CSF/plasma contrast from 04 should not silently vanish here.
by_fluid <- dcast(sb[, .(assays = uniqueN(exposure_id)), by = .(tissue, fluid)],
                  tissue ~ fluid, value.var = "assays", fill = 0L)
summary_dt <- merge(summary_dt, by_fluid, by.x = "cell_type", by.y = "tissue", all.x = TRUE)
setorder(summary_dt, -coloc_events)
fwrite(summary_dt, out_summary)

# ---- 2. Cell-type specificity ----------------------------------------------------
# The headline: in how many of the 7 cell types does a given (assay, gene) colocalize?
# n_cell_types == 1 is the result bulk brain cannot produce.
spec <- sb[, .(n_cell_types = uniqueN(tissue),
               cell_types = paste(sort(unique(tissue)), collapse = "|"),
               PP4_max = round(max(PP4), 3)),
           by = .(exposure_id, gene, fluid, ENSG)]
setorder(spec, -n_cell_types, gene)
fwrite(spec, out_spec)

spec_tab <- spec[, .(pairs = .N), by = n_cell_types][order(n_cell_types)]
spec_tab[, pct := round(100 * pairs / sum(pairs), 1)]
message("Cell-type breadth per (assay, gene):")
print(spec_tab)

# Which cell type carries the single-cell-type-specific pairs?
spec1 <- spec[n_cell_types == 1][, .(specific_pairs = .N), by = .(cell_type = cell_types)]
summary_dt <- merge(summary_dt, spec1, by = "cell_type", all.x = TRUE)
summary_dt[is.na(specific_pairs), specific_pairs := 0L]
fwrite(summary_dt, out_summary)

# ---- 3. Single-nucleus vs bulk brain ---------------------------------------------
# Per (assay, gene): does MetaBrain bulk also support it? A pair found only in
# SingleBrain is the concrete argument for the resource -- but note the two differ in
# more than resolution (different donors, different eGene power), so "bulk missed it"
# is a statement about these two datasets, not proof of a cell-type-restricted effect.
sb_pairs <- unique(sb[, .(exposure_id, gene, ENSG, fluid)])
mb_pairs <- unique(sig[resource == "MetaBrain", .(exposure_id, gene, ENSG)])
sb_pairs[, in_bulk := paste(exposure_id, ENSG) %in% mb_pairs[, paste(exposure_id, ENSG)]]
vs_bulk <- sb_pairs[, .(pairs = .N), by = .(fluid, in_bulk)]
fwrite(vs_bulk, out_vs_bulk)
message("SingleBrain pairs also seen in MetaBrain bulk:")
print(dcast(vs_bulk, fluid ~ in_bulk, value.var = "pairs", fill = 0L))

fwrite(sb[, .(exposure_id, gene, fluid, cell_type = tissue, ENSG, GENE,
              cs_id, lead_snp, o_lead, chr, pos, PP4)], out_evidence)

# No figure (2026-11-13): the per-cell-type counts are the SingleBrain block of
# 04_supp_coloc.png (CSF and plasma), and the text gives the single-cell-type share.

for (f in c(out_summary, out_spec, out_vs_bulk, out_evidence)) {
  if (!file.exists(f)) stop("Expected output missing: ", f)
  message("Saved: ", f)
}
