#!/usr/bin/env Rscript
# trans-pQTL hotspots (candidate master-regulator loci) for BF2.
# Fine-mapped *trans* credible sets (is_cis == FALSE) are clustered into loci by
# merging signals within 1 Mb; a locus is a hotspot when it acts in trans on >= 5
# distinct target proteins. Each locus lead is annotated with a regulator using a
# 3-tier priority: the gene the lead variant sits in (VEP), else a co-located
# measured protein whose own cis-pQTL colocates, else the nearest gene (VEP).
# VEP genes are reused from the main pipeline's per-SNP annotation
# (code/03_pqtl_mapping -> annotation/protein/variant_annotation.tsv.gz,
# SYMBOL/NEAREST columns); this script does not run VEP itself.
# Outputs:
#   01b_hotspot_loci.csv     - all trans loci (one row per locus)
#   01b_hotspot_summary.csv  - hotspots only (>= 5 proteins), ranked
#   01b_hotspot_replication.csv - BF1 replication of the trans signals in each hotspot
# The hotspot panel drawn from these tables (01_pqtl_architecture.png c) is in
# lib/hotspot_panels.R and is drawn by 01c_plot.R.
# Usage:
#   Rscript code/09_statistical_analysis/01b_hotspot_data.R [loci_csv] [summary_csv]

suppressPackageStartupMessages(library(data.table))

DEFAULT_OUT <- "output"

args <- commandArgs(trailingOnly = TRUE)
loci_csv <- if (length(args) >= 1) args[1] else
  file.path(DEFAULT_OUT, "01b_hotspot_loci.csv")
summary_csv <- if (length(args) >= 2) args[2] else
  file.path(DEFAULT_OUT, "01b_hotspot_summary.csv")

COHORT <- "bf2"
WINDOW <- 1e6          # merge trans signals within 1 Mb into one locus
INFO_DIR <- "info/bf_protein"
FINEMAP_CSV <- file.path(INFO_DIR, "finemap_lead_snps.csv")
KEPT_CSV <- file.path(DEFAULT_OUT, "00_kept_assays.csv")
REPL_TSV <- file.path(INFO_DIR, "replication", "bf1_targeted_replication.tsv.gz")
CS_REPL_TSV <- file.path(INFO_DIR, "replication", "cs_replication.tsv.gz")

source("code/09_statistical_analysis/lib/cs_signals.R")
source("code/09_statistical_analysis/lib/common.R")   # assay_platform, assay_fluid, HOTSPOT_MIN

autosomes <- as.character(1:22)

# --- Load fine-mapped credible sets, restrict to deduped analyzed assays --------
fm <- fread(FINEMAP_CSV, select = c(
  "cohort", "protein_assay", "ensg", "gene", "chr", "pos", "snp", "cs_id",
  "max_pip", "is_cis"
))
# One signal fine-mapped in several overlapping windows would otherwise be counted
# once per window when the trans targets of a locus are tallied.
fm <- distinct_signals(fm[cohort == COHORT])
fm[, is_cis := as.logical(is_cis)]
stopifnot(!anyNA(fm$is_cis))

# Deduped allowlist (one highest-N assay per protein per fluid); brings the fluid
# label so protein counts are not inflated by repeated assays.
kept <- fread(KEPT_CSV, select = c("cohort", "protein_assay", "fluid"))
fm <- merge(fm, kept, by = c("cohort", "protein_assay"))  # inner join = filter

fm[, chr := as.character(chr)]
fm <- fm[chr %in% autosomes]
fm[, chr := as.integer(chr)]
stopifnot(nrow(fm) > 0)

# --- Cluster trans signals into loci (1 Mb single-linkage merge) ----------------
trans <- fm[is_cis == FALSE]
stopifnot(nrow(trans) > 0)
setorder(trans, chr, pos)
trans[, new_locus := chr != shift(chr, fill = -1L) |
        (pos - shift(pos, fill = 0L)) > WINDOW]
trans[, locus_id := cumsum(new_locus)]

cis <- fm[is_cis == TRUE]

# Candidate local regulator: distinct cis-target genes whose cis credible set
# lies within WINDOW of the locus's trans lead SNP.
cis_gene_near <- function(lchr, lpos) {
  g <- cis[chr == lchr & abs(pos - lpos) <= WINDOW, unique(gene)]
  if (length(g)) paste(sort(g), collapse = ";") else ""
}

# Representative variant: the lead shared by the most distinct target proteins (each
# protein once, whatever its fluid or assay), ties by PIP. It is the most frequently
# shared lead, not one that carries most of the locus: a large hotspot (ABO) holds
# several distinct leads. The highest-PIP lead of a single protein, used until
# 2026-10-24, moved with the trans threshold (ABO's fell 30 kb outside ABO). Every
# signal keeps its own lead, credible set and PIP; only the locus summary uses this.
loci <- trans[, {
  shared <- .SD[, .(n_prot = uniqueN(ensg), pip = max(max_pip), pos = pos[1L]), by = snp]
  setorder(shared, -n_prot, -pip)
  lead <- shared[1L]
  ug <- sort(unique(gene))
  .(
    chr = chr[1],
    start = min(pos),
    end = max(pos),
    span_kb = round((max(pos) - min(pos)) / 1e3, 1),
    n_trans_cs = .N,
    # per fluid: a hotspot is called in one fluid, so its signals are counted in that fluid
    n_trans_cs_csf = sum(fluid == "CSF"),
    n_trans_cs_plasma = sum(fluid == "Plasma"),
    n_proteins = uniqueN(ensg),
    n_csf = uniqueN(ensg[fluid == "CSF"]),
    n_plasma = uniqueN(ensg[fluid == "Plasma"]),
    lead_snp = lead$snp,
    lead_pos = lead$pos,
    lead_n_proteins = lead$n_prot,
    target_genes = paste(ug, collapse = ";"),
    # per fluid, for the pathway enrichment of each fluid's hotspot (01b_hotspot_enrichment.R)
    target_genes_csf = paste(sort(unique(gene[fluid == "CSF"])), collapse = ";"),
    target_genes_plasma = paste(sort(unique(gene[fluid == "Plasma"])), collapse = ";")
  )
}, by = locus_id]

loci[, candidate_cis_gene := mapply(cis_gene_near, chr, lead_pos)]
# Hotspots are called separately in each fluid, not on the pooled target count. A
# pooled call is not readable off the figure: a locus with 4 CSF and 4 plasma targets
# would be flagged while never reaching the threshold in either panel, and ABO would
# be drawn orange in the CSF facet at y = 4. Per-fluid, an orange point always means
# the y-value next to it clears HOTSPOT_MIN.
loci[, `:=`(is_hotspot_csf = n_csf >= HOTSPOT_MIN,
            is_hotspot_plasma = n_plasma >= HOTSPOT_MIN)]
loci[, is_hotspot := is_hotspot_csf | is_hotspot_plasma]
loci[, hotspot_fluid := fcase(
  is_hotspot_csf & is_hotspot_plasma, "Both",
  is_hotspot_csf, "CSF",
  is_hotspot_plasma, "Plasma",
  default = "")]
setorder(loci, -n_proteins, chr, start)

# --- Regulator annotation, reusing the postgwas VEP pipeline --------------------
# code/03_pqtl_mapping VEP-annotates every lead / credible-set SNP; we reuse its
# per-SNP gene symbols (SYMBOL = gene the lead sits in, e.g. ABO/GNPTAB/APOE;
# NEAREST = nearest gene for intergenic leads) instead of re-running VEP here.
# candidate_cis_gene (a measured protein whose own cis-pQTL colocates) stays as
# independent, mechanistic supporting evidence.
ANNO_TSV <- file.path(GWAS_BASE, COHORT,
                      "annotation/protein/variant_annotation.tsv.gz")
# Without VEP the regulator labels silently fall back to the cis-pQTL tier and
# change the hotspot labels in 01_pqtl_architecture.png b, so a missing or symbol-less annotation is an error.
if (!file.exists(ANNO_TSV) ||
    !all(c("SNP", "SYMBOL", "NEAREST", "Extra") %in% names(fread(ANNO_TSV, nrows = 0)))) {
  stop("VEP SYMBOL/NEAREST not available (re-run code/03_pqtl_mapping with ",
       "--symbol): ", ANNO_TSV)
}
anno <- unique(fread(ANNO_TSV, select = c("SNP", "SYMBOL", "NEAREST", "Extra"),
                     colClasses = "character"), by = "SNP")
# VEP records where the symbol came from. A Clone_based_* source means the "gene"
# is an unannotated transcript named after the BAC clone it was sequenced on
# (RP11-649A16.1 and the like) -- a coordinate, not a candidate regulator -- so the
# symbol is dropped and the locus falls through to the next tier below. Read from
# SYMBOL_SOURCE rather than pattern-matched on the name, which is what VEP means.
anno[, sym_source := sub(".*SYMBOL_SOURCE=([^;]*).*", "\\1", Extra)]
anno[grepl("^Clone_based", sym_source), SYMBOL := ""]
# A pseudogene is not a candidate regulator (the ABO locus was once named after LCN1P2),
# so the gene a lead sits in is dropped when VEP calls its transcript a pseudogene.
if (!any(grepl("BIOTYPE=", anno$Extra)))
  stop("VEP output has no BIOTYPE (rerun code/03_pqtl_mapping with --biotype): ", ANNO_TSV)
anno[, biotype := fifelse(grepl("BIOTYPE=", Extra), sub(".*BIOTYPE=([^;]*).*", "\\1", Extra), "")]
message("VEP symbols dropped as pseudogenes: ", anno[SYMBOL != "" & grepl("pseudogene", biotype), .N])
anno[grepl("pseudogene", biotype), SYMBOL := ""]
anno[, c("Extra", "sym_source", "biotype") := NULL]
setnames(anno, c("SNP", "SYMBOL", "NEAREST"),
         c("lead_snp", "vep_symbol", "vep_nearest"))
for (col in c("vep_symbol", "vep_nearest"))
  anno[is.na(get(col)) | get(col) == "NA", (col) := ""]

# Locus gene, chosen apart from the representative variant: the gene (VEP SYMBOL) that
# the trans leads of the most distinct target proteins sit in, so a single lead in a
# pseudogene next to ABO cannot rename the locus. Ties go to the representative
# variant's gene, then alphabetically. A positional label, not a causal assignment.
votes <- merge(trans[, .(locus_id, lead_snp = snp, ensg)], anno[, .(lead_snp, vep_symbol)],
               by = "lead_snp", all.x = TRUE)
votes <- votes[!is.na(vep_symbol) & vep_symbol != "", .(n = uniqueN(ensg)), by = .(locus_id, vep_symbol)]
loci <- merge(loci, anno, by = "lead_snp", all.x = TRUE, sort = FALSE)
for (col in c("vep_symbol", "vep_nearest"))
  loci[is.na(get(col)), (col) := ""]
votes <- merge(votes, loci[, .(locus_id, rep_symbol = vep_symbol)], by = "locus_id")
votes[, is_rep := vep_symbol == rep_symbol]
setorder(votes, locus_id, -n, -is_rep, vep_symbol)
locus_gene <- votes[, .(locus_gene = vep_symbol[1L], locus_gene_n = n[1L]), by = locus_id]
loci <- merge(loci, locus_gene, by = "locus_id", all.x = TRUE)
loci[is.na(locus_gene), `:=`(locus_gene = "", locus_gene_n = 0L)]
setorder(loci, -n_proteins, chr, start)
message("Merged VEP genes from: ", ANNO_TSV)

# Regulator label priority: locus gene (above) > co-located cis-pQTL protein >
# nearest gene of the representative variant. Kept in `regulator` for the figure/table,
# with its source in `regulator_source`.
loci[, regulator_source := fifelse(locus_gene != "", "VEP gene",
                            fifelse(candidate_cis_gene != "", "cis-pQTL", "nearest"))]
loci[, regulator := fifelse(
  locus_gene != "", locus_gene,
  fifelse(candidate_cis_gene != "", tstrsplit(candidate_cis_gene, ";")[[1]],
          vep_nearest))]
hs_votes <- merge(loci[is_hotspot == TRUE, .(locus_id, regulator, lead_snp)],
                  votes[, .(locus_id, gene = vep_symbol, n)], by = "locus_id")[order(locus_id, -n)]
hs_votes[, locus_id := NULL]
loci[, locus_id := NULL]

fwrite(loci, loci_csv)
message("Saved: ", loci_csv, " (", nrow(loci), " trans loci)")

hotspots <- loci[is_hotspot == TRUE]
fwrite(hotspots, summary_csv)
# The runner-up genes behind each hotspot label, so a close vote is visible in the log.
message("Hotspot gene votes (distinct target proteins whose lead sits in the gene; top 3):")
print(hs_votes[, head(.SD, 3), by = .(regulator, lead_snp)])
message("Saved: ", summary_csv, " (", nrow(hotspots), " hotspots)")

# --- BF1 replication of the trans signals inside each hotspot -------------------
# Hotspots are built entirely from trans signals, the stratum that replicates least
# well overall, so each hotspot gets its own replication rate. Strict exact-assay trans tests are matched to a hotspot by
# genomic position of the BF2 sentinel.
# Same two-tier matching as the sentinel replication (01a_pqtl_data.R): exact assay where available,
# otherwise the same gene in the same fluid on a different assay version. Restricting
# this table to exact assays would leave every CSF-only hotspot untestable, because
# BF1 CSF was run on different Olink panels, and would use a different rule from the
# sentinel replication panel. Replication is called at concordant direction plus BF1
# P < 0.05, as there: these sentinels are pre-specified in BF2, so no multiplicity
# correction applies.
# The verdict is read, not recomputed -- see the note in 01a. The raw table is
# still what supplies fluid and is_cis, which the hotspot table splits on.
cs_rep <- fread(CS_REPL_TSV, select = c("trait", "cs_id", "cs_lead_snp", "is_representative",
                                        "bf1_status"))
repl <- fread(REPL_TSV)
repl[, trait_key := paste(cohort, bf2_assay, sep = "__")]
# Through cs_id, not SNP: a proxy-rescued sentinel carries the proxy variant in
# SNP, so an SNP join drops it. Same fix as in 01a.
repl <- merge(repl, unique(cs_rep[, .(trait_key = trait, cs_id, cs_lead_snp)]),
              by = c("trait_key", "cs_id"), all.x = TRUE)
# Distinct signals only, as in 01a: a credible set resolved again from an
# overlapping window would otherwise count twice in a hotspot's denominator.
rep_leads <- unique(cs_rep[is_representative == TRUE, .(trait_key = trait, cs_lead_snp)])
rep_leads <- merge(rep_leads, reported_leads()[, .(trait_key = trait, cs_lead_snp)],
                   by = c("trait_key", "cs_lead_snp"))   # reported signals only, as in 01a
repl <- merge(repl, rep_leads, by = c("trait_key", "cs_lead_snp"))
status <- unique(cs_rep[, .(trait_key = trait, cs_lead_snp, bf1_status)],
                 by = c("trait_key", "cs_lead_snp"))
repl <- merge(repl, status, by = c("trait_key", "cs_lead_snp"), all.x = TRUE)
if (anyNA(repl$bf1_status)) {
  stop(sum(is.na(repl$bf1_status)), " submitted sentinels have no verdict in ",
       basename(CS_REPL_TSV), " -- rerun code/05_replication/04_credible_set_status.R")
}
repl <- repl[bf1_status %in% c("replicated", "tested, not replicated")]
repl[, replicated := bf1_status == "replicated"]
repl <- repl[is_cis == FALSE]
repl[, fluid_lab := fifelse(fluid == "csf", "CSF", "Plasma")]
# One row per hotspot *per fluid*, matching how the hotspot map defines them: a
# locus is a hotspot in CSF, in plasma, or in both, and a locus that is a hotspot in
# both is two rows here. Pooling the two fluids into one rate would ask a question the
# hotspot definition never asks, and would let plasma signals decide the rate of a
# CSF-only hotspot.
hs_long <- rbind(
  loci[is_hotspot_csf == TRUE][, .(lead_snp, chr, start, end, regulator,
                                   fluid = "CSF", n_proteins = n_csf)],
  loci[is_hotspot_plasma == TRUE][, .(lead_snp, chr, start, end, regulator,
                                      fluid = "Plasma", n_proteins = n_plasma)]
)
hs_repl <- hs_long[, {
  # Renamed first: a bare chr/start/end inside repl[...] would resolve to a repl
  # column of that name if one ever appeared (the data.table trap in CLAUDE.md).
  h_fluid <- .BY$fluid; h_chr <- chr; h_start <- start; h_end <- end
  s <- repl[fluid_lab == h_fluid & CHR == h_chr & BP >= h_start & BP <= h_end]
  .(n_proteins = n_proteins,
    n_testable = nrow(s),
    n_repl     = sum(s$replicated, na.rm = TRUE))
}, by = .(lead_snp, chr, regulator, fluid)]
hs_repl[, fluid := factor(fluid, levels = c("CSF", "Plasma"))]
hs_repl[, frac_repl := fifelse(n_testable > 0, n_repl / n_testable, NA_real_)]
# Row order is the panel's order (lib/hotspot_panels.R sorts the same way).
setorder(hs_repl, fluid, frac_repl, n_testable, na.last = FALSE)
repl_csv <- file.path(dirname(summary_csv), "01b_hotspot_replication.csv")
fwrite(hs_repl[, .(lead_snp, chr, regulator, fluid, n_proteins, n_testable, n_repl,
                   frac_repl)], repl_csv)

stopifnot(file.exists(loci_csv), file.exists(summary_csv), file.exists(repl_csv))
message("01b_hotspot_data.R done: ", nrow(hotspots), " hotspots across ",
        uniqueN(hotspots$chr), " chromosomes")
