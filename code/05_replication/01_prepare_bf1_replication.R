#!/usr/bin/env Rscript
library(data.table)

source(path.expand("code/lib/common.R"))   # REPO_ROOT, GWAS_BASE
repo_root = REPO_ROOT

args = commandArgs(trailingOnly = TRUE)
out_dir = if (length(args) >= 1L) path.expand(args[1]) else file.path(repo_root, "info/bf_protein/replication")
bf1_clean_dir = if (length(args) >= 2L) path.expand(args[2]) else file.path(GWAS_BASE, "bf1/clean/protein")
bf2_clean_dir = if (length(args) >= 3L) path.expand(args[3]) else file.path(GWAS_BASE, "bf2/clean/protein")

annot_file = file.path(repo_root, "info/bf_protein/assay_annot.csv")
lead_file = file.path(repo_root, "info/bf_protein/finemap_lead_snps.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fluid_from_tag = function(x) {
  fifelse(grepl("csf", x, ignore.case = TRUE), "csf",
          fifelse(grepl("plasma", x, ignore.case = TRUE), "plasma", NA_character_))
}

annot = fread(annot_file)
# Fine-mapped credible-set leads, not clump leads. The manuscript reports credible
# sets, so the replication set has to be the same unit: for cis signals only ~45% of
# credible-set leads are the clump lead, and a clump can carry several credible sets.
# hc_p5e8 keeps the canonical 5e-8 discovery threshold.
leads = fread(lead_file)[cohort == "bf2" & hc_p5e8 == TRUE]
setnames(leads, c("snp", "chr", "pos", "p", "ensg"),
         c("SNP", "CHR", "BP", "P", "ensembl_gene_id"))
leads = leads[, .(cohort, protein_assay, SNP, CHR, BP, P, gene, ensembl_gene_id,
                  is_cis, cs_id, cs_size, cs_purity, max_pip)]

# Reconstruct the manuscript discovery set. In CSF, choose the largest-N BF2
# assay per ENSG among assays with N > 1,000. Plasma Explore HT assay identities
# are retained (including distinct NPPB and NT-proBNP analytes for the same ENSG).
bf2_csf = annot[
  cohort == "bf2" & grepl("csf", tag, ignore.case = TRUE) & N > 1000 &
    !is.na(ensembl_gene_id)
]
setorder(bf2_csf, ensembl_gene_id, -N, protein_assay)
bf2_csf = bf2_csf[!duplicated(ensembl_gene_id)]
bf2_csf[, fluid := "csf"]

bf2_plasma = annot[
  cohort == "bf2" & tag == "ExploreHT_o7plasma" & N > 1000 &
    !is.na(ensembl_gene_id)
]
bf2_plasma[, fluid := "plasma"]

bf2_selected = rbindlist(list(bf2_csf, bf2_plasma), use.names = TRUE, fill = TRUE)[,
  .(bf2_assay = protein_assay, bf2_tag = tag, bf2_annotation_n = N, fluid)
]
setnames(leads, c("protein_assay", "P"), c("bf2_assay", "bf2_lead_p"))
candidates = merge(leads, bf2_selected, by = "bf2_assay", all = FALSE)

observed = candidates[, .(leads = .N, cis = sum(is_cis), trans = sum(!is_cis)), by = fluid]
print(observed)
if (nrow(observed) != 2L || any(observed$leads < 1000L)) {
  stop("BF2 credible-set discovery set looks wrong")
}

bf1_files = list.files(bf1_clean_dir, pattern = "^bf1__.*\\.tsv\\.gz$")
bf1_current = sub("\\.tsv\\.gz$", "", sub("^bf1__", "", bf1_files))
bf1_annot = unique(annot[
  cohort == "bf1" & protein_assay %chin% bf1_current & !is.na(ensembl_gene_id),
  .(bf1_assay = protein_assay, bf1_tag = tag, ensembl_gene_id, bf1_annotation_n = N)
])
bf1_annot[, fluid := fluid_from_tag(bf1_tag)]
bf1_annot = bf1_annot[!is.na(fluid)]

# All same-gene/same-fluid assay mappings are considered, then one is selected
# deterministically. Exact assay identity defines the prespecified primary
# replication analysis. A greatest-N same-gene fallback is retained only as an
# explicitly exploratory cross-assay analysis.
candidates = merge(
  candidates,
  bf1_annot,
  by = c("ensembl_gene_id", "fluid"),
  all.x = TRUE,
  allow.cartesian = TRUE
)
candidates[, exact_assay := !is.na(bf1_assay) & bf2_assay == bf1_assay]
setorder(candidates, bf2_assay, SNP, -exact_assay, -bf1_annotation_n, bf1_assay, na.last = TRUE)
candidates = candidates[!duplicated(paste(bf2_assay, SNP))]
candidates[, assay_match_type := fcase(
  is.na(bf1_assay), "not_measured_same_fluid",
  exact_assay, "exact_assay",
  default = "same_gene_fluid"
)]
candidates[, `:=`(
  bf1_assay_available = !is.na(bf1_assay),
  strict_replication_eligible = exact_assay,
  exploratory_cross_assay_eligible = !exact_assay & !is.na(bf1_assay)
)]
setorder(candidates, fluid, bf2_assay, bf2_lead_p, SNP)
candidates[, replication_id := sprintf("BF1REP%05d", .I)]

candidates[, `:=`(snp_role = "lead", member_pip = max_pip, member_rank = 0L)]

# Proxy pool. When the credible-set lead is absent from the BF1 clean summary
# statistics the signal is not lost: any other member of the same credible set is
# in high LD with it, so the highest-PIP member present in BF1 tests the same
# signal. Members are carried here as extra rows so the extraction step picks them
# up unchanged; 03 selects at most one row per signal (lead first, then the
# highest-PIP available member) and flags which signals fell back to a proxy.
# Always on: the Methods describe this substitution, and the lead-only branch that
# used to stand here never ran, because the members file is tracked in git. So the
# file must exist and match this run's credible sets (01a reruns after every fine-mapping
# run). Sets are keyed by (assay, cs_id): a lead can head credible sets of several windows.
member_file = file.path(repo_root, "info/bf_protein/finemap_cs_members.csv.gz")
if (!file.exists(member_file)) stop("missing ", member_file, " -- run 01a_extract_cs_members.sbatch")
# Without this check a fine-mapping rerun would silently take its proxies from the
# previous run. Compare content, not mtimes -- git checkout/pull rewrites mtimes. 01a
# writes one set per (assay, cs_id) of the hc_p5e8 leads, so the key sets must match.
lead_keys = unique(leads[, paste(bf2_assay, cs_id)])
mem_keys = unique(fread(member_file, select = c("protein_assay", "cs_id"))[
  , paste(protein_assay, cs_id)])
n_miss = length(setdiff(lead_keys, mem_keys)); n_extra = length(setdiff(mem_keys, lead_keys))
if (n_miss || n_extra)
  stop(member_file, " does not match ", lead_file, " (", n_miss, " leads missing, ",
       n_extra, " stale sets) -- re-run 01a_extract_cs_members.sbatch")
members = fread(member_file)[member_rank > 1L]
# The candidate's own credible set, by cs_id (not by lead SNP, which several sets share).
proxies = merge(
  candidates[, .(replication_id, bf2_assay, cs_id)],
  members[, .(bf2_assay = protein_assay, cs_id, member_snp, member_pip, member_rank)],
  by = c("bf2_assay", "cs_id")
)
if (!nrow(proxies)) stop("No credible-set proxies matched the discovery signals")
proxies = merge(
  proxies[, .(replication_id, SNP = member_snp, member_pip, member_rank)],
  candidates[, !c("SNP", "member_pip", "member_rank", "snp_role"), with = FALSE],
  by = "replication_id"
)
proxies[, snp_role := "proxy"]
# A member that is itself another signal's lead must not be duplicated inside one assay.
proxies = proxies[!candidates[, .(bf2_assay, SNP)], on = .(bf2_assay, SNP)]
candidates = rbindlist(list(candidates, proxies), use.names = TRUE, fill = TRUE)
setorder(candidates, fluid, bf2_assay, bf2_lead_p, replication_id, member_rank)
message("Proxy pool: ", nrow(proxies), " member variants for ",
        uniqueN(proxies$replication_id), " of ",
        uniqueN(candidates$replication_id), " signals")

m_planned = candidates[snp_role == "lead", sum(strict_replication_eligible)]
message("Replication universe: candidates=", candidates[snp_role == "lead", .N],
        ", exact-assay planned=", m_planned)
if (m_planned < 1L) stop("No exact-assay tests in the replication universe")
candidates[, `:=`(
  correction_family_n = m_planned,
  bonferroni_alpha = 0.05 / m_planned
)]

# BF2 statistics for every candidate, also where BF1 lacks the assay (code review
# 2026-11-01, C60, kept): the manuscript's attrition check (01a_pqtl_data.R) compares the
# frequency and P of untested sentinels, the not-measured ones included.
bf2_requests = unique(candidates[, .(
  cohort = "bf2", assay = bf2_assay,
  clean_path = file.path(bf2_clean_dir, paste0("bf2__", bf2_assay, ".tsv.gz"))
)])
bf1_requests = unique(candidates[bf1_assay_available == TRUE, .(
  cohort = "bf1", assay = bf1_assay,
  clean_path = file.path(bf1_clean_dir, paste0("bf1__", bf1_assay, ".tsv.gz"))
)])
requests = rbindlist(list(bf2_requests, bf1_requests))
setorder(requests, cohort, assay)
requests[, request_id := .I]
if (anyDuplicated(requests[, .(cohort, assay)])) stop("Duplicated extraction request")
if (any(!file.exists(requests$clean_path))) {
  stop("Missing clean files for ", sum(!file.exists(requests$clean_path)), " requests")
}

candidate_file = file.path(out_dir, "bf1_replication_candidates.tsv.gz")
# An output, not a tmp/ manifest: 02 extracts by it and 03 checks completeness and
# freshness against it, so it is kept beside the candidates it was built from.
request_file = file.path(out_dir, "bf1_replication_extract_requests.tsv.gz")
fwrite(candidates, candidate_file, sep = "\t")
fwrite(requests, request_file, sep = "\t")

print(candidates[snp_role == "lead", .(
  signals = .N,
  strict_eligible = sum(strict_replication_eligible),
  exact_assay = sum(assay_match_type == "exact_assay"),
  cross_assay = sum(assay_match_type == "same_gene_fluid")
), by = fluid])
message("Prepared ", candidates[snp_role == "lead", .N], " BF2 signals, ", m_planned,
        " strict exact-assay BF1 tests, and ", nrow(requests), " extraction requests")
message("Bonferroni alpha fixed at ", signif(0.05 / m_planned, 5))
