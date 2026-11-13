#!/usr/bin/env Rscript
# Extract credible-set members (top-K by PIP) for the BF2 signals used in the BF1
# targeted replication: the proxy pool used when a credible-set lead is absent from the
# BF1 clean summary stats.
#
# Members come from the *.finemapped.tsv.gz tables, one set per (assay, cs_id) of the
# hc_p5e8 leads -- the reported credible set itself. They used to come from the stored
# SuSiE objects, matched by lead SNP across every window's fit: a re-clump left old fits
# behind (1,021 were found and removed on 2026-11-01), and a lead found from several
# windows took whichever fit gave it the highest PIP (code review 2026-11-01, C38/C52).
#
# Usage: Rscript 01a_extract_cs_members.R [finemapping_dir] [out_file]
suppressPackageStartupMessages(library(data.table))

source(path.expand("code/lib/common.R"))   # REPO_ROOT, GWAS_BASE
repo_root = REPO_ROOT

args = commandArgs(trailingOnly = TRUE)
fm_dir = if (length(args) >= 1L && nzchar(args[1])) path.expand(args[1]) else
  file.path(GWAS_BASE, "bf2/finemapping/protein")
out_file = if (length(args) >= 2L && nzchar(args[2])) path.expand(args[2]) else
  file.path(repo_root, "info/bf_protein/finemap_cs_members.csv.gz")
TOPK = 10L

lead_file = file.path(repo_root, "info/bf_protein/finemap_lead_snps.csv")
leads = unique(fread(lead_file)[cohort == "bf2" & hc_p5e8 == TRUE, .(protein_assay, cs_id, snp)])
if (anyDuplicated(leads[, .(protein_assay, cs_id)])) stop("two leads for one (assay, cs_id) in ", lead_file)
message("Signals: ", nrow(leads), " credible sets across ", uniqueN(leads$protein_assay), " BF2 assays")

members = rbindlist(lapply(unique(leads$protein_assay), function(a) {
  f = file.path(fm_dir, paste0("bf2__", a, ".finemapped.tsv.gz"))
  if (!file.exists(f)) stop("missing finemapped table: ", f)
  d = fread(f, select = c("cs_id", "lead_snp", "snp", "pip"))
  d = d[leads[protein_assay == a, .(cs_id)], on = "cs_id", nomatch = NULL]
  if (!nrow(d)) return(NULL)
  # The lead first (rank 1), then by PIP: a lead can tie on PIP with another member.
  d[, is_lead := snp == lead_snp]
  setorder(d, cs_id, -is_lead, -pip, snp)
  d[, `:=`(member_rank = seq_len(.N), cs_n_members = .N), by = cs_id]
  d[member_rank <= TOPK, .(protein_assay = a, cs_id, cs_lead_snp = lead_snp, member_snp = snp,
                           member_pip = pip, member_rank, cs_n_members)]
}))
if (!nrow(members)) stop("No credible-set members extracted")

# Every lead must find its set; rank 1 is the listed lead and has the set's top PIP.
got = unique(members[, .(protein_assay, cs_id)])
if (nrow(got) != nrow(leads)) stop(nrow(leads) - nrow(got), " leads have no credible set in their finemapped table")
chk = merge(leads, members[, .(top_pip = max(member_pip), lead = member_snp[member_rank == 1L],
                                lead_pip = member_pip[member_rank == 1L]), by = .(protein_assay, cs_id)],
            by = c("protein_assay", "cs_id"))
n_bad = chk[snp != lead | lead_pip < top_pip, .N]
if (n_bad) stop(n_bad, " sets whose rank-1 member is not the listed lead or not of top PIP")
stopifnot(!anyDuplicated(members, by = c("protein_assay", "cs_id", "member_rank")))

setorder(members, protein_assay, cs_id, member_rank)
fwrite(members, out_file)
message("Wrote ", nrow(members), " members for ", nrow(got), " credible sets -> ", out_file)
print(members[, .(sets = uniqueN(paste(protein_assay, cs_id))), by = .(has_proxy = member_rank > 1L)])
