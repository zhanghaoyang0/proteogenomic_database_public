#!/usr/bin/env Rscript
# The pQTL novelty scan (docs/pipeline.md): tier C and
# tier A for every BioFINDER-2 credible set, against the reference sentinel lists
# built by 07a_novelty_refs.R.
#
#   CSF    vs Cruchaga et al. CSF SomaScan 7k (N=3,506) and Western et al. 2026 CSF
#             Olink Explore HT (N=2,524), pooled: assayed or reported by either counts
#   Plasma vs UKB-PPP plasma Olink 3k        (N=54,219)
#
# Western 2026 lists index pQTLs only, so at its loci a secondary signal of ours can look
# unreported (trans, read per locus); cis is read per protein, where this does not matter.
#
# The two arms are not symmetric and must not be read as if they were. UKB-PPP is 34x
# our plasma N on the same platform family and used a stricter threshold, so plasma is
# a negative control on the matching pipeline, not a discovery arm. CSF is where a
# novelty claim can stand.
#
# Tiers, reported apart (pooling them is how novelty counts get inflated):
#   C  the reference never assayed this protein  -> coverage, not novelty
#   A  assayed, but no reference sentinel for this protein within 1 Mb at r2 >= 0.8
#      -> previously unreported
#   unresolved  no pair reaches r2 0.8, but at least one reference sentinel in the
#      window has no r2 (either SNP missing from 1000G EUR) -> not counted as A
#   -  otherwise known
# For cis, the main call is per protein's primary association (novelty_primary_cis, below),
# where B = protein reported, primary association not matched. The references' own
# summary statistics are not used here (07c flags B associations Western detected).
#
# Matching is deliberately UKB-PPP's own novelty rule: same protein, sentinel within
# 1 Mb, r2 >= 0.8 in 1000G EUR - the same panel the clumping used.
#
# Usage (Slurm, reads ~1,000 summary-statistic files on SLURM_CPUS_PER_TASK cores):
#   Rscript code/09_statistical_analysis/07b_novelty_scan.R

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/cs_signals.R")
source("code/09_statistical_analysis/lib/common.R")   # add_novelty

OUT_DIR  <- "output"
TMP_DIR  <- "tmp/proteogenomic_database"
S1       <- file.path(OUT_DIR, "S1_pqtl_finemapping.csv.gz")
KEPT     <- file.path(OUT_DIR, "00_kept_assays.csv")
ANNOT    <- "info/bf_protein/assay_annot.csv"
REF_DIR  <- path.expand("~/db/shared_data/pqtl_ref/parsed")
BFILE    <- path.expand("~/db/shared_data/ref_geno/1kg_v3/EUR")
WINDOW   <- 1e6   # +/- 1 Mb, UKB-PPP's definition
R2_MATCH <- 0.8   # UKB-PPP's LD threshold for "same signal"

for (f in c(S1, KEPT, ANNOT, file.path(REF_DIR, c("ukbppp_panel.csv", "ukbppp_sentinels.csv",
                                            "cruchaga_panel.csv", "cruchaga_sentinels.csv",
                                            "western_olink_panel.csv", "western_olink_sentinels.csv"))))
  if (!file.exists(f)) stop("missing input: ", f)
if (!file.exists(paste0(BFILE, ".bim"))) stop("missing LD panel: ", BFILE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- our credible-set sentinels ---------------------------------------------
s1 <- fread(S1)
sent <- s1[is_cs_lead == TRUE,
           .(cs_key, fluid, protein_assay, gene, ensembl_gene_id, is_cis, snp, chr, pos)]
if (uniqueN(sent$cs_key) != nrow(sent)) stop("more than one lead per credible set")

annot <- unique(fread(ANNOT)[cohort == "bf2", .(protein_assay, uniprot_id)])
sent <- merge(sent, annot, by = "protein_assay", all.x = TRUE)
sent[uniprot_id == "", uniprot_id := NA_character_]
message("credible sets: ", nrow(sent), "; without a UniProt accession: ",
        sent[is.na(uniprot_id), .N])

# ---- reference panels and sentinels -----------------------------------------
# Cruchaga aptamers and Western Olink complexes can tag more than one accession
# ("P29459|P29460"; 07a writes both with "|"); split so a protein counts as assayed if
# it appears in any of them.
split_uniprot <- function(d) {
  d <- copy(d)
  d[, uniprot_id := strsplit(as.character(uniprot_id), "[|]")]
  d[, .(uniprot_id = unlist(uniprot_id)), by = setdiff(names(d), "uniprot_id")]
}

cr_panel <- split_uniprot(fread(file.path(REF_DIR, "cruchaga_panel.csv")))
cr_sent  <- split_uniprot(fread(file.path(REF_DIR, "cruchaga_sentinels.csv")))
wo_panel <- split_uniprot(fread(file.path(REF_DIR, "western_olink_panel.csv")))
wo_sent  <- split_uniprot(fread(file.path(REF_DIR, "western_olink_sentinels.csv")))

# A complex is matched as one assay against UKB-PPP (same Olink platform), but the
# two UKB tables write it differently: the panel (ST3) as "P29459;P29460", sometimes
# repeating an accession ("Q96QH8;Q96QH8"); the sentinels (ST16 protein IDs) and our
# assay_annot as "P29459_P29460". norm_ids() gives all three one spelling: accessions
# deduplicated, sorted, joined by "_". Before 2026-10-28 the panel was read as written,
# so every complex of ours counted as "not assayed by the reference" (tier C).
# norm_ids() is in lib/common.R.
uk_panel <- fread(file.path(REF_DIR, "ukbppp_panel.csv"))[, uniprot_id := norm_ids(uniprot_id)]
uk_sent  <- fread(file.path(REF_DIR, "ukbppp_sentinels.csv"))[, uniprot_id := norm_ids(uniprot_id)]

# ref_name says which study a row is from, so each call can say who assayed or reported it.
ref <- rbind(
  cr_sent[, .(fluid = "CSF", ref_name = "SomaScan", uniprot_id, ref_rsid = rsid, ref_chr = chr,
              ref_pos = pos_hg19, ref_cis_trans = cis_trans)],
  wo_sent[, .(fluid = "CSF", ref_name = "Olink", uniprot_id, ref_rsid = rsid, ref_chr = chr,
              ref_pos = pos_hg19, ref_cis_trans = cis_trans)],
  uk_sent[, .(fluid = "Plasma", ref_name = "UKB-PPP", uniprot_id, ref_rsid = rsid, ref_chr = chr,
              ref_pos = pos_hg19, ref_cis_trans = cis_trans)])
ref[ref_rsid == "", ref_rsid := NA_character_]

panel <- rbind(cr_panel[, .(fluid = "CSF", ref_name = "SomaScan", uniprot_id)],
               wo_panel[, .(fluid = "CSF", ref_name = "Olink", uniprot_id)],
               uk_panel[, .(fluid = "Plasma", ref_name = "UKB-PPP", uniprot_id)])
panel <- unique(panel[!is.na(uniprot_id) & uniprot_id != ""])

ref_study <- c(CSF = "Cruchaga CSF SomaScan 7k + Western CSF Olink HT",
               Plasma = "UKB-PPP plasma Olink 3k")
sent[, ref_study := ref_study[fluid]]
# The accession(s) each credible set is looked up under. Plasma: the whole complex,
# as above. CSF: each subunit, because Cruchaga's aptamers are split the same way, so
# a complex of ours counts as assayed if SomaScan measured any of its subunits.
keys <- rbind(
  sent[fluid == "Plasma" & !is.na(uniprot_id), .(cs_key, fluid, match_id = norm_ids(uniprot_id))],
  sent[fluid == "CSF" & !is.na(uniprot_id),
       .(match_id = unique(strsplit(uniprot_id, "_")[[1]])), by = .(cs_key, fluid)])
assayed <- keys[paste(fluid, match_id) %in% panel[, paste(fluid, uniprot_id)], unique(cs_key)]
sent[, ref_assayed := cs_key %in% assayed]
sent[is.na(uniprot_id), ref_assayed := NA]
assayed_by <- merge(keys, unique(panel[, .(fluid, match_id = uniprot_id, ref_name)]),
                    by = c("fluid", "match_id"))[
  , .(ref_assayed_by = paste(sort(unique(ref_name)), collapse = ";")), by = cs_key]

# ---- candidate reference sentinels: same protein, same chromosome, +/-1 Mb ---
# Two subunits of a CSF complex can reach the same split Cruchaga row; count it once.
# Deduplicate on the reference row, not on its values: UKB-PPP's ST9 and ST16 can both
# list one variant, and n_ref_in_window has always counted those as two.
ref[, ref_row := .I]
cand <- merge(merge(keys[cs_key %in% assayed], sent[, .(cs_key, snp, chr, pos)], by = "cs_key"),
              ref, by.x = c("fluid", "match_id"), by.y = c("fluid", "uniprot_id"),
              allow.cartesian = TRUE)
# A reference sentinel without an hg19 position (9 Cruchaga rows whose liftover failed,
# none with an rsID) may lie in the window: keep it, and the missing r2 below marks the
# credible set unresolved instead of novel.
message("reference sentinels without an hg19 position: ", ref[is.na(ref_pos), .N])
cand <- unique(cand[ref_chr == chr & (is.na(ref_pos) | abs(ref_pos - pos) <= WINDOW)],
               by = c("cs_key", "ref_row"))
cand[, c("match_id", "ref_row") := NULL]
message("credible sets with >=1 reference sentinel for the same protein within 1 Mb: ",
        uniqueN(cand$cs_key))

# ---- r2 in 1000G EUR --------------------------------------------------------
# One plink pass over the union of our sentinels and the candidate reference sentinels,
# rather than a pairwise call each. Self-pairs (identical rsID) never appear in the
# plink output, so they are set to 1 afterwards.
cand_ld <- cand[!is.na(ref_rsid) & ref_rsid != snp]
snp_union <- unique(c(sent$snp, cand_ld$ref_rsid))
extract_f <- file.path(TMP_DIR, "07b_ld_snps.txt")
index_f   <- file.path(TMP_DIR, "07b_ld_index.txt")
ld_prefix <- file.path(TMP_DIR, "07b_ld")
writeLines(snp_union, extract_f)
writeLines(unique(sent$snp), index_f)
# The prefix is fixed, so clear the last run's output: a failed plink must not leave an
# old .ld file to be read as this run's.
unlink(paste0(ld_prefix, c(".ld", ".log")))

plink_cmd <- paste(
  "module load PLINK/1.9b_6.21-x86_64 >/dev/null 2>&1;",
  "plink --bfile", BFILE, "--extract", extract_f, "--ld-snp-list", index_f,
  "--r2 --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0.05",
  "--out", ld_prefix)
status <- system2("bash", c("-lc", shQuote(plink_cmd)), stdout = NULL, stderr = NULL)
if (status != 0 || !file.exists(paste0(ld_prefix, ".ld")))
  stop("plink --r2 failed (exit ", status, "); see ", ld_prefix, ".log")
ld <- fread(paste0(ld_prefix, ".ld"))[, .(SNP_A, SNP_B, R2)]
ld <- unique(rbind(ld, ld[, .(SNP_A = SNP_B, SNP_B = SNP_A, R2)]))
message("LD pairs from plink: ", nrow(ld))

cand <- merge(cand, ld, by.x = c("snp", "ref_rsid"), by.y = c("SNP_A", "SNP_B"),
              all.x = TRUE)
cand[ref_rsid == snp, R2 := 1]

# A pair plink did not report is r2 < 0.05 only when both SNPs are in the panel. If
# either is missing (no rsID, a non-rsID lead, not in 1000G EUR) its r2 is unknown;
# counting it as "no match" would manufacture novelty, so it is flagged instead.
panel_ids <- fread(paste0(BFILE, ".bim"), header = FALSE, select = 2L)[[1]]
cand[, ld_missing := is.na(R2) & !(snp %in% panel_ids & ref_rsid %in% panel_ids)]
cand[is.na(R2) & !ld_missing, R2 := 0]
message("candidate pairs without r2 (SNP not in the LD panel): ", cand[ld_missing == TRUE, .N],
        " of ", nrow(cand))

best <- cand[order(-R2, na.last = TRUE),
             .(n_ref_in_window       = .N,
               n_ref_window_no_ld    = sum(ld_missing),
               ref_best_rsid         = ref_rsid[1],
               ref_best_pos          = ref_pos[1],
               ref_best_cis_trans    = ref_cis_trans[1],
               ref_best_r2           = R2[1],
               ref_best_study        = ref_name[1]), by = cs_key]

# protein-level cis view: did the reference report any cis signal for this protein,
# LD aside? This is the robust reading for cis, where the claim is about the protein.
prot_cis <- unique(ref[ref_cis_trans %like% "cis", .(fluid, match_id = uniprot_id, ref_name)])
prot_cis <- merge(keys, prot_cis, by = c("fluid", "match_id"))[
  , .(ref_protein_cis = TRUE, ref_protein_cis_by = paste(sort(unique(ref_name)), collapse = ";")),
  by = cs_key]

out <- merge(sent, best, by = "cs_key", all.x = TRUE)
out <- merge(out, prot_cis, by = "cs_key", all.x = TRUE)
out <- merge(out, assayed_by, by = "cs_key", all.x = TRUE)
if (nrow(out) != nrow(sent)) stop("a merge multiplied credible sets: ", nrow(sent), " -> ", nrow(out))
out[is.na(ref_protein_cis), ref_protein_cis := FALSE]
out[is.na(n_ref_in_window), `:=`(n_ref_in_window = 0L, n_ref_window_no_ld = 0L)]

# ---- tiers ------------------------------------------------------------------
out[, novelty_tier := fifelse(
  is.na(uniprot_id), "unmappable: no UniProt accession",
  fifelse(ref_assayed == FALSE, "C: protein not assayed by the reference",
  fifelse(!is.na(ref_best_r2) & ref_best_r2 >= R2_MATCH, "known: matched reference sentinel",
  fifelse(n_ref_window_no_ld > 0, "unresolved: LD unavailable for a reference sentinel",
          "A: no matching reference sentinel"))))]
out[, novelty_protein_cis := fifelse(
  is.na(uniprot_id), NA_character_,
  fifelse(ref_assayed == FALSE, "C: protein not assayed by the reference",
  fifelse(ref_protein_cis, "known: reference reports a cis-pQTL for this protein",
          "A: no cis-pQTL reported for this protein")))]

out[, signal := fifelse(is_cis, "cis", "trans")]

# One vocabulary for both fluids and both signal classes. cis is read at the protein
# level and trans at the locus level - that is the resolution each carries its
# information at - but the labels are the same, so a reader never has to hold two
# tier systems in their head at once.
short_tier <- function(x)
  fifelse(is.na(x), "unmappable: no UniProt accession",
  fifelse(grepl("^C:", x), "C: reference did not assay this protein",
  fifelse(grepl("^A:", x), "A: not reported by the reference",
  fifelse(grepl("^known", x), "known: reported by the reference",
  fifelse(grepl("^unresolved", x), "unresolved: LD unavailable",
          "unmappable: no UniProt accession")))))

out[, novelty := short_tier(fifelse(signal == "cis", novelty_protein_cis, novelty_tier))]
out[, novelty_resolution := fifelse(signal == "cis", "protein", "locus")]

# ---- primary cis association ------------------------------------------------
# The cis call above is protein level, so a second signal for a protein with a reported
# cis-pQTL counts as reported. novelty_primary_cis asks instead whether each protein's
# PRIMARY cis association (lowest lead P among its distinct cis signals) was reported,
# with UKB-PPP's rule as written: a reference sentinel for the same protein within 1 Mb
# matches if it has r2 >= 0.8 (1000G EUR, GRCh37) with ANY variant at P < 5e-8 in our
# BF2 summary statistics within 1 Mb of our lead, or with a member of our credible set
# (the signal itself, whatever its P). Matching on our lead alone misses a
# reported signal whenever the two studies pick different leads in one LD block.
# Secondary cis signals keep only the protein-level call.
#   A          protein without a reported cis-pQTL (as novelty_protein_cis)
#   B          protein with a reported cis-pQTL, primary association not matched
#   known      matched; unresolved: not matched, and a sentinel in the window has no r2
#   C          protein not assayed by the reference
# 07c adds western_cs_min_p, which flags a B association that is genome-wide significant
# in Western's full summary statistics although Western did not list it.
P_SIG <- 5e-8
SUMS <- file.path(GWAS_BASE, "bf2/clean/protein")
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
kept <- fread(KEPT)
dist <- distinct_signals(s1[is_cs_lead == TRUE & is_cis == TRUE &
                              paste(cohort, protein_assay) %in% kept[, paste(cohort, protein_assay)],
                            .(cs_key, cohort, protein_assay, cs_id, fluid, is_cis, p)])
setorder(dist, fluid, protein_assay, p, cs_key)
prim <- dist[, .SD[1], by = .(fluid, protein_assay)]
if (any(prim$cohort != "bf2")) stop("primary cis signals outside BF2: ", paste(unique(prim$cohort), collapse = ", "))
out[, is_primary_cis := fifelse(is_cis, cs_key %in% prim$cs_key, NA)]

# Only primaries of proteins with a reported cis-pQTL need the test; A and C are decided.
todo <- out[is_primary_cis %in% TRUE & startsWith(novelty_protein_cis, "known") & cs_key %in% cand$cs_key]
bim <- fread(paste0(BFILE, ".bim"), header = FALSE, select = c(1L, 2L, 4L), col.names = c("chr", "rsid", "pos"))
sig_snps <- parallel::mclapply(split(todo, by = "cs_key"), function(d) {
  f <- file.path(SUMS, sprintf("bf2__%s.tsv.gz", d$protein_assay))
  if (!file.exists(f)) stop("missing summary statistics: ", f)
  g <- fread(f, select = c("CHR", "POS", "P"), showProgress = FALSE)[
    CHR == d$chr & abs(POS - d$pos) <= WINDOW & P < P_SIG]
  data.table(cs_key = d$cs_key, snp = unique(bim[chr == d$chr][g, on = .(pos = POS), nomatch = 0L, rsid]))
}, mc.cores = NC)
if (any(vapply(sig_snps, inherits, TRUE, "try-error"))) stop("reading summary statistics failed: ",
  paste(unlist(sig_snps[vapply(sig_snps, inherits, TRUE, "try-error")]), collapse = "; "))
sig_snps <- unique(rbind(rbindlist(sig_snps), s1[cs_key %in% todo$cs_key, .(cs_key, snp)]))
message("primary cis signals tested with the significant-variant rule: ", nrow(todo),
        "; variants (P < ", P_SIG, " or credible set) per signal, median ", sig_snps[, .N, by = cs_key][, median(N)])

pcand <- cand[cs_key %in% todo$cs_key, .(cs_key, ref_rsid)]
refs <- unique(pcand[ref_rsid %in% panel_ids, ref_rsid])
sig_f <- file.path(TMP_DIR, "07b_sig_snps.txt"); sig_i <- file.path(TMP_DIR, "07b_sig_index.txt")
sig_prefix <- file.path(TMP_DIR, "07b_sig_ld")
writeLines(unique(c(sig_snps$snp, refs)), sig_f); writeLines(refs, sig_i)
unlink(paste0(sig_prefix, c(".ld", ".log")))
status <- system2("bash", c("-lc", shQuote(paste(
  "module load PLINK/1.9b_6.21-x86_64 >/dev/null 2>&1;",
  "plink --bfile", BFILE, "--extract", sig_f, "--ld-snp-list", sig_i,
  "--r2 --ld-window-kb 1000 --ld-window 999999 --ld-window-r2", R2_MATCH, "--out", sig_prefix))),
  stdout = NULL, stderr = NULL)
if (status != 0 || !file.exists(paste0(sig_prefix, ".ld")))
  stop("plink --r2 (significant variants) failed (exit ", status, "); see ", sig_prefix, ".log")
ldp <- fread(paste0(sig_prefix, ".ld"))[, .(ref_rsid = SNP_A, snp = SNP_B)]
# identical rsID is a match even when the variant is not in the LD panel
self <- unique(na.omit(pcand$ref_rsid))
ldp <- unique(rbind(ldp, ldp[, .(ref_rsid = snp, snp = ref_rsid)], data.table(ref_rsid = self, snp = self)))
hit <- merge(merge(pcand, ldp, by = "ref_rsid", allow.cartesian = TRUE), sig_snps, by = c("cs_key", "snp"))
hit <- hit[order(ref_rsid), .(primary_ref_rsid = ref_rsid[1]), by = cs_key]
unres <- pcand[is.na(ref_rsid) | !ref_rsid %in% panel_ids, unique(cs_key)]

out <- merge(out, hit, by = "cs_key", all.x = TRUE)
out[, novelty_primary_cis := NA_character_]
out[is_primary_cis %in% TRUE, novelty_primary_cis := fcase(
  is.na(uniprot_id), "unmappable: no UniProt accession",
  !startsWith(novelty_protein_cis, "known"), novelty_protein_cis,   # A or C
  !is.na(primary_ref_rsid), "known: matched reference sentinel",
  cs_key %in% unres, "unresolved: LD unavailable for a reference sentinel",
  default = "B: protein reported, primary association not matched")]
setorder(out, fluid, -is_cis, cs_key)

# ---- write ------------------------------------------------------------------
cols <- c("cs_key", "fluid", "signal", "protein_assay", "gene", "ensembl_gene_id",
          "uniprot_id", "snp", "chr", "pos", "ref_study", "ref_assayed",
          "n_ref_in_window", "n_ref_window_no_ld", "ref_best_rsid", "ref_best_pos",
          "ref_best_cis_trans", "ref_best_r2", "novelty_tier", "ref_protein_cis",
          "novelty_protein_cis", "novelty", "novelty_resolution",
          # which study (SomaScan / Olink for CSF, UKB-PPP for plasma), ";"-joined
          "ref_assayed_by", "ref_protein_cis_by", "ref_best_study",
          "is_primary_cis", "novelty_primary_cis", "primary_ref_rsid")
res <- out[, ..cols]
fwrite(res, file.path(OUT_DIR, "07_novelty_credible_sets.csv.gz"))

# The manuscript-facing table: counts by fluid x cis/trans x tier, in the single
# vocabulary above. The locus-level read of cis is kept alongside because it is the
# more conservative number and a reviewer will ask for it.
# S1 (and res) carry every credible set of every assay, so the per-set labels stay
# there; the counts are distinct signals of the analysed assays, the unit every other
# signal count in the paper uses (README).
sig <- s1[is_cs_lead == TRUE &
            paste(cohort, protein_assay) %in% kept[, paste(cohort, protein_assay)],
          .(cs_key, cohort, protein_assay, cs_id)]
sig <- distinct_signals(sig)
res_sig <- res[cs_key %in% sig$cs_key]
if (nrow(res_sig) != nrow(sig))
  stop("distinct signals without a novelty call: ", nrow(sig) - nrow(res_sig))
summ <- rbind(
  res_sig[, .(resolution = "reported", N = .N), by = .(fluid, signal, tier = novelty)],
  res_sig[signal == "cis", .(resolution = "cis, locus level", N = .N),
          by = .(fluid, signal, tier = short_tier(novelty_tier))])
setorder(summ, fluid, signal, resolution, tier)

# ---- fold the call back into S1 ---------------------------------------------
# Tiers are final, so S1 can
# now carry the novelty call: credible-set level, joined on cs_key by add_novelty()
# (lib/common.R), which 06 also calls when it rebuilds S1.
# na.strings: 06 writes NA as an empty field and "" as a quoted one; reading the empty
# field as "" would make fwrite quote it, so S1 changed on every first 07b run (r023).
s1_new <- add_novelty(fread(S1, na.strings = c("", "NA")),
                      file.path(OUT_DIR, "07_novelty_credible_sets.csv.gz"))
if (is.null(s1_new)) stop("some S1 credible-set leads did not receive a novelty call")
fwrite(s1_new, S1)

for (f in "07_novelty_credible_sets.csv.gz") {
  p <- file.path(OUT_DIR, f)
  if (!file.exists(p)) stop("expected output missing: ", p)
}
print(dcast(summ[resolution == "reported"], fluid + signal + tier ~ .,
            value.var = "N", fun.aggregate = sum))
message("primary cis associations:")
print(dcast(res[is_primary_cis %in% TRUE], novelty_primary_cis ~ fluid, value.var = "cs_key", fun.aggregate = length))
message("07b done")
