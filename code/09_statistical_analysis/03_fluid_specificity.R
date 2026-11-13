#!/usr/bin/env Rscript
# Shared versus fluid-specific genetic regulation, at two levels:
#   1. protein level  - does a protein carry a fine-mapped pQTL in CSF, plasma, both, neither
#   2. credible-set level - does a credible set in one fluid share >=1 variant with any
#      credible set of the same protein in the other fluid
# Then test whether UniProt protein properties distinguish the protein-level classes.
#
# Restricted to the deduped analysed set (00_dedup_assays.R): one highest-N assay per
# protein per fluid, so each protein contributes once per fluid.
#
# Usage: Rscript code/09_statistical_analysis/03_fluid_specificity.R

suppressPackageStartupMessages(library(data.table))

OUT_DIR <- "output"
KEPT <- file.path(OUT_DIR, "00_kept_assays.csv")
FEATURES <- file.path(OUT_DIR, "02_uniprot_feature_data.csv.gz")
source("code/09_statistical_analysis/lib/common.R")   # GWAS_BASE, read_finemap_cs, check_finemap_cs
FM_DIR <- file.path(GWAS_BASE, "bf2/finemapping/protein")

for (f in c(KEPT, FEATURES)) if (!file.exists(f)) stop("missing input: ", f)
if (!dir.exists(FM_DIR)) stop("missing fine-mapping directory: ", FM_DIR)

source("code/09_statistical_analysis/lib/cs_signals.R")

kept <- fread(KEPT)

# ---- credible sets for every kept assay -------------------------------------
cs <- rbindlist(Map(read_finemap_cs, kept$cohort, kept$protein_assay,
                    MoreArgs = list(fm_dir = FM_DIR, select = c("cs_id", "snp", "is_cis"))),
                use.names = TRUE)
if (nrow(cs) == 0) stop("no credible sets read from ", FM_DIR)
check_finemap_cs(cs, kept, FM_DIR)
cs <- cs[, .(cohort, protein_assay, cs_id, snp, is_cis)]
# The credible-set-level share below is a count of signals, so the duplicates that
# overlapping +/-500 kb windows create have to go first (finemap/03_dedup_credible_sets.R).
# The protein-level classes are unaffected: dedup never removes a protein's last set.
# Units (a protein's pQTL, a credible set to classify) are the reported signals; the
# other-fluid lookup below uses every distinct credible set, as replication does, so a
# reported trans signal matched by a sub-study-wide set in the other fluid is shared.
cs_any <- add_cs_signal(cs)[is_representative == TRUE]
cs_any <- merge(cs_any, kept[, .(cohort, protein_assay, fluid, platform, ensembl_gene_id, gene)],
                by = c("cohort", "protein_assay"))
cs_any[, cs_key := paste(fluid, ensembl_gene_id, cs_id, sep = "|")]
cs <- distinct_signals(cs)
cs <- merge(cs, kept[, .(cohort, protein_assay, fluid, platform, ensembl_gene_id, gene)],
            by = c("cohort", "protein_assay"))
cs[, cs_key := paste(fluid, ensembl_gene_id, cs_id, sep = "|")]

n_missing <- kept[!paste(cohort, protein_assay) %in%
                    cs[, unique(paste(cohort, protein_assay))], .N]
message("assays with no credible-set file or no credible set: ", n_missing,
        " of ", nrow(kept))

# ---- protein level ----------------------------------------------------------
has <- dcast(kept[, .(ensembl_gene_id, gene, fluid,
                      hit = paste(cohort, protein_assay) %in%
                        cs[, unique(paste(cohort, protein_assay))])],
             ensembl_gene_id + gene ~ fluid, value.var = "hit",
             fun.aggregate = function(x) as.integer(any(x)), fill = NA_integer_)
both <- has[!is.na(CSF) & !is.na(Plasma)]
both[, class := fcase(
  CSF == 1L & Plasma == 1L, "shared",
  CSF == 1L & Plasma == 0L, "CSF-only",
  CSF == 0L & Plasma == 1L, "plasma-only",
  default = "neither"
)]
# A cis credible set in both fluids: the head-to-head denominator of the eQTL
# colocalization (04), which can only test cis loci. "shared" also counts proteins
# whose pQTL in one fluid is trans only, which can never colocalize there.
cis_assay <- cs[as.logical(is_cis) %in% TRUE, unique(paste(cohort, protein_assay))]
cis_has <- kept[paste(cohort, protein_assay) %in% cis_assay, .(n_fl = uniqueN(fluid)),
                by = ensembl_gene_id]
both[, cis_both := ensembl_gene_id %in% cis_has[n_fl == 2, ensembl_gene_id]]
stopifnot(both[cis_both == TRUE, all(class == "shared")])
fwrite(both, file.path(OUT_DIR, "03_fluid_specificity_protein.csv"))

# Analysis-set funnel, written so a set-definition panel can be drawn
# without hard-coding any of these counts. `has` is one row per protein assayed in at
# least one fluid, with NA where the protein was not assayed in that fluid, so the
# per-fluid assayed totals come straight off it.
sets <- data.table(
  step = c("assayed_csf", "assayed_plasma", "assayed_both",
           "pqtl_both", "pqtl_csf_only", "pqtl_plasma_only", "pqtl_neither", "cis_both"),
  n = c(has[!is.na(CSF), .N], has[!is.na(Plasma), .N], nrow(both),
        both[class == "shared", .N], both[class == "CSF-only", .N],
        both[class == "plasma-only", .N], both[class == "neither", .N],
        both[cis_both == TRUE, .N])
)
fwrite(sets, file.path(OUT_DIR, "03_analysis_sets.csv"))

cat("\n== protein level (assayed in both fluids: ", nrow(both), ") ==\n", sep = "")
print(both[, .N, by = class][order(-N)])

# ---- credible-set level -----------------------------------------------------
# Only proteins assayed in both fluids can be classified; a credible set is "shared"
# when it shares at least one variant with a credible set of the same protein in the
# other fluid.
cs2 <- cs[ensembl_gene_id %in% both$ensembl_gene_id]
cs2_any <- cs_any[ensembl_gene_id %in% both$ensembl_gene_id]
pairs <- merge(
  cs2_any[fluid == "CSF",    .(ensembl_gene_id, snp, csf_key = cs_key)],
  cs2_any[fluid == "Plasma", .(ensembl_gene_id, snp, pl_key = cs_key)],
  by = c("ensembl_gene_id", "snp"), allow.cartesian = TRUE
)
matched <- unique(c(pairs$csf_key, pairs$pl_key))

cs_lvl <- unique(cs2[, .(cs_key, fluid, ensembl_gene_id, gene, cs_id, is_cis)])
cs_lvl[, shared := cs_key %in% matched]
fwrite(cs_lvl, file.path(OUT_DIR, "03_fluid_specificity_credible_sets.csv"))

summ <- cs_lvl[, .(n = .N), by = .(fluid, is_cis, shared)]
summ[, pct_within_fluid_and_type := round(100 * n / sum(n), 1), by = .(fluid, is_cis)]
# Written once at the end, together with the feature test, as one summary file.

cat("\n== credible-set level (proteins assayed in both fluids) ==\n")
print(cs_lvl[, .(n_cs = .N, n_shared = sum(shared),
                 pct_shared = round(100 * mean(shared), 1)), by = .(fluid)][order(fluid)])
cat("\nby cis/trans:\n")
print(cs_lvl[, .(n_cs = .N, n_shared = sum(shared),
                 pct_shared = round(100 * mean(shared), 1)),
             by = .(fluid, is_cis)][order(fluid, -is_cis)])

# ---- protein properties vs protein-level class ------------------------------
fe <- fread(FEATURES)
# The same feature set 02 tested, taken from the columns it wrote rather than
# re-listed here, so the two analyses cannot drift apart. Continuous features are
# the log2_* ones (02 names them that way); the rest are binary.
FEAT <- setdiff(names(fe), c("fluid", "uniprot_id", "protein_assay", "gene",
                             "platform", "N", "n_cis_cs", "n_trans_cs"))
is_continuous <- function(f) startsWith(f, "log2_")
fe <- unique(fe[, c("gene", FEAT), with = FALSE], by = "gene")
m <- merge(both[class != "neither", .(gene, class)], fe, by = "gene")
m[, class := factor(class, levels = c("shared", "CSF-only", "plasma-only"))]

# Fixed seed: the simulated Fisher P values below otherwise change on every run, and
# Bonferroni calls near 0.05 with them.
set.seed(1)
res <- rbindlist(lapply(FEAT, function(f) {
  x <- m[[f]]
  if (uniqueN(x[!is.na(x)]) < 2) return(NULL)
  p <- if (is_continuous(f)) kruskal.test(x ~ m$class)$p.value
       else fisher.test(table(m$class, x), simulate.p.value = TRUE, B = 1e5)$p.value
  st <- m[, .(v = if (is_continuous(f)) median(get(f), na.rm = TRUE)
                  else mean(get(f), na.rm = TRUE)), by = class]
  data.table(feature = f, p = p,
             shared = st[class == "shared", v],
             csf_only = st[class == "CSF-only", v],
             plasma_only = st[class == "plasma-only", v])
}))
res[, p_bonferroni := pmin(1, p * .N)]
res[, significant := p_bonferroni < 0.05]
setorder(res, p)
# See the single summary file written at the end of this script.

cat("\n== protein properties vs class (n = ", nrow(m), ") ==\n", sep = "")
print(res, digits = 3)

# The credible-set counts and the feature test are both small derived tables from
# this one script, so they share a file with a `table` column rather than being two
# more names to keep track of.
out_summ <- file.path(OUT_DIR, "03_fluid_specificity_summary.csv")
fwrite(rbindlist(list(credible_set_counts = summ, protein_features = res),
                 idcol = "table", use.names = TRUE, fill = TRUE), out_summ)
message("\nWrote 03_fluid_specificity_{protein,credible_sets,summary}.csv to ", OUT_DIR)
