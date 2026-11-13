#!/usr/bin/env Rscript
# Section 12 tests (docs/pipeline.md), from the per-protein table that
# code/09_statistical_analysis/local/12_fluid_clinical.R writes. Summary-level only; runs anywhere.
#   G  genetic specificity: |z_het| of the protein's lead 10a signal (smallest discovery P)
#   P  clinical effect heterogeneity (text; "phenotypic specificity" before 2026-11-13): mean z_pheno^2 over the outcomes, z_pheno = (beta_CSF -
#      beta_plasma) / se from the joint model; also for baseline and longitudinal outcomes
# Primary: Spearman(G, P), with the usual asymptotic P (2026-11-13; a 10,000-fold
# permutation of G over proteins gave the same P, 0.451 vs 0.457). Secondary: effect-specific
# against shared (Wilcoxon on P); Spearman(G, P) partial on the CSF-plasma residual
# correlation. Reliability control: CSF-minus-plasma asymmetry in mean |z| (separate
# models), shared proteins as the reference. Examples chosen by the plan's rule.
# Writes output/12b_protein_summary.csv (one row per protein), output/12b_tests.csv.
# Usage (repo root): Rscript code/09_statistical_analysis/12b_fluid_clinical_summary.R [in_csv]

suppressPackageStartupMessages(library(data.table))
OUT_DIR <- "output"
args <- commandArgs(trailingOnly = TRUE)
IN <- if (length(args)) args[1] else file.path(OUT_DIR, "12_fluid_clinical.csv")

r <- fread(IN)
fh <- fread(file.path(OUT_DIR, "10a_fluid_heterogeneity.csv"))[apoe_region == FALSE & !is.na(p_het)]
p_lead <- fread(file.path(OUT_DIR, "S1_pqtl_finemapping.csv.gz"))[is_cs_lead == TRUE, .(cs_key, p_disc = p)]
n0 <- nrow(fh)
fh <- merge(fh, unique(p_lead, by = "cs_key"), by = "cs_key")
if (nrow(fh) != n0) stop("10a signals without a discovery P in S1")

# ---- G: lead signal per protein ------------------------------------------------------------
lead <- fh[order(ensembl_gene_id, p_disc)][, .SD[1], by = ensembl_gene_id]
lead[, G := abs(delta / se_delta)]
lead <- lead[, .(ensembl_gene_id, G, class, resid_r)]

# ---- P: clinical effect heterogeneity ----------------------------------------------------------------
dif <- r[model == "joint" & term == "difference"][, z := est / se]
P <- dif[, .(P = mean(z^2), P_baseline = mean(z[analysis == "baseline"]^2),
             P_long = mean(z[analysis != "baseline"]^2), n_outcomes = .N), by = .(ensembl_gene_id, gene)]
# Reliability: per protein, mean |z| of the separate CSF and plasma models over outcomes
sep <- r[model %in% c("CSF only", "plasma only")][, z := est / se]
asym <- dcast(sep[, .(mz = mean(abs(z))), by = .(ensembl_gene_id, model)], ensembl_gene_id ~ model, value.var = "mz")
setnames(asym, c("CSF only", "plasma only"), c("mean_absz_csf", "mean_absz_plasma"))
asym[, asym := mean_absz_csf - mean_absz_plasma]
# Overall direction of a protein's difference: share of outcomes with |z_CSF| > |z_plasma|
dirn <- dcast(sep, ensembl_gene_id + analysis + outcome ~ model, value.var = "z")
setnames(dirn, c("CSF only", "plasma only"), c("zc", "zp"))
dirn <- dirn[, .(csf_stronger = mean(abs(zc) > abs(zp), na.rm = TRUE)), by = ensembl_gene_id]

pr <- Reduce(function(a, b) merge(a, b, by = "ensembl_gene_id"), list(P, lead, asym, dirn))
if (nrow(pr) != nrow(P)) stop("proteins lost joining G, asymmetry or direction")
message(sprintf("proteins: %d (outcomes per protein: median %d)", nrow(pr), as.integer(median(pr$n_outcomes))))

# ---- tests ----------------------------------------------------------------------------------
sp <- function(x, y) cor(x, y, method = "spearman", use = "complete.obs")
sp_p <- function(x, y) cor.test(x, y, method = "spearman", exact = FALSE)$p.value
partial_sp <- function(x, y, z) { rx <- resid(lm(rank(x) ~ rank(z))); ry <- resid(lm(rank(y) ~ rank(z))); cor(rx, ry) }
tests <- rbindlist(lapply(c("P", "P_baseline", "P_long"), function(v) data.table(
  test = "Spearman(G, P)", P_def = v, estimate = sp(pr$G, pr[[v]]), p = sp_p(pr$G, pr[[v]]), n = nrow(pr))))
es <- pr[class == "effect-specific"]; sh <- pr[class == "shared"]
tests <- rbind(tests,
  data.table(test = "P: effect-specific vs shared (Wilcoxon; estimate = median ratio)", P_def = "P",
             estimate = median(es$P) / median(sh$P), p = wilcox.test(es$P, sh$P)$p.value, n = nrow(es) + nrow(sh)),
  data.table(test = "Spearman(G, P) partial on resid_r", P_def = "P",
             estimate = partial_sp(pr$G, pr$P, pr$resid_r), p = NA_real_, n = nrow(pr)),
  data.table(test = "CSF-minus-plasma |z| asymmetry: shared (reference)", P_def = "asym",
             estimate = mean(sh$asym), p = t.test(sh$asym)$p.value, n = nrow(sh)),
  data.table(test = "CSF-minus-plasma |z| asymmetry: effect-specific", P_def = "asym",
             estimate = mean(es$asym), p = t.test(es$asym)$p.value, n = nrow(es)),
  data.table(test = "asymmetry: effect-specific vs shared (Wilcoxon)", P_def = "asym",
             estimate = mean(es$asym) - mean(sh$asym), p = wilcox.test(es$asym, sh$asym)$p.value, n = nrow(es) + nrow(sh)))
print(tests[, .(test, P_def, estimate = signif(estimate, 3), p = signif(p, 3), n)])

# ---- examples, by the plan's rule -------------------------------------------------------------
# direction of the overall result: CSF stronger if the effect-specific asymmetry exceeds the
# shared reference, plasma stronger otherwise
csf_dir <- mean(es$asym) > mean(sh$asym)
ex1 <- es[(if (csf_dir) csf_stronger > 0.5 else csf_stronger < 0.5)][order(-P)][1]
ex2 <- sh[order(-(mean_absz_csf + mean_absz_plasma))][1]
pr[, example := fifelse(ensembl_gene_id == ex1$ensembl_gene_id, "effect-specific example",
                 fifelse(ensembl_gene_id == ex2$ensembl_gene_id, "shared contrast", ""))]
message(sprintf("examples: %s (effect-specific, P = %.1f, %s stronger), %s (shared contrast)",
                ex1$gene, ex1$P, if (csf_dir) "CSF" else "plasma", ex2$gene))

# ---- single-fluid associations at BH FDR < 0.05 per model x analysis x outcome ----------
# The counts Results 12 and 12c's panel give; one BH here, read by both (code review
# 2026-10-30, R18: the qmd and 12c each repeated it).
sf <- r[model %in% c("CSF only", "plasma only")]
sf[, fdr := p.adjust(p, "BH"), by = .(model, analysis, outcome)]
sf_n <- sf[, .(n_tested = .N, n_fdr = sum(fdr < 0.05)), by = .(model, analysis, outcome)]

f1 <- file.path(OUT_DIR, "12b_protein_summary.csv"); f2 <- file.path(OUT_DIR, "12b_tests.csv")
f3 <- file.path(OUT_DIR, "12b_single_fluid_counts.csv")
if (length(args)) { f1 <- file.path("tmp/proteogenomic_database", basename(f1)); f2 <- file.path("tmp/proteogenomic_database", basename(f2))
                    f3 <- file.path("tmp/proteogenomic_database", basename(f3)) }   # a pilot input writes to tmp/
fwrite(pr[order(-P)], f1); fwrite(tests, f2); fwrite(sf_n, f3)
for (f in c(f1, f2, f3)) if (!file.exists(f)) stop("not written: ", f)
message("12b done -> ", f1)
