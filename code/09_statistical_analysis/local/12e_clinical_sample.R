#!/usr/bin/env Rscript
# LOCAL (Mac). Sample description for the section 12 clinical models (12_fluid_clinical.R):
# per outcome, the CU participants with CSF and plasma levels of at least one study protein,
# their observations, follow-up and events. The per-protein model N is in
# 12_fluid_clinical.csv (column n); this table is the union over proteins.
#
# Writes aggregate counts only (no participant rows):
#   output/12e_clinical_sample.csv   one row per outcome
#
# Usage (repo root, Mac):
#   Rscript code/09_statistical_analysis/local/12e_clinical_sample.R

suppressPackageStartupMessages(library(data.table))
source("code/09_statistical_analysis/lib/local_data.R")
source("code/09_statistical_analysis/lib/local_outcomes.R")   # load_mri_cu, load_clinical

OUT <- "output/12e_clinical_sample.csv"

# ---- the 12_fluid_clinical protein set and participants with both levels --------------------
fh <- fread("output/10a_fluid_heterogeneity.csv")[apoe_region == FALSE & !is.na(p_het)]
kept <- fread("output/00_kept_assays.csv")[, .(protein_assay, tag)]
prot <- unique(fh[, .(gene, csf_assay, plasma_assay)])
n0 <- nrow(prot)
prot <- merge(prot, kept[, .(csf_assay = protein_assay, csf_tag = tag)], by = "csf_assay")
prot <- merge(prot, kept[, .(plasma_assay = protein_assay, plasma_tag = tag)], by = "plasma_assay")
if (nrow(prot) != n0) stop("assays missing from 00_kept_assays.csv")
n_res <- fread("output/12_fluid_clinical.csv")
if (uniqueN(n_res$gene) != uniqueN(prot$gene))
  stop("protein set differs from 12_fluid_clinical.csv: ", uniqueN(prot$gene), " vs ", uniqueN(n_res$gene))

cov <- load_covar()
pv <- load_proteins(unique(rbind(prot[, .(protein_assay = csf_assay, panel_tag = csf_tag)],
                                 prot[, .(protein_assay = plasma_assay, panel_tag = plasma_tag)])), cov)
both <- Reduce(`|`, lapply(seq_len(nrow(prot)), function(i)
  !is.na(pv[[prot$csf_assay[i]]]) & !is.na(pv[[prot$plasma_assay[i]]])))
paired <- pv$sid[both]
message(sprintf("participants with both levels of at least one of %d proteins: %d", nrow(prot), length(paired)))

# ---- outcomes, as in 12_fluid_clinical.R -----------------------------------------------------
M <- load_mri_cu()
K <- load_clinical()
cu_base <- K$base[cu == TRUE, sid]
q <- function(x) sprintf("%.1f (%.1f-%.1f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
row <- function(outcome, d, time_col = "years") {
  fu <- d[, .(n_obs = .N, fu = max(get(time_col))), by = sid]
  data.table(outcome = outcome, participants = nrow(fu), observations = nrow(d),
             with_followup = fu[n_obs >= 2, .N],
             obs_per_participant = q(fu$n_obs), followup_years = q(fu[n_obs >= 2, fu]))
}
# Complete cases on the model covariates, as lm / lmer / coxph drop them in 12: age, sex
# and APOE e4 for every model, education for mPACC and AD-CI. (mean_protein_<panel> is
# missing only with the protein itself, which `paired` already requires.)
adj <- merge(cov[, .(sid, age, sex)], K$base[, .(sid, e4_carrier, edu)], by = "sid", all.x = TRUE)
adj_ok <- adj[!is.na(age) & !is.na(sex) & !is.na(e4_carrier), sid]
edu_ok <- adj[!is.na(edu), sid]
mri <- M$mri[sid %in% paired & sid %in% adj_ok]
cog <- K$long[sid %in% paired & sid %in% cu_base & !is.na(mpacc) & sid %in% adj_ok & sid %in% edu_ok]
sv  <- K$cu[sid %in% paired & sid %in% adj_ok & sid %in% edu_ok & !is.na(ad_ci)]
tab <- rbindlist(list(
  row("MRI (HV, ILV, thickness composites)", mri),
  row("mPACC", cog),
  data.table(outcome = "AD-CI (Cox)", participants = nrow(sv), observations = NA_integer_,
             with_followup = nrow(sv), obs_per_participant = NA_character_,
             followup_years = q(sv$t_years), events = sum(sv$ad_ci))
), fill = TRUE)
rng <- n_res[model == "joint", .(n_range = sprintf("%d-%d", min(n), max(n))),
             by = .(outcome_group = fifelse(outcome == "mPACC", "mPACC",
                                    fifelse(outcome == "AD-CI", "AD-CI (Cox)", "MRI (HV, ILV, thickness composites)")))]
n0 <- nrow(tab); tab <- merge(tab, rng, by.x = "outcome", by.y = "outcome_group", all.x = TRUE)
stopifnot(nrow(tab) == n0, !anyNA(tab$n_range))
# The table is the union over proteins, so no single model may have more participants.
n_max <- as.integer(sub(".*-", "", tab$n_range))
if (any(tab$participants < n_max)) stop("a model has more participants than the sample table")
print(tab)
fwrite(tab, OUT)
stopifnot(file.exists(OUT))
message("12e done: ", OUT)
