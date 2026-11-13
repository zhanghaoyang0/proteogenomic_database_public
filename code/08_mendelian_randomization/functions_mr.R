# Per-pair MR helpers shared by 03_run_protein2disease.R and anything that needs the
# same instrument set afterwards (e.g. a Steiger check), so both see the SNPs the MR
# estimate was actually computed on.

suppressPackageStartupMessages(library(data.table))

# iv: one exposure's rows of mr_iv.tsv.gz (SNP, EA_e, NEA_e, BETA_e, SE_e, P_e).
# oc_full: the outcome's clean sumstats. Returns the IVs that survive harmonization,
# outcome beta and EAF aligned to EA_e, or NULL.
harmonize <- function(iv, oc_full) {
  oc <- oc_full[SNP %in% iv$SNP]
  if (!nrow(oc)) return(NULL)
  m <- merge(iv, oc[, .(SNP, EA_o = EA, NEA_o = NEA, EA_FREQ_o = EA_FREQ,
                        BETA_o = BETA, SE_o = SE, P_o = P, N_o = N)],
             by = "SNP")
  # A duplicated SNP on either side would count an instrument twice in IVW and Steiger.
  if (anyDuplicated(m$SNP)) stop("duplicated SNP after harmonizing: ", m$SNP[duplicated(m$SNP)][1])
  if (!nrow(m)) return(NULL)
  palind <- (m$EA_e == "A" & m$NEA_e == "T") | (m$EA_e == "T" & m$NEA_e == "A") |
            (m$EA_e == "C" & m$NEA_e == "G") | (m$EA_e == "G" & m$NEA_e == "C")
  maf <- pmin(m$EA_FREQ_o, 1 - m$EA_FREQ_o)
  m <- m[!(palind & maf > 0.42)]
  if (!nrow(m)) return(NULL)
  flip <- m$EA_o == m$NEA_e & m$NEA_o == m$EA_e
  same <- m$EA_o == m$EA_e & m$NEA_o == m$NEA_e
  m <- m[flip | same]
  if (!nrow(m)) return(NULL)
  flip2 <- m$EA_o == m$NEA_e & m$NEA_o == m$EA_e
  # Align the outcome to EA_e: the beta and the allele frequency (reported as EAF).
  m[flip2, `:=`(BETA_o = -BETA_o, EA_FREQ_o = 1 - EA_FREQ_o)]
  m
}

# Steiger directionality on a harmonized pair (Hemani 2017; the test TwoSampleMR's
# mr_steiger runs). Variance explained per SNP from z and N, r2 = z^2 / (z^2 + N - 2),
# summed over the IVs on each side, then the two correlations compared with Fisher's z.
# z rather than P because the strongest cis-pQTLs have P underflowing to 0. The outcome
# r2 is on the observed scale (total N, no case/control split), as get_r_from_pn gives.
# n_e: exposure GWAS N (one value); the outcome N is the mean of the IVs' N_o, as
# mr_steiger takes it.
steiger_test <- function(hm, n_e) {
  r2 <- function(b, se, n) { z2 <- (b / se)^2; z2 / (z2 + n - 2) }
  n_o  <- mean(hm$N_o)
  r2_e <- sum(r2(hm$BETA_e, hm$SE_e, n_e))
  r2_o <- sum(r2(hm$BETA_o, hm$SE_o, hm$N_o))
  # Summed r2 of 1 or more has no Fisher z (atanh(1) = Inf, above 1 NaN): no test, NA,
  # rather than a P computed from Inf or NaN. None in the 2026-10 run.
  if (!(r2_e < 1 && r2_o < 1))
    return(list(steiger_r2_exposure = r2_e, steiger_r2_outcome = r2_o,
                steiger_correct = NA, steiger_p = NA_real_))
  z <- (atanh(sqrt(r2_e)) - atanh(sqrt(r2_o))) / sqrt(1 / (n_e - 3) + 1 / (n_o - 3))
  list(steiger_r2_exposure = r2_e, steiger_r2_outcome = r2_o,
       steiger_correct = r2_e > r2_o, steiger_p = 2 * stats::pnorm(-abs(z)))
}
