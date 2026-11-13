#!/usr/bin/env Rscript
# rs1985157 at the LRRC25 Alzheimer's disease locus, aligned to one allele across the four
# signals the Results compare: Alzheimer's disease risk (PGC-ALZ3, the 09e AD track), CSF and
# plasma LRRC25 pQTLs (09e tracks) and the microglial LRRC25 eQTL (SingleBrain, MG).
# Reads 09e_lrrc25_tracks.csv.gz (run 09e_lrrc25_locus.R first) and the cleaned SingleBrain
# summary statistics; summary-level only, runs anywhere in seconds.
# Writes output/09e_lrrc25_alleles.csv: signal, ea, nea, beta, se, p, with every beta for
# the same effect allele (the AD risk allele), and whether rs1985157 is that signal's lead.
# Usage (repo root): Rscript code/09_statistical_analysis/09e_lrrc25_alleles.R
suppressPackageStartupMessages(library(data.table))
OUT_DIR <- "output"
SNP_ID <- "rs1985157"   # not SNP: that is a column name, and DT[SNP == SNP] matches every row
SB_MG <- path.expand("~/db/shared_data/gwas/singlebrain/clean/MG/ENSG00000175489.tsv.gz")   # LRRC25
for (f in c(file.path(OUT_DIR, "09e_lrrc25_tracks.csv.gz"), SB_MG)) if (!file.exists(f)) stop("missing input: ", f)

tr <- fread(file.path(OUT_DIR, "09e_lrrc25_tracks.csv.gz"))
need <- c("AD GWAS", "CSF pQTL", "Plasma pQTL")
stopifnot(all(need %in% tr$track))
tr_lead <- tr[, .SD[which.min(P)], by = track]
lead_of <- setNames(tr_lead$SNP, tr_lead$track)
x <- tr[SNP == SNP_ID & track %in% need, .(signal = track, ea = EA, nea = NEA, beta = BETA, se = SE, p = P)]
mg <- fread(SB_MG)
mg_lead <- mg$SNP[which.min(mg$P)]
x <- rbind(x, mg[SNP == SNP_ID, .(signal = "Microglia eQTL", ea = EA, nea = NEA, beta = BETA, se = SE, p = P)])
stopifnot(nrow(x) == 4, !anyNA(x$beta))
# Align every signal to the Alzheimer's disease risk allele.
ad <- x[signal == "AD GWAS"]
risk <- if (ad$beta > 0) ad$ea else ad$nea
x[ea != risk, `:=`(beta = -beta, nea = ea, ea = risk)]
stopifnot(all(x$ea == risk), x[signal == "AD GWAS", beta > 0])
leads <- c(lead_of, `Microglia eQTL` = mg_lead)
x[, is_lead := unname(leads[signal]) == SNP_ID]
x[, snp := SNP_ID]
fwrite(x[, .(snp, signal, ea, nea, beta, se, p, is_lead)], file.path(OUT_DIR, "09e_lrrc25_alleles.csv"))
print(x)
message("09e_lrrc25_alleles.R done")
