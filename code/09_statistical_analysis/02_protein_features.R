#!/usr/bin/env Rscript
# Which intrinsic protein properties predict carrying a fine-mapped pQTL.
#
# Analysis and plotting in one file. They used to be two scripts; the
# split existed only by analogy with 01a/01c, where it earns its keep because the
# hotspot analysis is slow. Here both halves run in about three seconds together, so
# the split bought nothing and cost a file.
#
# Nine structural UniProt features (length, signal peptide, sites, localization) are
# tested against the fine-mapped signals, in each fluid and for cis and trans
# separately, under two models: logistic on carrying any credible set, and negative
# binomial on how many. Annotation-depth counts (GO, Pfam, InterPro, Reactome,
# DrugBank) were dropped on 2026-11-13: they measure how much a protein has been
# studied more than what it is, and were hard to read. Standard errors are clustered by
# UniProt identifier so repeated assays of one protein are not treated as independent.
#
# Inputs:
#   00_kept_assays.csv          deduped analyzed set (run 00 first)
#   info/bf_protein/uniprot/    UniProt feature annotation
# Outputs:
#   02_uniprot_feature_results.csv   one row per fluid x signal type x model x feature
#   02_uniprot_feature_data.csv.gz   the per-assay frame, continuous features UNSCALED
#                                    (the models use within-fluid z-scores of them)
# Usage:
#   Rscript code/09_statistical_analysis/02_protein_features.R

suppressPackageStartupMessages({
  library(data.table)
  library(MASS)
})

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1) args[1] else
  "output/02_uniprot_feature_results.csv"
data_csv <- if (length(args) >= 2) args[2] else
  "output/02_uniprot_feature_data.csv.gz"

annot_file <- "info/bf_protein/assay_annot.csv"
feature_file <- "info/bf_protein/uniprot/features.csv"
finemap_file <- "info/bf_protein/finemap_lead_snps.csv"

annot <- fread(annot_file)[cohort == "bf2"]
features <- fread(feature_file)
# The count model (negative binomial on how many credible sets) is only meaningful
# on distinct signals: overlapping windows resolve one signal repeatedly, and that
# redundancy is not a property of the protein. The logistic model is unaffected,
# because dedup never removes a protein's last credible set.
source("code/09_statistical_analysis/lib/cs_signals.R")
source("code/09_statistical_analysis/lib/common.R")   # assay_platform, assay_fluid, CHR_LEN_HG19
finemap <- fread(
  finemap_file,
  select = c("cohort", "protein_assay", "cs_id", "is_cis")
)[cohort == "bf2"]
finemap <- distinct_signals(finemap)

annot[, fluid := assay_fluid(tag)]
annot[, platform := assay_platform(tag)]
# SomaLogic is plasma-only and trans-inflated (00_dedup_assays.R drops it too); exclude it so
# the CSF-vs-plasma comparison rests on a fair Olink (and CSF NuLISA) basis.
annot <- annot[platform != "SomaLogic"]
annot <- annot[!is.na(fluid) & !is.na(uniprot_id) & nzchar(uniprot_id)]
# Restrict to the deduped analyzed set (N>1000, one highest-N assay per protein
# per fluid; code/00_dedup_assays.R) so each protein is counted once per fluid.
kept_assays <- fread("output/00_kept_assays.csv")
annot <- annot[protein_assay %in% kept_assays$protein_assay]

cs_counts <- finemap[, .(
  n_cis_cs = sum(is_cis == TRUE),
  n_trans_cs = sum(is_cis == FALSE)
), by = protein_assay]

dat <- merge(annot, cs_counts, by = "protein_assay", all.x = TRUE)
dat[is.na(n_cis_cs), n_cis_cs := 0L]
dat[is.na(n_trans_cs), n_trans_cs := 0L]
dat <- merge(dat, features, by = "uniprot_id", all.x = TRUE)
# The sample-size covariate is the GWAS N (00b), not N, which counts every measured
# sample whether genotyped or not.
n0 <- nrow(dat)
dat <- merge(dat, fread("output/00b_gwas_qc.csv",
                        select = c("protein_assay", "gwas_n")), by = "protein_assay")
if (nrow(dat) != n0 || anyNA(dat$gwas_n)) stop("00b_gwas_qc.csv does not cover the analysed assays; run 00b")

# Protein length is heavy-tailed, so it enters on a log2 scale.
dat[, `:=`(
  log2_length = log2(length),
  has_transmembrane = as.integer(n_transmembrane > 0),
  has_domain = as.integer(n_domain > 0),
  has_disulfide = as.integer(n_disulfide > 0),
  has_glycosylation = as.integer(n_glycosylation > 0),
  has_phosphorylation = as.integer(n_phosphorylation > 0),
  sample_n_z = as.numeric(scale(gwas_n))
)]

# Structural features: properties of the molecule itself.
feature_specs <- data.table(
  feature = c(
    "log2_length", "signal_peptide", "has_transmembrane", "has_domain",
    "has_disulfide", "has_glycosylation", "has_phosphorylation",
    "is_secreted", "is_membrane"
  ),
  label = c(
    "Protein length", "Signal peptide", "Transmembrane region", "Domain",
    "Disulfide bond", "Glycosylation site", "Phosphorylation site",
    "Secreted protein", "Membrane protein"
  ),
  feature_type = c("continuous", rep("binary", 8))
)

# Standardize the continuous protein features within fluid, so their coefficients
# are per standard deviation and comparable to each other. Binary coefficients
# compare presence vs absence and remain directly interpretable.
cont_features <- feature_specs[feature_type == "continuous", feature]
# The exported frame keeps the raw values: 03 pools the two fluids per gene, and
# z-scores taken within each fluid are on two different scales.
dat_raw <- copy(dat)
dat[, (cont_features) := lapply(.SD, function(x) as.numeric(scale(x))),
    .SDcols = cont_features, by = fluid]
dat[, platform := droplevels(factor(platform))]

if (!exists("PANEL_EXPORT_ONLY")) {
dir.create(dirname(data_csv), recursive = TRUE, showWarnings = FALSE)
# Columns follow feature_specs, so a feature added above lands here without a
# second edit.
fwrite(
  dat_raw[, c("fluid", "uniprot_id", "protein_assay", "gene", "platform", "N",
          "n_cis_cs", "n_trans_cs", feature_specs$feature), with = FALSE],
  data_csv
)
}

# Cluster-robust sandwich covariance for Bernoulli GLMs and NB2 models.
cluster_vcov <- function(fit, cluster, model_name) {
  X <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  mu <- fitted(fit)

  if (model_name == "Presence (logistic)") {
    denom <- rep(1, length(mu))
    weight <- mu * (1 - mu)
  } else {
    theta <- fit$theta
    denom <- 1 + mu / theta
    weight <- mu / denom
  }

  score <- X * as.numeric((y - mu) / denom)
  score_cluster <- rowsum(score, group = cluster, reorder = FALSE)
  bread <- tryCatch(
    solve(crossprod(X, X * as.numeric(weight))),
    error = function(e) qr.solve(crossprod(X, X * as.numeric(weight)))
  )
  meat <- crossprod(score_cluster)

  n <- nrow(X)
  p <- ncol(X)
  g <- nrow(score_cluster)
  correction <- if (g > 1L && n > p) (g / (g - 1)) * ((n - 1) / (n - p)) else 1
  correction * bread %*% meat %*% bread
}

extract_term <- function(
    fit, term, model_name, n_total, n_positive, n_clusters, cluster) {
  beta <- coef(fit)
  if (!term %in% names(beta)) return(NULL)
  vc <- tryCatch(cluster_vcov(fit, cluster, model_name), error = function(e) NULL)
  if (is.null(vc) || !term %in% rownames(vc)) return(NULL)
  variance <- unname(vc[term, term])
  if (!is.finite(variance) || variance <= 0) return(NULL)
  se <- sqrt(variance)
  z <- beta[term] / se
  data.table(
    model = model_name,
    estimate = unname(beta[term]),
    std_error = unname(se),
    p = 2 * pnorm(abs(z), lower.tail = FALSE),
    n_total = n_total,
    n_positive = n_positive,
    n_clusters = n_clusters
  )
}

fit_one <- function(d, count_col, feature, feature_type) {
  keep <- complete.cases(d[, c(count_col, feature, "sample_n_z", "platform"),
                           with = FALSE])
  x <- d[keep]
  if (nrow(x) < 100L || uniqueN(x[[feature]]) < 2L) return(NULL)

  x[, any_signal := as.integer(get(count_col) > 0L)]
  if (uniqueN(x$any_signal) < 2L) return(NULL)

  rhs_terms <- c(feature, "sample_n_z")
  if (uniqueN(x$platform) > 1L) rhs_terms <- c(rhs_terms, "platform")
  rhs <- paste(rhs_terms, collapse = " + ")
  logistic <- tryCatch(
    glm(as.formula(paste("any_signal ~", rhs)), data = x, family = binomial()),
    error = function(e) NULL
  )
  nb <- tryCatch(
    glm.nb(as.formula(paste(count_col, "~", rhs)), data = x),
    error = function(e) NULL
  )

  term <- feature
  n_positive <- sum(x[[count_col]] > 0L)
  n_clusters <- uniqueN(x$uniprot_id)
  out <- rbindlist(list(
    if (!is.null(logistic)) {
      extract_term(
        logistic, term, "Presence (logistic)", nrow(x), n_positive,
        n_clusters, x$uniprot_id
      )
    },
    if (!is.null(nb)) {
      extract_term(
        nb, term, "Count (negative binomial)", nrow(x), n_positive,
        n_clusters, x$uniprot_id
      )
    }
  ), fill = TRUE)
  if (!nrow(out)) return(NULL)
  out[, feature_type := feature_type]
  out
}

results <- list()
for (fluid_i in c("CSF", "Plasma")) {
  d <- dat[fluid == fluid_i]
  for (signal_i in c("cis", "trans")) {
    count_col <- paste0("n_", signal_i, "_cs")
    for (i in seq_len(nrow(feature_specs))) {
      r <- fit_one(
        d, count_col, feature_specs$feature[i], feature_specs$feature_type[i]
      )
      if (is.null(r)) next
      r[, `:=`(
        fluid = fluid_i,
        signal_type = signal_i,
        feature = feature_specs$feature[i],
        feature_label = feature_specs$label[i]
      )]
      results[[length(results) + 1L]] <- r
    }
  }
}
res <- rbindlist(results, fill = TRUE)
# fit_one / extract_term return NULL on a failed or degenerate fit. Every cell of the
# grid is expected, and a lost one would also shrink the Bonferroni denominator below.
n_expected <- 2L * 2L * nrow(feature_specs) * 2L   # fluid x signal type x feature x model
if (nrow(res) != n_expected) {
  got <- res[, paste(fluid, signal_type, feature, model)]
  all <- CJ(f = c("CSF", "Plasma"), s = c("cis", "trans"), x = feature_specs$feature,
            m = c("Presence (logistic)", "Count (negative binomial)"))[, paste(f, s, x, m)]
  stop(n_expected - nrow(res), " feature models not fitted, e.g. ", setdiff(all, got)[1])
}

res[, `:=`(
  conf_low = estimate - qnorm(0.975) * std_error,
  conf_high = estimate + qnorm(0.975) * std_error
)]
# Bonferroni correction across all tested feature associations: one global
# threshold so the CSV flag and the plot's dashed line use a single criterion.
bonf_alpha <- 0.05
res[, p_bonferroni := pmin(1, p * .N)]
res[, significant := p_bonferroni < bonf_alpha]

setcolorder(res, c(
  "fluid", "signal_type", "model", "feature", "feature_label",
  "feature_type", "estimate", "std_error", "conf_low", "conf_high",
  "p", "p_bonferroni", "significant", "n_total", "n_positive", "n_clusters"
))
setorder(res, fluid, signal_type, model, p_bonferroni, feature)

# PANEL_EXPORT_ONLY means another script is sourcing this one for the panel object;
# it should not have file side effects. See the same guard on the ggsave below.
if (!exists("PANEL_EXPORT_ONLY")) {
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(res, out_csv)
  message("Analysis assays: ", nrow(dat))
  message("Unique fluid-specific UniProt targets: ",
          uniqueN(dat[, paste(fluid, uniprot_id)]))
  message("Saved results: ", out_csv)
  message("Saved analysis data: ", data_csv)
}

# ---- Figure -------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
source("code/09_statistical_analysis/lib/fig_theme.R")   # theme_pgd, save_fig

# The analysis half above resolves its paths relative to the repo root; the figure
# goes next to the tables it is drawn from, whatever out_csv the caller passed.
OUT_DIR <- dirname(out_csv)

COLS <- FLUID_COL
FLUIDS <- c("CSF", "Plasma")



base_theme <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -2)
  )

# ---- Figure: which proteins carry a fine-mapped pQTL -------------------------
# The estimates from 02_protein_features.R, shown as -log10 p with the direction in the
# symbol. p rather than the effect size: the two models are on different scales (log
# odds and a negative-binomial rate ratio) and cannot share an axis.
P_CAP <- 20
# `res` is the table written above, used in memory rather than read back, so that
# sourcing this file for the panel alone (PANEL_EXPORT_ONLY) cannot pick up a stale
# CSV from an earlier run.
res[, neglog10p := pmin(-log10(p), P_CAP)]
res[, direction := factor(fifelse(estimate >= 0, "positive", "negative"),
                          levels = c("positive", "negative"))]
# Ordered by the single most significant point anywhere in the grid, so the strongest
# feature anchors the left edge.
res[, feature_label := factor(
  feature_label,
  levels = res[, .(min_p = min(p, na.rm = TRUE)), by = feature_label][
    order(min_p), feature_label])]
res[, signal_type := factor(signal_type, levels = c("cis", "trans"))]
res[, fluid := factor(fluid, levels = FLUIDS)]
res[, model := factor(model,
  levels = c("Count (negative binomial)", "Presence (logistic)"),
  labels = c("Number of credible sets", "Presence of any credible set"))]

# Colour is the fluid, as everywhere in the manuscript (CSF blue, plasma orange); cis
# full, trans a light tint of the same hue; triangle orientation, direction. One plot per
# fluid under a fluid title spanning its two model columns, so "CSF" and "Plasma" are
# each written once (2026-11-12; the former all-grey version repeated them per column).
tint <- function(col) grDevices::colorRampPalette(c(col, "white"))(5)[3]
feat_col <- c(CSF.cis = FLUID_COL[["CSF"]], CSF.trans = tint(FLUID_COL[["CSF"]]),
              Plasma.cis = FLUID_COL[["Plasma"]], Plasma.trans = tint(FLUID_COL[["Plasma"]]))
res[, col_key := paste(fluid, signal_type, sep = ".")]
feat_panel <- function(fl, first) {
  ggplot(res[fluid == fl], aes(feature_label, neglog10p, colour = col_key,
                               fill = col_key, shape = direction)) +
    geom_hline(yintercept = -log10(0.05 / nrow(res)), linetype = "dashed",
               colour = "grey45", linewidth = 0.35) +
    geom_point(size = 1.7, stroke = 0.25, position = position_dodge(width = 0.6)) +
    facet_grid(. ~ model, labeller = labeller(model = label_wrap_gen(16))) +
    scale_x_discrete(limits = rev) +
    coord_flip() +
    scale_shape_manual(values = c(positive = 24, negative = 25), name = "Direction") +
    scale_colour_manual(values = feat_col, guide = "none") +
    scale_fill_manual(values = feat_col, name = "pQTL type",
                      breaks = paste0(fl, c(".cis", ".trans")),
                      labels = c(expression(italic("cis")*"-pQTL"), expression(italic("trans")*"-pQTL"))) +
    # One x range in all four columns, so CSF and plasma read on the same scale.
    scale_y_continuous(limits = c(0, P_CAP), breaks = c(0, 5, 10, 15, P_CAP),
                       labels = c("0", "5", "10", "15", paste0("\u2265", P_CAP))) +
    labs(x = NULL, y = expression("\u2212" * log[10] ~ italic(P)), title = fl) +
    # One pQTL-type key for both fluids, in neutral greys (dark cis, light trans); the
    # hue is the fluid, named by the column titles.
    guides(fill = if (first) guide_legend(order = 1, override.aes = list(
             shape = 24, size = 1.8, fill = c("grey25", "grey75"), colour = c("grey25", "grey75")))
           else "none",
           shape = guide_legend(order = 2, override.aes = list(fill = "grey30", size = 1.8))) +
    base_theme +
    theme(panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          axis.text = element_text(size = 6)) +
    theme_pgd() +
    theme(plot.title = element_text(size = 7.2, hjust = 0.5),
          panel.border = element_blank(), axis.line = element_line(linewidth = .3, colour = "grey30"),
          axis.text.y = if (first) element_text(size = 6) else element_blank(),
          axis.ticks.y = if (first) element_line() else element_blank())
}
p_feat <- (feat_panel("CSF", TRUE) | feat_panel("Plasma", FALSE)) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# The same estimates turned on their side for Supplementary Fig. (09f_plot_novelty.R):
# features along x, the two models as rows and the fluids as columns, short enough for
# panel b to sit beside it. p_feat above stays for conferences/202609_adpd_2027.
# One title per panel, fluid and model together, instead of row and column strips.
res_h <- copy(res)
res_h[, panel := factor(paste(fluid, as.integer(model)),
  levels = c("CSF 1", "Plasma 1", "CSF 2", "Plasma 2"),
  labels = c("CSF, number of pQTL signals", "Plasma, number of pQTL signals",
             "CSF, presence of any pQTL signal", "Plasma, presence of any pQTL signal"))]
stopifnot(!anyNA(res_h$panel))
# Fuller labels than the short ones above, ordered by each feature's strongest P.
H_LAB <- c(log2_length = "Sequence length", has_domain = "Annotated domain")
res_h[, lab_h := fifelse(feature %in% names(H_LAB), H_LAB[feature], as.character(feature_label))]
h_ord <- res_h[, .(min_p = min(p, na.rm = TRUE)), by = lab_h][order(min_p), lab_h]
res_h[, lab_h := factor(lab_h, levels = h_ord)]
stopifnot(uniqueN(res_h$lab_h) == nrow(feature_specs))
p_feat_h <- ggplot(res_h, aes(lab_h, neglog10p, colour = col_key, fill = col_key,
                              shape = direction)) +
  geom_hline(yintercept = -log10(0.05 / nrow(res)), linetype = "dashed",
             colour = "grey45", linewidth = 0.35) +
  geom_point(size = 1.5, stroke = 0.25, position = position_dodge(width = 0.6)) +
  facet_wrap(~ panel, ncol = 2) +
  scale_shape_manual(values = c(positive = 24, negative = 25), name = "Direction") +
  scale_colour_manual(values = feat_col, guide = "none") +
  # Key names all four colours, fluid and pQTL type together.
  scale_fill_manual(values = feat_col, name = NULL, breaks = names(feat_col),
                    labels = c(expression("CSF" ~ italic("cis")*"-pQTL"), expression("CSF" ~ italic("trans")*"-pQTL"),
                               expression("Plasma" ~ italic("cis")*"-pQTL"), expression("Plasma" ~ italic("trans")*"-pQTL"))) +
  scale_y_continuous(limits = c(0, P_CAP), breaks = c(0, 5, 10, 15, P_CAP),
                     labels = c("0", "5", "10", "15", paste0("≥", P_CAP))) +
  labs(x = NULL, y = expression("−" * log[10] ~ italic(P))) +
  guides(fill = guide_legend(order = 1, nrow = 1, override.aes = list(
           shape = 24, size = 1.8, colour = unname(feat_col))),
         shape = guide_legend(order = 2, override.aes = list(fill = "grey30", size = 1.8))) +
  base_theme +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        strip.background = element_blank(),
        axis.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme_pgd() +
  theme(strip.text = element_text(size = 7), legend.position = "bottom")

# No figure is written here: 09f_plot_novelty.R draws p_feat_h as the supplementary
# figure, and conferences/202609_adpd_2027/code/prepare_figures.R takes p_feat (both source this
# file with PANEL_EXPORT_ONLY).
