#!/usr/bin/env Rscript
# Previously unmatched primary CSF cis associations (plotting only).
#   Fig. 3 (09_novelty.png)
#     a  minor allele frequency x effect of every primary CSF cis association, coloured
#        by its 09a assoc_class in four display groups; inset, the counts per group (the
#        colour key)
#     b  LRRC25 locus
#     c  rs1985157-C across Alzheimer's disease risk, CSF and plasma LRRC25 and microglial
#        LRRC25 expression (09e_lrrc25_alleles.csv)
#   Supplementary Fig. 1 (09_supp_replication.png)
#     a  BF2 vs BF1 effect sizes of all credible-set sentinels (lib/repl_panel.R, 01a table)
#     b  BF2 vs Olink CSF study effects (07c), unmatched and BF1-untestable CSF signals
#     c  effects in amyloid-negative vs amyloid-positive participants (lib/repl_panel.R, 10b)
#   Supplementary protein-features figure (09_supp_protein_features.png): UniProt features
#     vs carrying a pQTL (p_feat_h, sourced from 02_protein_features.R without its file
#     side effects)
# Run 09e_lrrc25_locus.R once for the regional input tables.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
input <- out <- "output"
source("code/09_statistical_analysis/lib/fig_theme.R")
source("code/09_statistical_analysis/lib/cs_signals.R")
source("code/09_statistical_analysis/lib/repl_panel.R")   # make_repl, read_repl
FONT <- FIG_FAMILY
WINE <- NOVEL_COL[["novel"]]
MATCH <- NOVELTY_CLS_COL[["known"]]
# 09a classes in four display groups: the three that could not be called reported or
# unreported (detected in Western, unresolved, not measured) share one light grey. Matched
# is a light blue drawn underneath, so the red and light-red unmatched points stand out.
CLS <- data.table(
  assoc_class = c("known", "A", "B", "B, detected in Western", "unresolved", "C"),
  label = c("Matched", "No reported cis-pQTL", "Reported cis-pQTL, other signal",
            rep("Not classified", 3)))
DISP_COL <- c(Matched = NOVELTY_CLS_COL[["known"]], `No reported cis-pQTL` = WINE,
              `Reported cis-pQTL, other signal` = NOVELTY_CLS_COL[["B"]], `Not classified` = "grey88")
GROUP_COL <- c(Unmatched = WINE, Matched = MATCH)
theme_set(theme_compact())
update_geom_defaults("text", list(family = FONT))

# ---- primary CSF cis associations (09a), outside the MHC --------------------------
cand <- fread(file.path(input, "09a_novel_candidates.csv"))
pri <- cand[is_cis == TRUE & in_mhc == FALSE & assoc_class != "secondary"]
stopifnot(all(pri$assoc_class %in% CLS$assoc_class), !anyDuplicated(pri$protein_assay))
res <- fread(file.path(input, "09e_novel_vs_known.csv"))
pri[, grp := fcase(assoc_class %in% c("A", "B"), "Unmatched", assoc_class == "known", "Matched")]
stopifnot(pri[grp == "Unmatched", .N] == res[metric == "abs_beta", novel_n],
          pri[grp == "Matched", .N] == res[metric == "abs_beta", known_n])
s1 <- fread(file.path(input, "S1_pqtl_finemapping.csv.gz"),
            select = c("cs_key", "is_cs_lead", "bf1_beta"))[is_cs_lead == TRUE]
n0 <- nrow(pri)
pri <- merge(pri, s1[, .(cs_key, bf1_beta)], by = "cs_key", all.x = TRUE)
rp <- unique(fread("info/bf_protein/replication/bf1_targeted_replication.tsv.gz",
                   select = c("cohort", "bf2_assay", "SNP", "bf2_eaf", "bf2_beta"))[cohort == "bf2"],
             by = c("bf2_assay", "SNP"))
pri <- merge(pri, rp[, .(protein_assay = bf2_assay, snp = SNP, bf2_eaf, bf2_beta)],
             by = c("protein_assay", "snp"), all.x = TRUE)
stopifnot(nrow(pri) == n0, all(abs(pri$beta - pri$bf2_beta) < 1e-3, na.rm = TRUE))
pri[, maf := pmin(bf2_eaf, 1 - bf2_eaf)]
message(sprintf("a: %d of %d primary associations with a frequency", pri[!is.na(maf), .N], n0))

# ---- a: MAF x effect, inset counts -------------------------------------------------
pri <- merge(pri, CLS, by = "assoc_class", all.x = TRUE)
stopifnot(nrow(pri) == n0, !anyNA(pri$label))
cnt <- pri[, .N, by = label]
stopifnot(sum(cnt$N) == n0, setequal(cnt$label, names(DISP_COL)))
cnt[, row := factor(label, rev(names(DISP_COL)))]
fwrite(pri[, .N, by = .(assoc_class, label)], file.path(out, "09f_primary_classes.csv"))
inset <- ggplot(cnt, aes(N, row, fill = label)) +
  geom_col(width = .72) +
  geom_text(aes(label = N), hjust = -.2, size = pt_size(6.5), family = FONT) +
  scale_fill_manual(values = DISP_COL, guide = "none") +
  # cis in italics, as in the text.
  scale_y_discrete(labels = function(x) lapply(x, function(l) if (grepl("cis-pQTL", l))
    bquote(.(sub(" cis-pQTL.*", "", l)) ~ italic("cis") * "-pQTL" * .(sub(".*cis-pQTL", "", l))) else l)) +
  scale_x_continuous(limits = c(0, 820), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6.5),
        plot.background = element_rect(fill = "white", colour = NA))
sc <- pri[!is.na(maf)]
sc[, top := label %in% c("No reported cis-pQTL", "Reported cis-pQTL, other signal")]
setorder(sc, top)
f3a <- ggplot(sc, aes(maf, beta, fill = label)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = .25) +
  geom_point(aes(size = top, alpha = top), shape = 21, colour = "white", stroke = .2) +
  scale_fill_manual(values = DISP_COL, guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.3, `TRUE` = 1.9), guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = .8, `TRUE` = 1), guide = "none") +
  scale_x_continuous(limits = c(0, .5), breaks = seq(0, .5, .1), expand = c(.01, 0)) +
  labs(tag = "a", x = "Minor allele frequency", y = "Effect (s.d. per effect allele)") +
  inset_element(inset, left = .38, bottom = .74, right = 1, top = 1.02, align_to = "panel")

# ---- b: LRRC25 locus -----------------------------------------------------------------
tr <- fread(file.path(input, "09e_lrrc25_tracks.csv.gz"))
co <- fread(file.path(input, "09e_lrrc25_coloc.csv"))
lead_pos <- unique(tr[lead_variant == TRUE, POS])
stopifnot(length(lead_pos) == 1L, all(tr$P > 0),
          all(co[status == "ok" & primary == TRUE, PP4] > .99))
tr[, `:=`(mlp = -log10(P), ld = cut(r2, LD_BREAKS, LD_LABELS, right = FALSE))]
tr[, track := factor(track, c("AD GWAS", "CSF pQTL", "Plasma pQTL"))]
setorder(tr, r2, na.last = FALSE)
lr_csf <- co[grepl("o2csf", exposure_id) & primary == TRUE, PP4]
lr_plasma <- co[grepl("o7plasma", exposure_id) & primary == TRUE, PP4]
# PP4 above 0.99 shown as "> 0.99": more digits only separate values that mean the same
pp4_lab <- function(x) if (x > 0.99) "PP4 > 0.99" else sprintf("PP4 = %.2f", x)
track_text <- data.table(track = factor(c("AD GWAS", "CSF pQTL", "Plasma pQTL"), levels(tr$track)),
                        pp4_label = c("", pp4_lab(lr_csf), pp4_lab(lr_plasma)))
xr <- (lead_pos + c(-150000, 150000))/1e6
f3b <- ggplot(tr, aes(POS/1e6, mlp)) +
  geom_hline(yintercept = -log10(5e-8), colour = "grey68", linetype = "dotted", linewidth = .3) +
  geom_vline(xintercept = lead_pos/1e6, colour = LEAD_LINE_COL, linewidth = .35, alpha = .7) +
  geom_point(aes(colour = ld), size = 1, alpha = .8, stroke = 0) +
  geom_point(data = tr[lead_variant == TRUE], shape = 23, size = 2.1,
             fill = LEAD_COL, colour = "grey20", stroke = .3) +
  geom_text(data = tr[lead_variant == TRUE & track == "AD GWAS"],
            label = "rs1985157", nudge_x = .006, nudge_y = .8, hjust = 0,
            size = pt_size(6), family = FONT) +
  geom_text(data = track_text[pp4_label != ""],
            aes(x = xr[2] - .003, y = Inf, label = pp4_label), hjust = 1, vjust = 1.2,
            size = pt_size(6.5), family = FONT, inherit.aes = FALSE) +
  scale_colour_manual(values = setNames(LD_COLS, LD_LABELS), drop = FALSE,
                      na.value = "grey78", name = "LD r²", breaks = LD_LABELS) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.8, alpha = 1))) +
  # track names centred above each track (strips), not inside its top-left corner
  facet_wrap(~ track, ncol = 1, scales = "free_y", labeller = as_labeller(
    c(`AD GWAS` = "AD GWAS", `CSF pQTL` = "CSF LRRC25 pQTL", `Plasma pQTL` = "Plasma LRRC25 pQTL"))) +
  scale_x_continuous(limits = xr, breaks = scales::breaks_pretty(4), expand = c(0, 0)) +
  scale_y_continuous(breaks = function(l) pretty(c(0, l[2]*.75), 3),
                     expand = expansion(mult = c(.04, .35))) +
  labs(tag = "b", x = "Chromosome 19 (Mb, GRCh37)", y = "−log₁₀ P") +
  theme(strip.text = element_text(size = 6.5, hjust = 0.5, margin = margin(b = 1)),
        panel.spacing.y = unit(4, "pt"),
        legend.position = "bottom", legend.title = element_text(size = 6.5),
        legend.key.width = unit(5, "pt"), legend.spacing.x = unit(1, "pt"),
        legend.margin = margin(0, 0, 0, 0))

# c: rs1985157-C (the Alzheimer's disease risk allele) across the four signals of b and the
# microglial eQTL (09e_lrrc25_alleles.R), one row each with its own x axis and unit: the
# odds ratio of disease and the per-allele betas of protein level and expression are not
# on one scale, so only the direction against the null line is compared across rows.
al <- fread(file.path(input, "09e_lrrc25_alleles.csv"))
stopifnot(nrow(al) == 4, uniqueN(al$ea) == 1, all(al$is_lead))
al[, `:=`(est = beta, lo = beta - 1.96 * se, hi = beta + 1.96 * se, null = 0)]
al[signal == "AD GWAS", `:=`(est = exp(beta), lo = exp(lo), hi = exp(hi), null = 1)]
AL_LAB <- c(`AD GWAS` = "Alzheimer's disease risk (OR)", `CSF pQTL` = "CSF LRRC25 protein (\u03b2)",
            `Plasma pQTL` = "Plasma LRRC25 protein (\u03b2)", `Microglia eQTL` = "Microglial LRRC25 expression (\u03b2)")
AL_COL <- c(`AD GWAS` = "grey15", `CSF pQTL` = FLUID_COL[["CSF"]], `Plasma pQTL` = FLUID_COL[["Plasma"]],
            `Microglia eQTL` = "grey45")
al[, signal := factor(signal, levels = names(AL_LAB))]
f3c <- ggplot(al, aes(est, 0, colour = signal)) +
  geom_vline(aes(xintercept = null), colour = "grey60", linetype = "dashed", linewidth = .3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .45) +
  geom_point(size = 1.8) +
  facet_wrap(~ signal, ncol = 1, scales = "free_x", labeller = as_labeller(AL_LAB)) +
  scale_colour_manual(values = AL_COL, guide = "none") +
  scale_x_continuous(breaks = scales::breaks_pretty(3), expand = expansion(mult = .25)) +
  labs(tag = "c", title = sprintf("%s-%s", al$snp[1], al$ea[1]), x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        strip.text = element_text(size = 6.5, hjust = 0, margin = margin(b = 0)),
        plot.title = element_text(size = 7, hjust = 0), panel.spacing.y = unit(2, "pt"))
# a above c on the left, the locus (b) on the right.
fig3 <- (wrap_elements(full = (f3a / f3c) + plot_layout(heights = c(2.4, 1))) | wrap_elements(full = f3b)) +
  plot_layout(widths = c(100, 83))
save_fig(fig3, file.path(out, "figures"), "09_novelty", 125, dpi = 350, raster_device = "ragg")

# ---- supplementary: replication and protein features -------------------------------
wn <- fread(file.path(input, "07c_western_replication.csv"))[cs_key %in% pri[grp == "Unmatched", cs_key]]
stopifnot(nrow(wn) == pri[grp == "Unmatched", .N])
rep_bf <- res[metric == "bf1_replicated"]
rd <- data.table(study = c("BioFINDER-1", "BioFINDER-1", "External Olink"),
                 group = c("Unmatched", "Matched", "Unmatched"),
                 k = c(rep_bf$novel_k, rep_bf$known_k, sum(wn$western_status == "replicated")),
                 n = c(rep_bf$novel_n, rep_bf$known_n,
                       sum(wn$western_status %in% c("replicated", "tested, not replicated"))),
                 x = c(.87, 1.21, 2.0))
rd[, `:=`(pct = 100*k/n, label = paste0(k, "/", n))]
stopifnot(all(rd$k <= rd$n), all(rd$n > 0))
fwrite(rd, file.path(out, "09f_replication.csv"))
# b: BF2 vs Olink CSF study (Western) effects of the CSF signals whose replication there the
# text reports, the unmatched primary cis associations and the other signals BioFINDER-1
# could not test. rd (above) keeps the counts the text checks; the BF1 rates are in the text.
# The study's effects are on its own scale, so there is no identity line.
wb <- fread(file.path(input, "07c_western_replication.csv"))
wb[, grp := fcase(cs_key %in% pri[grp == "Unmatched", cs_key], "Unmatched",
                  startsWith(bf1_status, "not testable"), "Other")]
wb <- wb[!is.na(grp) & western_status %in% c("replicated", "tested, not replicated")]
stopifnot(wb[grp == "Unmatched", .N] == rd[study == "External Olink", n],
          wb[grp == "Unmatched", sum(western_status == "replicated")] == rd[study == "External Olink", k],
          !anyNA(wb[, .(beta_bf2, beta_western)]))
ws <- wb[, .(r = cor(beta_bf2, beta_western), pct = 100 * mean(western_status == "replicated")),
         keyby = grp]
# Key below the panel, as in a, so the long group names fit on one line with r and %.
W_COL <- c(Unmatched = WINE, Other = "grey60")
W_LAB <- c(Unmatched = "Unmatched to prior CSF reports", Other = "Other signals not testable in BioFINDER-1")
w_keys <- setNames(sprintf('"%s"~~italic(r) == "%.2f"*","~~"%.0f%%"', W_LAB[ws$grp], ws$r, ws$pct), ws$grp)
setorder(wb, -grp)   # unmatched drawn on top
s_a <- ggplot(wb, aes(beta_western, beta_bf2, colour = grp)) +
  geom_hline(yintercept = 0, linewidth = .2, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey70") +
  geom_point(size = 1.2, alpha = .75, shape = 16) +
  scale_colour_manual(values = W_COL, breaks = c("Unmatched", "Other"),
                      labels = function(x) parse(text = w_keys[x]), name = "Correlation, % replicated",
                      guide = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(tag = "b", title = "CSF", x = "Olink CSF study effect (replication)",
       y = "BioFINDER-2 effect (discovery)") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 7),
        axis.text = element_text(size = 6), plot.title = element_text(size = 7, hjust = .5),
        # Right-aligned: the key is about as wide as the whole column.
        legend.position = "bottom", legend.direction = "vertical", legend.justification = "right",
        legend.text = element_text(size = 6), legend.title = element_text(size = 6, hjust = 0),
        legend.key.size = unit(.6, "lines"), legend.key.spacing.y = unit(0, "pt"),
        aspect.ratio = 1) + legend_box

# BioFINDER-1 replication (a, two squares) and the Olink CSF study (b) in the first row,
# amyloid-negative vs amyloid-positive (c, two squares) in the second: equal squares on a
# 3-column grid so that they align, the last cell left empty.
rp <- make_repl(read_repl(input), key_below = TRUE)
am <- make_amyloid(input)
s_fig <- ((rp[[1]] + labs(tag = "a")) + rp[[2]] + s_a +
            (am[[1]] + labs(tag = "c")) + am[[2]] + plot_spacer()) +
  plot_layout(ncol = 3, widths = c(1, 1, 1)) &
  theme_pgd() & theme(plot.tag = element_text(size = 8, face = "bold"))
save_fig(s_fig, file.path(out, "figures"), "09_supp_replication", 160, dpi = 350, raster_device = "ragg",
         max_h_mm = FIG_H_MAX_SUPP_MM)

# Protein features: all proteins by whether they carry a pQTL (from 02, sourced in its own
# environment so its `res` and the like stay out of this script).
fe <- new.env()
assign("PANEL_EXPORT_ONLY", TRUE, envir = fe)
sys.source("code/09_statistical_analysis/02_protein_features.R", envir = fe)
save_fig(fe$p_feat_h, file.path(out, "figures"), "09_supp_protein_features", 95, dpi = 350, raster_device = "ragg",
         max_h_mm = FIG_H_MAX_SUPP_MM)
message("09f_plot_novelty.R done")
