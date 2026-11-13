# BioFINDER-2 vs BioFINDER-1 effect sizes of the credible-set sentinels (Supplementary
# Fig. 1a), drawn from 01a_fig1_replication.csv. In a lib because 09f_plot_novelty.R puts it
# above the replication of previously unmatched associations in the same figure.
# Needs data.table, ggplot2 and lib/fig_theme.R (FLUID_COL), loaded by the caller.

read_repl <- function(out_dir) {
  f <- file.path(out_dir, "01a_fig1_replication.csv")
  if (!file.exists(f)) stop("Missing ", f, ": run 01a_pqtl_data.R first.")
  d <- data.table::fread(f)
  d[, cis_lab := factor(cis_lab, levels = c("cis", "trans"))]
  d
}

# BF2 vs BF1 sentinel effect sizes, one square per fluid, each with its own key (cis and
# trans of that fluid). The two nested test sets (exact assay, cross-assay) are pooled and
# the distinction is left to the caption; 01a_pqtl_data.R defines them and the replication
# rule. Colour is fluid, shape is cis (triangle) or trans (circle), as in Fig. 1.
# Returns a two-plot patchwork, the y title on the left square only.
make_repl <- function(d, key_below = FALSE) {
  r <- copy(d)
  lim <- max(abs(c(r$bf2_beta, r$bf1_beta)), na.rm = TRUE)
  # Per fluid x cis/trans: correlation of the effects and the share replicated (same
  # direction and BF1 P < 0.05, the number the text quotes). These are the key labels.
  st <- r[, .(rho = cor(bf2_beta, bf1_beta, use = "complete.obs"),
              pct = 100 * mean(passed, na.rm = TRUE)), keyby = .(fluid_lab, cis_lab)]
  stopifnot(nrow(st) == 4, !anyNA(st$rho), !anyNA(st$pct))
  one <- function(fl) {
    rf <- r[fluid_lab == fl]
    sf <- st[fluid_lab == fl]
    keys <- setNames(sprintf('italic("%s")~~italic(r) == "%.2f"*","~~"%.0f%%"',
                             sf$cis_lab, sf$rho, sf$pct), as.character(sf$cis_lab))
    col <- FLUID_COL[[fl]]
    # BF1 on x: read left to right as "how did it hold up in BF1"; r is symmetric.
    ggplot(rf, aes(bf1_beta, bf2_beta)) +
      geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
      geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
      geom_point(aes(shape = cis_lab), colour = col, fill = col, size = 1.4, stroke = 0.35,
                 alpha = 0.7) +
      scale_shape_manual(values = c(cis = 24, trans = 21), labels = function(x) parse(text = keys[x]),
                         name = "Correlation, % replicated", drop = FALSE,
                         guide = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
      coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
      # The cohort names alone do not say which side is which; discovery/replication does.
      labs(title = fl, x = "BioFINDER-1 effect (replication)", y = "BioFINDER-2 effect (discovery)") +
      theme_bw(base_size = 7) +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 7),
            axis.text = element_text(size = 6), plot.title = element_text(size = 7, hjust = 0.5),
            # Key inside the panel, bottom right (the points lie along the identity line,
            # so that corner is empty), or below it with key_below.
            legend.position = if (key_below) "bottom" else "inside",
            legend.position.inside = c(0.98, 0.02), legend.direction = "vertical",
            legend.justification = if (key_below) "center" else c(1, 0),
            legend.text = element_text(size = 6), legend.title = element_text(size = 6, hjust = 0),
            legend.key.size = unit(0.6, "lines"), legend.key.spacing.y = unit(0, "pt")) +
      legend_box
  }
  # y title on the left square only; each square keeps its x title.
  one("CSF") | (one("Plasma") + labs(y = NULL))
}

# BF2 effects in amyloid-negative vs amyloid-positive participants (Supplementary Fig. 1c),
# from 10b_pqtl_context.csv, in the style of make_repl: one square per fluid, signed
# effects, identity line, shape cis/trans, key below with r per class. APOE +- 1 Mb is left
# out, as in 10b. The signals with an FDR-significant interaction are marked (named in the
# caption, not on the plot).
make_amyloid <- function(out_dir, fdr_max = 0.05) {
  f <- file.path(out_dir, "10b_pqtl_context.csv")
  if (!file.exists(f)) stop("Missing ", f, ": run code/09_statistical_analysis/local/10b on the Mac first.")
  r <- data.table::fread(f)[modifier == "abeta" & apoe_region == FALSE &
                              !is.na(beta_pos) & !is.na(beta_neg)]
  r[, cis_lab := factor(fifelse(is_cis, "cis", "trans"), levels = c("cis", "trans"))]
  lim <- max(abs(c(r$beta_neg, r$beta_pos)))
  st <- r[, .(rho = cor(beta_neg, beta_pos)), keyby = .(fluid, cis_lab)]
  stopifnot(nrow(st) == 4, !anyNA(st$rho), r[fdr_int < fdr_max, .N] >= 1)
  one <- function(fl) {
    rf <- r[fluid == fl]
    sf <- st[fluid == fl]
    keys <- setNames(sprintf('italic("%s")~~italic(r) == "%.2f"', sf$cis_lab, sf$rho),
                     as.character(sf$cis_lab))
    col <- FLUID_COL[[fl]]
    hit <- rf[fdr_int < fdr_max]
    ggplot(rf, aes(beta_neg, beta_pos)) +
      geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
      geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
      geom_point(aes(shape = cis_lab), colour = col, fill = col, size = 1.4, stroke = 0.35,
                 alpha = 0.7) +
      geom_point(data = hit, aes(shape = cis_lab), colour = "black", fill = "red3", size = 2,
                 stroke = 0.4, show.legend = FALSE) +
      scale_shape_manual(values = c(cis = 24, trans = 21), labels = function(x) parse(text = keys[x]),
                         name = "Correlation", drop = FALSE,
                         guide = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
      coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
      labs(title = fl, x = expression("A" * beta * "-negative effect"),
           y = expression("A" * beta * "-positive effect")) +
      theme_bw(base_size = 7) +
      theme(panel.grid = element_blank(), axis.title = element_text(size = 7),
            axis.text = element_text(size = 6), plot.title = element_text(size = 7, hjust = 0.5),
            legend.position = "bottom", legend.direction = "vertical",
            legend.text = element_text(size = 6), legend.title = element_text(size = 6, hjust = 0),
            legend.key.size = unit(0.6, "lines"), legend.key.spacing.y = unit(0, "pt")) +
      legend_box
  }
  one("CSF") | (one("Plasma") + labs(y = NULL))
}
