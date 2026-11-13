#!/usr/bin/env Rscript
# Regulatory-context panel definitions for Fig. 4 (03c_plot_annotation.R), from exported summaries.
#   a  brain context of colocalization for the head-to-head proteins: one column per
#      protein, one row per MetaBrain region / SingleBrain cell type (plus GTEx
#      peripheral as reference), tile = colocalizes in both fluids / CSF only / plasma
#      only; margin = share of each state per row
#   b  where GTEx support lies by protein-level pQTL detection class (04)
#   c-e  GRN (shared), ROR1 (CSF only) and THY1 (plasma only) loci from 10e
# Sourced by 03c_plot_annotation.R, which draws the figure; this file writes nothing (its
# standalone composition, 10_regulatory_context.png, was dropped on 2026-11-13).
# The GTEx panel is deliberately protein-level: current coloc exports cannot identify
# tested negatives per SuSiE signal (see docs/issues.md). Never join coloc cs_id
# to S1 cs_id: the former is a window index, the latter a credible-set identifier.
# For the same reason a blank tile in a means no retained colocalization, not a
# tested negative.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
# Reads exported summaries from output/; it needs nothing from 10_plot.R.
IN_DIR <- OUT_DIR <- "output"
source("code/09_statistical_analysis/lib/fig_theme.R")   # theme_pgd, save_fig, pt_size, axes_fig, LD colours
source("code/09_statistical_analysis/lib/common.R")      # COLOC_PP4
source("code/09_statistical_analysis/lib/eqtl_coloc_scope.R")   # CELLTYPE_LAB

# Axis text is set here; titles, legends and tags come from theme_pgd().
locus_axes <- axes_fig + theme(plot.title = element_text(size = 7))   # axes_fig plus the panel-title size, not a copy

# ---- b: tissue class of GTEx support by protein-level class ------------------------------
# Among proteins colocalizing anywhere in GTEx: brain only, brain and peripheral, or
# peripheral only (brain-only = n - n_peripheral; both = n_brain + n_peripheral - n).
gradient <- fread(file.path(IN_DIR, "04_eqtl_support_by_specificity.csv"))[
  denominator == "colocalizing anywhere in GTEx"]
stopifnot(nrow(gradient) == 3L, all(gradient$n_proteins > 0))
gb <- gradient[, .(`Brain only` = n_proteins - n_peripheral,
                   `Brain + peripheral` = n_brain + n_peripheral - n_proteins,
                   `Peripheral only` = n_proteins - n_brain), by = .(class, n_proteins, fisher_p)]
gb <- melt(gb, id.vars = c("class", "n_proteins", "fisher_p"), variable.name = "where", value.name = "n")
# The three parts sum to n_proteins by construction; a negative part means the
# brain/peripheral counts disagree with n_proteins.
stopifnot(all(gb$n >= 0))
gb[, `:=`(class = factor(class, levels = c("plasma-only", "shared", "CSF-only"),
                         labels = c("Plasma only", "Both fluids", "CSF only")),
          where = factor(where, levels = rev(c("Brain only", "Brain + peripheral", "Peripheral only"))))]
# Upright bars, from plasma-only to CSF-only, so the panel stays narrow beside b.
# Brain is the beta red (BRAIN_COL); blue and orange are taken by CSF and plasma.
WHERE_COL <- c(`Brain only` = BRAIN_COL[["specific"]], `Brain + peripheral` = BRAIN_COL[["any"]],
               `Peripheral only` = "grey80")
n_by_class <- unique(gb[, .(class, n_proteins)])
x_labs <- setNames(paste0(sub(" .*", "", as.character(n_by_class$class)), "\n(", n_by_class$n_proteins, ")"),
                   as.character(n_by_class$class))
bridge <- ggplot(gb, aes(class, n / n_proteins, fill = where)) +
  geom_col(width = 0.62) +
  # Legend below the bars: the two "only" classes on one row, the mixed class below.
  scale_fill_manual(values = WHERE_COL, name = NULL,
                    breaks = c("Brain only", "Peripheral only", "Brain + peripheral")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_x_discrete(labels = x_labs) +
  scale_y_continuous(labels = scales::percent, breaks = c(0, .5, 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.02), clip = "off") +
  labs(x = NULL, y = "Proteins with GTEx support") +
  theme_bw(base_size = 7) + theme(panel.grid = element_blank(), panel.border = element_blank(),
    axis.line.y = element_line(linewidth = .3)) + axes_fig +
  theme(axis.text.x = element_text(size = 6.5, lineheight = .9))

# ---- a: brain-context heatmap ------------------------------------------------------------
# Head-to-head proteins (cis-pQTL in both fluids, 03 cis_both), colocalizations as in 04
# (retained, cis, cognate, PP4 > 0.8).
h2h <- fread(file.path(IN_DIR, "03_fluid_specificity_protein.csv"))[cis_both == TRUE, unique(gene)]
ev <- fread(file.path(IN_DIR, "04_protein_tissue_pair_evidence.csv.gz"))
sig <- unique(ev[is_cis == TRUE & cognate == TRUE & coloc_retained == TRUE & PP4 > COLOC_PP4 &
                 gene %in% h2h, .(gene, fluid, resource, tissue_class, tissue)])
# Must reproduce the head-to-head counts 04 exports.
db <- fread(file.path(IN_DIR, "04_tissue_dumbbell.csv"))[analysis == "head_to_head"]
chk <- merge(db, sig[, .(n2 = uniqueN(gene)), by = .(fluid, resource, tissue_class, tissue)],
             by = c("fluid", "resource", "tissue_class", "tissue"), all.x = TRUE)
chk[is.na(n2), n2 := 0L]
stopifnot(nrow(chk) == nrow(db), all(chk$n == chk$n2))

STATES <- c(`Both fluids` = "grey35", `CSF only` = FLUID_COL[["CSF"]], `Plasma only` = FLUID_COL[["Plasma"]])
MB <- c(cortex = "Cortex", cerebellum = "Cerebellum", hippocampus = "Hippocampus",
        basalganglia = "Basal ganglia", spinalcord = "Spinal cord")
SB <- CELLTYPE_LAB[c("Ext", "IN", "Ast", "OD", "OPC", "MG", "End")]   # row order of the panel
ROWS <- c("GTEx peripheral", MB, SB)
GROUP <- setNames(c("Reference", rep("MetaBrain", length(MB)), rep("SingleBrain", length(SB))), ROWS)
GRP_LV <- c("Reference", "MetaBrain", "SingleBrain")
brain <- sig[resource %in% c("MetaBrain", "SingleBrain")]
brain[, ctx := fifelse(resource == "MetaBrain", MB[tissue], SB[tissue])]
stopifnot(!anyNA(brain$ctx))
periph <- sig[resource == "GTEx" & tissue_class == "Peripheral"][, ctx := "GTEx peripheral"]
st <- rbind(brain, periph)[, .(state = fifelse(uniqueN(fluid) == 2L, "Both fluids",
                                fifelse(fluid[1] == "CSF", "CSF only", "Plasma only"))), by = .(gene, ctx)]

# Columns: the brain-colocalizing proteins, grouped by their most frequent state across
# the 12 brain rows (ties to "Mostly shared"), then ordered row by row so blocks read
# left to right.
genes <- unique(brain$gene)
gs <- st[ctx != "GTEx peripheral", .(nb = sum(state == "Both fluids"), nc = sum(state == "CSF only"),
                                     np = sum(state == "Plasma only")), by = gene]
gs[, grp := fifelse(nb >= pmax(nc, np), "Mostly shared", fifelse(nc > np, "Mostly CSF",
             fifelse(np > nc, "Mostly plasma", "Mostly shared")))]
code <- merge(gs, dcast(st[gene %in% genes], gene ~ factor(ctx, levels = ROWS),
                        value.var = "state", drop = FALSE), by = "gene")
rank_cols <- intersect(c(MB, SB), names(code))
for (r in rank_cols) set(code, j = r, value = match(code[[r]], names(STATES), nomatch = 9L))
code[, grp := factor(grp, levels = c("Mostly shared", "Mostly CSF", "Mostly plasma"))]
setorderv(code, c("grp", rank_cols))
gene_lv <- code$gene
stopifnot(length(gene_lv) == length(genes))

tiles <- st[gene %in% genes]
tiles[, `:=`(gene = factor(gene, levels = gene_lv), ctx = factor(ctx, levels = rev(ROWS)),
             grp = factor(GROUP[as.character(ctx)], levels = GRP_LV))]
grid <- CJ(gene = factor(gene_lv, levels = gene_lv), ctx = factor(ROWS, levels = rev(ROWS)))
grid[, grp := factor(GROUP[as.character(ctx)], levels = GRP_LV)]
breaks_at <- cumsum(table(code$grp))[1:2] + 0.5
glab <- data.table(x = c(0.5, breaks_at) + diff(c(0.5, breaks_at, length(gene_lv) + 0.5)) / 2, y = -Inf,
                   label = paste0(levels(code$grp), " (", as.integer(table(code$grp)), ")"),
                   grp = factor("SingleBrain", levels = GRP_LV))
# Only the proteins the Results name are labelled; per-protein colocalizations are in S3.
NAMED <- c("C1QA", "CD14", "CNTN2", "CTSH", "FAS", "GRN", "PILRB")
stopifnot(all(NAMED %in% gene_lv))

heat <- ggplot() +
  geom_tile(data = grid, aes(gene, ctx), fill = "grey95", colour = NA, height = 0.78) +
  geom_tile(data = tiles, aes(gene, ctx, fill = state), colour = NA, height = 0.78) +
  geom_vline(xintercept = breaks_at, colour = "black", linewidth = 0.3) +
  scale_fill_manual(values = STATES, name = NULL) +
  scale_x_discrete(breaks = NAMED, position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(grp = c(Reference = "", MetaBrain = "MetaBrain", SingleBrain = "SingleBrain"))) +
  geom_text(data = glab, aes(x, y, label = label), vjust = 1.5, size = pt_size(6), inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 7) +
  theme(axis.text.x.top = element_text(size = 5.5, angle = 90, hjust = 0, vjust = 0.5, face = "italic"),
        axis.ticks.x = element_line(linewidth = 0.2), axis.text.y = element_text(size = 6),
        panel.grid = element_blank(), panel.spacing.y = unit(2, "pt"),
        strip.text.y.left = element_text(size = 6, angle = 90), strip.background = element_blank(),
        strip.placement = "outside")

# Margin: share of colocalizing proteins in each state, over the same head-to-head set.
# The GTEx row counts every head-to-head protein colocalizing there, including those
# with no brain colocalization and so no column (the caption says so).
mg <- st[, .N, by = .(ctx, state)]
mg[, `:=`(tot = sum(N), frac = N / sum(N)), by = ctx]
mg[, `:=`(ctx = factor(ctx, levels = rev(ROWS)), state = factor(state, levels = rev(names(STATES))),
          grp = factor(GROUP[as.character(ctx)], levels = GRP_LV))]
tot <- unique(mg[, .(ctx, grp, tot)])
bars <- ggplot(mg, aes(frac, ctx, fill = state)) +
  geom_col(width = 0.62) +
  geom_text(data = tot, aes(x = 1.02, y = ctx, label = tot), inherit.aes = FALSE,
            hjust = 0, size = pt_size(5.5)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.2, colour = "grey50") +
  scale_fill_manual(values = STATES, guide = "none") +
  scale_x_continuous(labels = scales::percent, breaks = c(.5, 1), expand = c(0, 0)) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  labs(x = "% of colocalizing proteins", y = NULL) +
  theme_bw(base_size = 7) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),
        panel.spacing.y = unit(2, "pt"), strip.text.y = element_blank(),
        strip.background = element_blank(), plot.margin = margin(l = 14, r = 12))

# ---- c-e: loci (10e) -------------------------------------------------------------------
# LD bins, colours and lead colour from lib/fig_theme.R, shared with Fig. 6 (05g). The y
# axis runs to at least 10 and carries the genome-wide line, so a flat track reads as flat.
lt <- fread(file.path(IN_DIR, "10e_locus_tracks.csv.gz"))
stopifnot(setequal(unique(lt$gene), c("GRN", "ROR1", "THY1")))
lt[, `:=`(mlp = -log10(pmax(p, 1e-300)),
          ld = cut(r2, LD_BREAKS, LD_LABELS, right = FALSE),
          # missing r2 (variant not in the LD panel) sits in the background with low LD
          point_class = fifelse(is.na(r2) | r2 < 0.2, "Low LD", "Other LD"))]

locus <- function(g, tag, legend = FALSE, show_y_title = TRUE) {
  d <- lt[gene == g]
  d[, track := factor(track, levels = c(setdiff(unique(track), c("CSF pQTL", "Plasma pQTL")),
                                        "CSF pQTL", "Plasma pQTL"))]
  setorder(d, r2, na.last = FALSE)                     # high-LD points drawn on top
  lead_pos <- d[lead == TRUE, pos[1]]
  top <- d[lead == TRUE & p < 5e-8]                    # diamond only where the lead is significant
  floor_y <- d[, .(pos = lead_pos, mlp = 10), by = track]   # invisible: y axis to at least 10
  track_lab <- unique(d[, .(track)])
  track_lab[, label := sub(" ", "\n", as.character(track), fixed = TRUE)]
  ggplot(d, aes(pos / 1e6, mlp)) +
    geom_blank(data = floor_y) +
    GW_LINE +
    # full height in every track, flat ones included, so the lead lines up across them
    geom_vline(xintercept = lead_pos / 1e6, colour = LEAD_LINE_COL, linewidth = 0.15, alpha = 0.55) +
    # Left-aligned at the first variant, a fixed step inside the axis: an hjust offset
    # scales with the label width and put "Cerebellum" against the top y tick.
    geom_label(data = track_lab, aes(label = label), x = min(d$pos) / 1e6, y = Inf,
               hjust = 0, vjust = 1, size = pt_size(6), lineheight = .85,
               fill = alpha("white", .88), linewidth = 0, label.padding = unit(.8, "pt"),
               inherit.aes = FALSE) +
    geom_point(aes(colour = ld, size = point_class, alpha = point_class), stroke = 0) +
    geom_point(data = top, shape = 23, size = 1.65, fill = LEAD_COL, colour = "black", stroke = 0.2) +
    scale_colour_manual(values = setNames(LD_COLS, LD_LABELS), breaks = LD_LABELS,
                        na.value = "grey80", name = expression(italic(r)^2), drop = FALSE) +
    scale_size_manual(values = c(`Low LD` = 0.95, `Other LD` = 1.1), guide = "none") +
    scale_alpha_manual(values = c(`Low LD` = 0.65, `Other LD` = 0.85), guide = "none") +
    guides(colour = if (legend) guide_legend(
      nrow = 1, title.position = "left", title.hjust = 0.5,
      override.aes = list(size = 1.5, alpha = 1)
    ) else "none") +
    facet_grid(track ~ ., scales = "free_y") +
    scale_y_continuous(breaks = function(l) pretty(l, n = 2),
                       expand = expansion(mult = c(.02, .22))) +
    labs(x = paste0("Chromosome ", d$chr[1], " (Mb)"),
         y = if (show_y_title) expression(-log[10]*italic(P)) else NULL, tag = tag,
         title = bquote(italic(.(g))*": "*.(d$label[1]))) +
    theme_classic(base_size = 7) +
    theme(panel.grid = element_blank(), strip.text = element_blank(), strip.background = element_blank(),
          panel.spacing.y = unit(2, "pt"), axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3)) +
    locus_axes
}
