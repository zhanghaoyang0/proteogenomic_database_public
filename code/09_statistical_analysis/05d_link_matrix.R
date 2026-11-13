#!/usr/bin/env Rscript
# Draft panel for the redesigned protein-disease figure: a matrix of jointly supported
# protein-disease links with MR and colocalization evidence in each fluid side by side,
# so which fluid carries a link -- and why the other does not -- reads off the pattern.
# The panel shows only that difference; drug and novelty findings are described in the
# text and other panels, not marked here.
#
# Rows (agreed 2026-10-27): every protein-disease pair that colocalizes (PP4 > 0.8) in
# either fluid for a neurological or psychiatric FinnGen endpoint (the disease area
# closest to a CSF study). Coloc, not MR, is the entry rule: MR alone can ride on LD
# (the APOE region), and the MR column then shows which of these MR also supports.
# One block per disease with the disease as its heading, one row per protein.
# All links in the same layout belong in an Extended Data figure.
#
# Symbols, one column per fluid (after the COLOC / MR symbol matrix of the brain
# cell-type eQTL figure the author pointed to), coloured by fluid:
#   circle            colocalization, PP4 > 0.8
#   triangle up/down  MR estimate positive / negative at FDR < 0.05 (FDR only, the
#                     threshold the text uses)
#   both symbols = joint support in that fluid; empty = tested, not significant
#   grey x            protein not measured in that fluid
#   grey -            measured but no cis instrument, or the pair could not be tested
# The symbols are explained in the figure legend, not inside the figure.
#
# Reads only tables written by 05a and 00 (analysed assays), never the secure MR/coloc
# files. output/05a_disease_pair_evidence.csv.gz is 05a's gitignored pair table: rerun 05a if missing.
#
# Panel a of 05_protein_disease.png: 05c_plot.R sources this file with PANEL_EXPORT_ONLY
# set and takes p_matrix. Run on its own it also saves a preview.
# Outputs:
#   output/05d_link_matrix_rows.csv  one row per drawn link and its four cells (read by
#                                    the manuscript for the coloc-only pairs)
#   tmp/proteogenomic_database/05d_link_matrix.png  standalone preview, repo root (not run when sourced)
# Usage: Rscript code/09_statistical_analysis/05d_link_matrix.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source("code/09_statistical_analysis/lib/fig_theme.R")
source("code/09_statistical_analysis/lib/common.R")   # COLOC_PP4
COLOC_LAB <- sprintf("PP4 > %g", COLOC_PP4)       # the coloc cell label
PGD <- "."
OUT_DIR <- "tmp/proteogenomic_database"
dir.create(OUT_DIR, showWarnings = FALSE)

pair_file <- file.path(PGD, "output", "05a_disease_pair_evidence.csv.gz")
if (!file.exists(pair_file)) stop("missing ", pair_file, "; rerun 05a_disease_data.R")

pair  <- fread(pair_file, select = c("gene", "outcome_id", "disease", "disease_category",
                                     "fluid", "mr_ok", "mr_sig", "ivw_b", "ivw_p",
                                     "coloc_tested", "coloc_sig", "joint_sig"))

# Proteins measured per fluid: the analysed assays (00_dedup_assays.R), as 05a uses them.
kept <- fread(file.path(PGD, "output", "00_kept_assays.csv"))
measured <- unique(kept[platform != "SomaLogic", .(gene, fluid)])
stopifnot(measured[, uniqueN(gene), by = fluid][, all(V1 > 1000)])

# ---- rows ----------------------------------------------------------------------
links <- unique(pair[disease_category %in% c("Neurological", "Psychiatric") &
                       coloc_sig %in% TRUE, .(gene, outcome_id, disease)])
stopifnot(nrow(links) > 0)

short_disease <- function(x) {
  x <- sub(", including avohilmo", "", x)
  x <- sub("^MS-disease / ", "", x)
  x <- sub("Other and unspecified polyneuropathies, also in other diseases",
           "Polyneuropathy", x)
  x <- sub("Alzheimer’s disease, wide definition", "Alzheimer's disease (wide definition)", x)
  # FinnGen writes "Alzheimer disease"; the manuscript uses "Alzheimer's disease".
  x <- sub("Alzheimer disease", "Alzheimer's disease", x, fixed = TRUE)
  x <- sub("Multiple Sclerosis", "Multiple sclerosis", x)
  x
}
links[, sd := short_disease(disease)]
# One row per protein and disease: FinnGen variants of one endpoint (wide, AD dementia)
# collapse to one row, the definition with the most evidence (fluids with significant
# MR plus fluids with coloc), ties to the one without a qualifier.
links[, base := sub(" \\(.*\\)$", "", sd)]
links[base == "Dementia in Alzheimer's disease", base := "Alzheimer's disease"]
score <- unique(pair[, .(gene, outcome_id, fluid, mr = mr_sig %in% TRUE,
                         co = coloc_sig %in% TRUE)])[
  , .(score = sum(mr) + sum(co)), by = .(gene, outcome_id)]
links <- merge(links, score, by = c("gene", "outcome_id"))
links[, keep := seq_len(.N) == order(-score, sd != base, nchar(sd))[1], by = .(gene, base)]
message("collapsed to one definition: ",
        paste(links[keep == FALSE, paste(gene, sd)], collapse = "; "))
links <- links[keep == TRUE][, c("keep", "score") := NULL]
links[, key := paste(gene, outcome_id)]
# the row shows the protein; a qualifier the base heading drops stays on the row
qualifier <- function(sd, base) fcase(
  sd == base, "",
  grepl("\\)$", sd), sub("^.*( \\(.*\\))$", "\\1", sd),
  sd == "Dementia in Alzheimer's disease", " (AD dementia)",
  default = paste0(" (", sd, ")"))
links[, label := paste0(gene, qualifier(sd, base))]

# ---- cells ---------------------------------------------------------------------
# One assay per protein and fluid is analysed (00_kept_assays), but keep the
# strongest row if a complex maps one gene twice: joint first, then smallest MR P.
setorder(pair, -joint_sig, ivw_p)
ev <- unique(pair[paste(gene, outcome_id) %in% links$key],
             by = c("gene", "outcome_id", "fluid"))
cells <- merge(links[, .(fluid = c("CSF", "Plasma")), by = .(gene, outcome_id)],
               ev, by = c("gene", "outcome_id", "fluid"), all.x = TRUE)
n0 <- nrow(cells)
cells[, is_measured := paste(gene, fluid) %in% measured[, paste(gene, fluid)]]
stopifnot(nrow(cells) == n0, cells[!is.na(mr_ok), all(is_measured)])

cells[, mr := fcase(
  !is_measured, "not measured",
  is.na(mr_ok) | !mr_ok, "not testable",
  mr_sig %in% TRUE & ivw_b > 0, "MR positive",
  mr_sig %in% TRUE & ivw_b < 0, "MR negative",
  default = "tested, not significant")]
cells[, coloc := fcase(
  !is_measured, "not measured",
  is.na(coloc_tested) | !coloc_tested, "not testable",
  coloc_sig %in% TRUE, COLOC_LAB,
  default = "tested, not significant")]
cells[, joint := joint_sig %in% TRUE]
# the joint flag is defined by these two cells; check they agree
stopifnot(cells[joint == TRUE, all(grepl("^MR", mr) & coloc == COLOC_LAB)])
cells[, key := paste(gene, outcome_id)]

# One row per drawn mark: each fluid has a Coloc and an MR sub-column (centres at
# COLOC_X and MR_X); a grey x or - fills both when the fluid cannot test the link.
COLOC_X <- -0.5
MR_X    <- 0.5
long <- rbind(
  cells[coloc == COLOC_LAB, .(key, fluid, x = COLOC_X, sym = "Colocalization")],
  cells[mr %in% c("MR positive", "MR negative"), .(key, fluid, x = MR_X, sym = mr)])
blank <- cells[mr %in% c("not measured", "not testable"),
               .(key, fluid,
                 sym = fifelse(mr == "not measured", "Not measured", "No cis instrument"))]
blank <- blank[, .(x = c(COLOC_X, MR_X)), by = .(key, fluid, sym)]

# ---- row order -------------------------------------------------------------------
# Diseases: the neurodegenerative ones first, then the rest by number of links.
# Within a disease: colocalizing in both fluids, then CSF only, then plasma only.
cells[, has_coloc := coloc == COLOC_LAB]
pat <- dcast(cells, key ~ fluid, value.var = "has_coloc")
pat[, pattern := fcase(CSF & Plasma, 1L, CSF, 2L, default = 3L)]
links <- merge(links, pat[, .(key, pattern)], by = "key")
FIRST <- c("Alzheimer's disease", "Frontotemporal dementia", "Parkinson's disease",
           "Multiple sclerosis")
rest <- links[!base %in% FIRST, .N, by = base][order(-N, base), base]
links[, base := factor(base, c(intersect(FIRST, links$base), rest))]
setorder(links, base, pattern, gene)

# Rows top to bottom: each disease heading and its proteins.
rows <- rbind(
  links[, .(text = c(as.character(base[1]), label),
            type = c("heading", rep("link", .N)), k = c(NA, key)), by = base
        ][, base := NULL])
rows[, y := rev(seq_len(.N))]
ykey <- rows[type == "link", .(key = k, y)]
long  <- merge(long, ykey, by = "key")
blank <- merge(blank, ykey, by = "key")
# every link draws at least its coloc circle in one fluid
stopifnot(uniqueN(long$key) == nrow(links))

# The drawn rows, one per link, for the manuscript: which cells carry evidence.
rows_out <- dcast(cells[key %in% links$key], key ~ fluid, value.var = c("coloc", "mr", "joint"))
rows_out <- merge(links[, .(key, gene, outcome_id, disease = as.character(base), label)],
                  rows_out, by = "key")[, key := NULL]
stopifnot(nrow(rows_out) == nrow(links))
fwrite(rows_out, file.path(PGD, "output", "05d_link_matrix_rows.csv"))

# ---- plot ----------------------------------------------------------------------
# Three facets share the rows: the labels (drawn as text, so headings can be bold)
# and two sub-columns of cells per fluid. LABEL_W is the label column's width in cells.
LABEL_W <- 5.4
COLS <- c(" ", "CSF", "Plasma")
FLUID_COLS <- FLUID_COL
# 4 and 45 are R's x and "-" glyphs, so the two cannot-test marks share the legend
# "No evidence" is never drawn: an empty cell says it. It is in the scale only so the
# legend can show an empty cell as its key (shape 22, a white square with a cell border).
SHAPES <- c("Colocalization" = 16, "MR positive" = 17, "MR negative" = 25,
            "Not measured" = 4, "No cis instrument" = 45, "No evidence" = 22)
BLANK_COL <- "grey45"
long[, `:=`(sym = factor(sym, names(SHAPES)), col = factor(fluid, COLS))]
blank[, `:=`(sym = factor(sym, names(SHAPES)), col = factor(fluid, COLS))]
rows[, col := factor(" ", COLS)]
# Column titles: the fluid on its own row above, Coloc / MR on the first disease
# heading's row, so no row is left empty.
heads <- data.table(col = factor(c("CSF", "Plasma"), COLS), x = 0, y = max(rows$y) + 1,
                    text = c("CSF", "Plasma"), face = "plain")
heads <- rbind(heads, CJ(col = factor(c("CSF", "Plasma"), COLS), x = c(COLOC_X, MR_X))[
  , `:=`(y = max(rows$y), text = fifelse(x == COLOC_X, "Coloc", "MR"), face = "plain")])
# one bordered cell per link, fluid and test, as in the reference matrix
grid <- CJ(y = ykey$y, col = factor(c("CSF", "Plasma"), COLS), x = c(COLOC_X, MR_X))
xlim <- data.table(col = factor(rep(COLS, each = 2), COLS),
                   x = c(0, LABEL_W, -1, 1, -1, 1), y = max(rows$y))

# R draws a triangle around the point's circle, not its own middle, so an up
# triangle sits high and a down one low; nudge each back level with the circles
# (offset about 0.39 x the symbol radius, ~0.08 of a 3 mm row at size 1.6).
TRI_NUDGE <- 0.08
pt <- function(d, dy) geom_point(data = d, aes(x, y, shape = sym, colour = fluid, fill = fluid),
                                 size = 1.6, stroke = 0.3, position = position_nudge(y = dy))
p_matrix <- ggplot() +
  geom_blank(data = xlim, aes(x, y)) +
  geom_text(data = rows, aes(fifelse(type == "heading", 0, 0.3), y, label = text,
                             fontface = fifelse(type == "heading", "bold", "plain")),
            hjust = 0, size = pt_size(6), family = FIG_FAMILY) +
  geom_text(data = heads, aes(x, y, label = text, fontface = face),
            size = pt_size(7), family = FIG_FAMILY) +
  geom_tile(data = grid, aes(x, y), width = 1, height = 1,
            fill = NA, colour = "grey30", linewidth = 0.2) +
  pt(long[sym == "Colocalization"], 0) +
  pt(long[sym == "MR positive"], -TRI_NUDGE) +
  pt(long[sym == "MR negative"], TRI_NUDGE) +
  # R's "-" glyph is drawn far smaller than its x, so it gets a larger size
  geom_point(data = blank, aes(x, y, shape = sym,
                               size = fifelse(sym == "Not measured", 1.4, 3)),
             colour = BLANK_COL, stroke = 0.4,
             # TRUE: key every shape level, also "No evidence", which no layer draws
             # (with NA a layer keys only the values in its data); override.aes in
             # guides() sets each key's look
             show.legend = TRUE) +
  scale_size_identity() +
  facet_grid(. ~ col, scales = "free_x", space = "free_x") +
  # Two rows of three under the matrix (2026-11-12, was above the column headers): evidence
  # marks on the first, the cannot-test marks and the empty cell on the second; at the
  # bottom like the keys of b-d, so the matrix starts at the top of the figure.
  scale_shape_manual(values = SHAPES, name = NULL, drop = FALSE,
                     labels = c("Colocalization", "MR (+)", "MR (\u2212)", "Not measured",
                                "No instrument", "No evidence")) +
  scale_colour_manual(values = FLUID_COLS, guide = "none") +
  scale_fill_manual(values = FLUID_COLS, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, max(rows$y) + 1.5), expand = c(0, 0)) +
  guides(shape = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(
    colour = c(rep("black", 3), rep(BLANK_COL, 2), "grey30"),
    fill = c(rep("black", 3), rep(BLANK_COL, 2), "white"), size = c(1.6, 1.6, 1.6, 1.4, 3, 2.2),
    stroke = c(0.3, 0.3, 0.3, 0.4, 0.4, 0.3)))) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 7) +
  theme(panel.spacing.x = unit(2, "pt"),
        # theme_void() leaves the background transparent
        plot.background = element_rect(fill = "white", colour = NA),
        legend.position = "bottom", legend.justification = "center", legend.location = "plot",
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        legend.text = element_text(size = 5.5))
# theme_pgd() styles strips, so blank them after it: the column titles are text rows
p_matrix <- p_matrix + theme_pgd() + theme(strip.text = element_blank())

if (!exists("PANEL_EXPORT_ONLY"))
  save_fig(p_matrix, OUT_DIR, "05d_link_matrix", height_mm = 17 + 3 * nrow(rows), width_mm = 66)
message("05d_link_matrix done: ", nrow(links), " links over ", uniqueN(links$base),
        " diseases")
