# Shared look of the main figures. Nature's figure guide: Helvetica or Arial, text 5-7 pt,
# panel letters 8-pt bold lowercase, figures drawn at their final size (183 mm double
# column, at most 170 mm tall) and supplied as vector files. Each main-figure script
# finishes with
#   fig <- fig & theme_pgd()
#   save_fig(fig, OUT_DIR, "name", height_mm)
# theme_pgd() sets only elements no panel blanks (font family, legend, strips, titles,
# tags): adding an element_text over an element_blank would bring a hidden axis back
# (theme_void brain maps, a table panel). Axis text and titles therefore come from each
# script's theme_bw(base_size = 7) (axis text 5.6 pt) or explicit sizes within 5-7 pt.
# Needs ggplot2, loaded by the caller.

FIG_FAMILY   <- "Helvetica"
FIG_W_MM     <- 183
FIG_H_MAX_MM <- 170
# Supplementary figures are not placed beside text, so they may fill a page.
FIG_H_MAX_SUPP_MM <- 255

# geom_text / annotate size (mm) for a given point size.
pt_size <- function(pt) pt / ggplot2::.pt

theme_pgd <- function() ggplot2::theme(
  text          = ggplot2::element_text(family = FIG_FAMILY),
  legend.text   = ggplot2::element_text(size = 6),
  legend.title  = ggplot2::element_text(size = 7),
  strip.text    = ggplot2::element_text(size = 7),
  plot.title    = ggplot2::element_text(size = 7, face = "plain"),
  plot.subtitle = ggplot2::element_text(size = 6),
  plot.tag      = ggplot2::element_text(size = 8, face = "bold"),
  legend.key.size = ggplot2::unit(8, "pt")
)

# Italic P for plotmath labels (geom_text(parse = TRUE)); small values as m x 10^e.
p_expr <- function(p) vapply(p, function(x) {
  if (x >= 1e-3) return(sprintf('italic(P)=="%s"', formatC(x, digits = 2, format = "fg")))
  e <- floor(log10(x))
  sprintf('italic(P)==%.1f%%*%%10^%d', x / 10^e, e)
}, "")

# PNG only, at the final print size; the HTML manuscript shows it (no PDFs, 2026-10-29).
save_fig <- function(fig, out_dir, name, height_mm, width_mm = FIG_W_MM,
                     dpi = 300, raster_device = c("cairo", "ragg"), max_h_mm = FIG_H_MAX_MM) {
  raster_device <- match.arg(raster_device)
  if (height_mm > max_h_mm) stop(name, ": ", height_mm, " mm is taller than ", max_h_mm, " mm")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(out_dir, paste0(name, ".png"))
  # Cairo draws Unicode ("≥") and falls back to another font for glyphs Helvetica lacks ("→").
  if (raster_device == "ragg") {
    ggplot2::ggsave(f, fig, width = width_mm, height = height_mm, units = "mm",
                   dpi = dpi, device = ragg::agg_png, bg = "white")
  } else {
    ggplot2::ggsave(f, fig, width = width_mm, height = height_mm, units = "mm",
                   dpi = dpi, device = "png", type = "cairo")
  }
  if (!file.exists(f) || file.info(f)$size == 0) stop("not written: ", f)
  message("Saved: ", f, " (", width_mm, " x ", height_mm, " mm)")
}

# A key drawn inside a panel: white box with a thin black frame, as the in-panel labels of
# Fig. 6c,d, so points behind it do not show through and it reads as a key, not data.
legend_box <- ggplot2::theme(
  legend.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.25),
  legend.margin = ggplot2::margin(2, 3, 2, 3))

# Axis text and titles for the multi-panel figures (03c, 10d): 6 and 7 pt, inside the 5-7 pt range.
axes_fig <- ggplot2::theme(
  axis.text.x = ggplot2::element_text(size = 6), axis.text.y = ggplot2::element_text(size = 6),
  axis.title.x = ggplot2::element_text(size = 7), axis.title.y = ggplot2::element_text(size = 7))

# Every colour code of the figures is defined here and nowhere else, so a palette can be
# changed in one place and two figures cannot drift apart (check_fig_style.sh enforces it).
# Orange means plasma and blue means CSF everywhere, so no other palette uses them.
FLUID_COL <- c(CSF = "#0072B2", Plasma = "#D55E00")
# Diverging beta on the brain maps (lib/brain_maps.R) and the enrichment gradient of the
# hotspot gene sets (01c).
BETA_COL <- c(low = "#3B6BB0", mid = "white", high = "#E8735A")
# Paired-effect classes of single-fluid signals (10a); shared/unresolved are greys, divergent
# the beta red and concordant a light beta blue (as novel and matched in Fig. 3a). The
# order matters: 10_plot.R and 12c take names(CLASS_COL) as the factor levels, which set
# the stacking and drawing order.
# Mid and light grey, so the coloured classes stay visible in the Fig. 4a scatter too.
CLASS_COL <- c(shared = "grey60", present = "#A9BEDD", `effect-specific` = BETA_COL[["high"]],
               unresolved = "grey85")
# SomaScan match of a plasma Olink lead (13, Fig. 2c): the panel is all plasma, so a shared
# signal takes the plasma orange (full for the identical lead, light for the same locus) and
# the rest are greys, not the CSF blue. Names are 13's comparison levels, in stacking order.
PLATFORM_COL <- c(`exact lead` = FLUID_COL[["Plasma"]], `same locus` = "#EAAE7F",
                  `different locus` = "grey60", `no SomaScan signal` = "grey90")
# Brain eQTL support: the beta red (light for brain + peripheral), as for divergent signals.
BRAIN_COL <- c(specific = BETA_COL[["high"]], any = "#F4B9AC")
# Novel versus previously reported cis-pQTLs (09f): the beta red for novel, as on the brain
# maps of Fig. 6.
NOVEL_COL <- c(novel = BETA_COL[["high"]], reference = "#7570B3")
# Primary cis-signal classes (09f, Fig. 3): light beta blue for reported, red for unmatched
# (light red for another reported signal), in the order of 09f's CLS table (known, detected
# in Western, unresolved, B, A, C).
NOVELTY_CLS_COL <- c(known = CLASS_COL[["present"]], detected = "#7B8FCB", unresolved = "#C9CEDF",
                     B = BRAIN_COL[["any"]], A = NOVEL_COL[["novel"]], C = "#E8B04B")
# Alzheimer's disease stage at which a protein's CSF level changes (12e, Fig. 6a): muted
# warm reds, light to dark from early to late (Core between the two).
STAGE_COL <- c(Early = "#E9A17F", Core = "#CD6859", Late = "#B94B42")
# Drug evidence of supported pairs (05c, Fig. 5d), as the opacity of each bar's own colour
# (CSF only blue, plasma only orange, both grey): full for a drug for the same disease,
# half for an approved drug for another disease only, faint for neither (2026-11-13; plain
# greys looked uncoloured, a separate hue appeared nowhere else). Names are 05c's parts in
# stacking order.
DRUG_ALPHA <- c(`Same disease` = 1, `Other disease` = 0.5, Neither = 0.12)
# Most severe VEP consequence (03c). Grouped by hue: protein-altering alone in red, splice
# and UTR oranges, synonymous and non-coding exon greys, the non-coding positions blues
# that differ clearly in lightness.
VEP_COL <- c("Protein-altering" = "#b2182b", "Splice region" = "#e6550d", "5' UTR" = "#fdae6b",
             "3' UTR" = "#fdd0a2", "Synonymous" = "#969696", "Non-coding exon" = "#d9d9d9",
             "Up/downstream" = "#c6dbef", "Intronic" = "#6baed6", "Intergenic" = "#08519c")
# Lower saturation for the merged regulatory figure; retain the PAV highlight.
VEP_SOFT_COL <- vapply(VEP_COL, function(col) {
  z <- grDevices::col2rgb(col) / 255
  grDevices::rgb(z[1] * .7 + .3, z[2] * .7 + .3, z[3] * .7 + .3)
}, character(1))
VEP_SOFT_COL[["Protein-altering"]] <- VEP_COL[["Protein-altering"]]

# Compact multi-panel typography, shared by the novelty and integrated regulatory figures.
theme_compact <- function() ggplot2::theme_classic(base_size = 7, base_family = FIG_FAMILY) +
  ggplot2::theme(text = ggplot2::element_text(colour = "grey15"),
    axis.text = ggplot2::element_text(size = 6.2, colour = "grey15"),
    axis.title = ggplot2::element_text(size = 7),
    axis.line = ggplot2::element_line(linewidth = .25, colour = "grey45"),
    axis.ticks = ggplot2::element_line(linewidth = .25, colour = "grey45"),
    axis.ticks.length = ggplot2::unit(1.5, "pt"),
    plot.title = ggplot2::element_text(size = 7.2, face = "plain", hjust = 0,
                                      margin = ggplot2::margin(b = 4)),
    plot.subtitle = ggplot2::element_text(size = 6, colour = "grey40", margin = ggplot2::margin(b = 4)),
    plot.tag = ggplot2::element_text(size = 8, face = "bold"), plot.tag.position = "topleft",
    legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 6.5),
    legend.key.size = ggplot2::unit(6, "pt"), legend.key.spacing.x = ggplot2::unit(2, "pt"),
    legend.key.spacing.y = ggplot2::unit(1, "pt"), legend.spacing.x = ggplot2::unit(2, "pt"),
    legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = ggplot2::unit(2, "pt"),
    strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 6.5, face = "plain"),
    panel.grid = ggplot2::element_blank(), plot.margin = ggplot2::margin(3, 4, 3, 3))
# Protein-level tertiles in the clinical figures (12d, 12e): the beta ends, low blue and high
# red, with a grey middle, so Fig. 6 reads blue as lower and red as higher throughout.
LEVEL_COL <- c(Low = BETA_COL[["low"]], Middle = "grey60", High = BETA_COL[["high"]])
# Gene-set source in the hotspot enrichment panel (Supplementary Fig. 1b): the three GO
# branches as shades of one green, Reactome and UniProt keywords in hues no other figure uses
# for a meaning (blue/orange are the fluids, red/blue the beta classes).
SOURCE_COL <- c(`GO BP` = "#1B7837", `GO CC` = "#5AAE61", `GO MF` = "#A6DBA0",
                Reactome = "#E8B04B", UniProt = "#8C510A")

# Regional association plots (Fig. 5c-e, 10d; Fig. 6b,c, 05g): LD r2 in five bins, the
# LocusZoom scheme, and the purple lead / credible-set diamond. One definition, so the
# two figures cannot drift apart.
LD_BREAKS <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0001)
LD_LABELS <- c("0\u20130.2", "0.2\u20130.4", "0.4\u20130.6", "0.6\u20130.8", "0.8\u20131")
LD_COLS   <- c("#2C3E8C", "#7FB8DA", "#5CB85C", "#EEA236", "#D43F3A")
LEAD_COL  <- "#7B3294"
LEAD_LINE_COL <- "#E57373"
GW_LINE <- ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = 3, colour = "grey55",
                               linewidth = .3)

# The theme sets the font of titles and axes only; text layers (P labels, the region
# names lib/brain_maps.R adds) keep the default font unless given one (12d, 12e).
font_layers <- function(p) {
  for (l in p$layers) if (inherits(l$geom, c("GeomText", "GeomLabel", "GeomTextRepel", "GeomLabelRepel")))
    l$aes_params$family <- FIG_FAMILY
  p
}
