# Brain-map panels for 12d_plot_examples.R and 12e_plot_clinical_overview.R, styled after Lu et al.'s
# rs5848 Fig. 1: ggseg 1.x atlases (aseg with ventricles, Desikan), std beta in blue-white-
# salmon, names in black boxes for P < 0.05 and grey otherwise, dark grey not tested.
# Needs data.table, ggplot2, ggrepel, sf and ggseg (< 2.0) loaded by the caller.

# ---- 11a region names -> ggseg 1.x labels ----------------------------------------------
to_dk <- function(r) paste0(ifelse(grepl("_L$", r), "lh_", "rh_"), sub("_[LR]$", "", r))
to_aseg <- function(r) {
  out <- ifelse(grepl("^CC_", r), r, gsub("_", "-", r))
  out <- sub("^3rd-Ventricle$", "x3rd-ventricle", out)
  sub("^4th-Ventricle$", "x4th-ventricle", out)
}

# ---- style --------------------------------------------------------------------------------
NOT_TESTED <- "grey40"
theme_brain_fig <- theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0.5, margin = margin(b = 2)),
        legend.position = "bottom", legend.key.height = unit(5, "pt"),
        legend.key.width = unit(30, "pt"), legend.title = element_text(size = 8),
        plot.tag = element_text(face = "bold", size = 11))

# One brain panel: fill by std beta, names at region centroids (taken from the built
# plot, so they follow ggseg's own layout). label_ns = FALSE names only significant regions.
# names: optional data.frame(label, text) -- name only these labels, with this text
# (the composite panel names each lobe once instead of every region in it).
# max_labels: at most this many significant regions are named (the smallest P), so a map
# with many significant regions stays readable.
# position: optional ggseg position_brain(), e.g. hemi ~ side for a 2 x 2 cortical block.
brain <- function(d, atlas, lim, title, fill_lab, label_ns = TRUE, names = NULL, max_labels = 8,
                  position = NULL) {
  d <- data.frame(label = d$label, fill = pmax(pmin(d$std_beta, lim), -lim), p = d$p)
  g <- ggplot() +
    (if (is.null(position)) geom_brain(atlas = atlas, data = d[, c("label", "fill")], mapping = aes(fill = fill),
               colour = "grey25", size = 0.15)
     else geom_brain(atlas = atlas, data = d[, c("label", "fill")], mapping = aes(fill = fill),
               colour = "grey25", size = 0.15, position = position)) +
    scale_fill_gradient2(low = BETA_COL[["low"]], mid = BETA_COL[["mid"]], high = BETA_COL[["high"]], midpoint = 0,   # lib/fig_theme.R
                         limits = c(-lim, lim), na.value = NOT_TESTED, name = fill_lab,
                         # A thin frame round the colour bar, so its pale middle reads as part
                         # of the bar (2026-11-12).
                         guide = guide_colourbar(theme = theme(
                           legend.frame = element_rect(colour = "grey30", linewidth = 0.25))))
  built <- ggplot_build(g)$data[[1]]
  if (!"label" %in% names(built)) stop("ggseg layer lost its labels; cannot place region names")
  # One name per region, in the view where it is largest (a region drawn in the lateral
  # and medial views would otherwise get one centroid between the two).
  want <- if (is.null(names)) d$label else names$label
  b <- st_as_sf(built[!is.na(built$label) & built$label %in% want, c("label", "side", "geometry")])
  b <- aggregate(b, by = list(lab = b$label, sd = b$side), FUN = function(x) x[1])
  b$area <- as.numeric(st_area(b))
  b <- b[order(b$lab, -b$area), ]; b <- b[!duplicated(b$lab), ]
  cen <- suppressWarnings(st_centroid(b))
  xy <- data.frame(label = cen$lab, st_coordinates(cen))
  reg <- if (is.null(names)) unique(as.data.frame(atlas)[, c("label", "region")])
         else data.frame(label = names$label, region = names$text)
  xy <- merge(merge(xy, reg, by = "label"), d, by = "label")
  xy$sig <- xy$p < 0.05
  if (sum(xy$sig) > max_labels) xy$sig <- xy$sig & rank(xy$p, ties.method = "first") <= max_labels
  g +
    { if (label_ns && any(!xy$sig)) geom_text_repel(data = xy[!xy$sig, ], aes(X, Y, label = region),
        size = 1.8, colour = "grey65", segment.size = 0.15, segment.colour = "grey70",
        min.segment.length = 0, max.overlaps = Inf, seed = 1) } +
    { if (any(xy$sig)) geom_label_repel(data = xy[xy$sig, ], aes(X, Y, label = region),
        size = 2.1, colour = "black", fill = "white", label.size = 0.15, label.padding = 0.12,
        segment.size = 0.2, min.segment.length = 0, max.overlaps = Inf, seed = 1) } +
    labs(title = title) + theme_brain_fig
}

