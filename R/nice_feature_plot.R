#' @title nice_feature_plot
#' @description nicer looking version of the Seurat function FeaturePlot.
#' @param seurat_obj Seurat object.
#' @param features genes to plot.
#' @param pt.size size of plotted points.
#' @param n_col if plotting multiple genes, number of columns in final plot_grid object.
#' @param reduction one of "umap" or "pca".
#' @param dims_plot number of dimensions to plot.
#' @param diverging_cmap if T, will fill points with diverging colours between blue and red. otherwise, will use a linear colourmap.
#' @param scale_data plot Z-scored data?
#' @param cutoffs expression level cutoffs for plotting.
#' @param order plot cells in order of expression?
#' @param cols optional: colourmap.
#' @export
#' @return UMAP or PCA plot.

nice_feature_plot <- function(seurat_obj, features, pt.size = 1.3, n_col = NULL, reduction = "umap", dims_plot = 1:2, diverging_cmap = F, scale_data = F, cutoffs = NA, order = F) {

  if (reduction == "umap") {
    xlab <- "UMAP 1"
    ylab <- "UMAP 2"
  } else if (reduction == "pca") {
    # determine % variability associated with each PC
    pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100

    # make x and y axis labels
    xlab <- paste0("PC", dims_plot[1], " ", round(pct[dims_plot[1]],2), "%")
    ylab <- paste0("PC", dims_plot[2], " ", round(pct[dims_plot[2]],2), "%")
  }

  if (diverging_cmap) {
    cols <- pals::brewer.rdbu(n = 100) %>% rev()
  } else if (is.null(cols)) {
    cols <- pals::virids(n = 100)
  } else {
    cols <- cols
  }

  slot <- ifelse(scale_data == T, "scale.data", "data")

  if (anyNA(cutoffs)) {
    min_cutoff <- NA
    max_cutoff <- NA
  } else {
    min_cutoff <- cutoffs[1]
    max_cutoff <- cutoffs[2]
  }

  if (length(features) == 1) {

    plot <- Seurat::FeaturePlot(seurat_obj,
                                features = features,
                                pt.size = pt.size,
                                reduction = reduction,
                                slot = slot,
                                min.cutoff = min_cutoff,
                                max.cutoff = max_cutoff,
                                order = order) +
      scale_colour_gradientn(colours = cols) +
      xlab(xlab) +
      ylab(ylab) +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_text(size = 20),
            axis.line = element_line(size = 1))
    return(plot)

  } else if (length(features) > 1) {

    plots <- Seurat::FeaturePlot(seurat_obj,
                                 features = features,
                                 pt.size = pt.size,
                                 combine = F,
                                 reduction = reduction,
                                 slot = slot,
                                 min.cutoff = min_cutoff,
                                 max.cutoff = max_cutoff,
                                 order = order)
    plots <- lapply(plots, function(x) x +
                      xlab(xlab) +
                      ylab(ylab) +
                      scale_colour_gradientn(colours = cols) +
                      theme(axis.ticks = element_blank(),
                            axis.text = element_blank(),
                            axis.title = element_text(size = 20),
                            axis.line = element_line(size = 1)))
    plots_use <- cowplot::plot_grid(plotlist = plots, ncol = n_col)
    return(plots_use)
  } else {
    stop("features are missing")
  }
}
