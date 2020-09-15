#' @title nice_spatial_feature_plot
#' @description nicer looking version of the Seurat function SpatialFeaturePlot.
#' @param seurat_obj Seurat object.
#' @param features genes to plot.
#' @param im_alpha set to 1 to plot the tissue image, 0 otherwise.
#' @param pt.size size of plotted points.
#' @param stroke linewidth to outline plotted points in black (default = no outline).
#' @param new_cmap if T, change the colourmap from Seurat default.
#' @param cmap colourmap to use if new_cmap = T.
#' @param n_col if plotting multiple genes, number of columns in final plot_grid object.
#' @param diverging_cmap if T, will fill points with diverging colours between blue and red. otherwise, will use a linear colourmap.
#' @param alpha if T, will fill points such that lower expression values are more transparent.
#' @param cols optional: gradient vector of colours to use for plotting.
#' @export
#' @return SpatialFeaturePlot.

nice_spatial_feature_plot <- function(seurat_obj, features, im_alpha = 0, pt_size = 1.4, stroke = 0, new_cmap = T, cmap = "inferno", n_col = NULL, diverging_cmap = T, alpha = T, cols = NULL) {

  if (length(features) == 1 & is.null(n_col)) {
    n_col <- 1
  }

  if (alpha) {
    alpha <- c(0.1,1)
  } else {
    alpha <- c(1,1)
  }

  if (new_cmap) {

    if (!is.null(cols)) {
      cols <- cols
    } else {
      if (diverging_cmap) {
        cols <- pals::brewer.rdbu(n = 100) %>% rev()
      } else {
        cols <- pals::magma(n = 100)
      }
    }

    plots <- Seurat::SpatialPlot(seurat_obj,
                                 features = features,
                                 alpha = alpha,
                                 image.alpha = im_alpha,
                                 stroke = stroke,
                                 pt.size.factor = pt_size,
                                 combine = F)
    plots <- lapply(plots, function(x) x +
                      scale_fill_gradientn(colours = cols) +
                      theme(legend.title = element_text(face = "bold", size = 16)))

  } else if (new_cmap == F) {

    plots <- Seurat::SpatialPlot(seurat_obj,
                                 features = features,
                                 alpha = alpha,
                                 image.alpha = im_alpha,
                                 stroke = stroke,
                                 pt.size.factor = pt_size,
                                 combine = F)
    plots <- lapply(plots, function(x) x +
                      theme(legend.title = element_text(face = "bold", size = 16)))
  }
  plots <- plot_grid(plotlist = plots, ncol = n_col)
  return(plots)
}
