#' @title nice_violin_plot
#' @description nicer looking version of the Seurat function VlnPlot.
#' @param seurat_obj Seurat object.
#' @param features genes to plot.
#' @param group_by how to group your data, i.e the X axis of the plot - usually a metadata column in the seurat object.
#' @param cols optional: colours to fill the violins with.
#' @param pt.size size of plotted points.
#' @param sort sort X axis in descending order of average expression?
#' @param n_col if multiple features, how many columns to create in the final plot_grid object.
#' @param plot_hline plot a dashed horizontal line at 0?
#' @export
#' @return violin plot.
#'
nice_violin_plot <- function(seurat_obj, features, group_by = NULL, cols = NULL, pt.size = 0.3, sort = T, n_col = NULL, plot_hline = T) {

  if (!is.null(group_by)) {
    if (is.null(cols)) {
      n_cols <- seurat_obj[[]] %>% dplyr::select(all_of(group_by)) %>% unique() %>% nrow()
      if (n_cols <= 12) {
        cols <- pals::tol(n_cols)
      }
    }
  }

  if (length(features) == 1) {
    plot <- Seurat::VlnPlot(seurat_obj,
                            group.by = group_by,
                            features = features,
                            pt.size = pt.size,
                            cols = cols,
                            sort = sort) +
      Seurat::NoLegend() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10))

    if (plot_hline == T) {
      plot <- plot + geom_hline(yintercept = 0, linetype = "dashed")
    }
    return(plot)

  } else {
    plotlist <- Seurat::VlnPlot(seurat_obj,
                                group.by = group_by,
                                features = features,
                                pt.size = pt.size,
                                sort = sort,
                                cols = cols,
                                combine = F)
    plotlist <- lapply(plotlist, function(x)
      x + theme(axis.title.x = element_blank(),
                axis.title.y = element_text(size = 10)) + NoLegend())

    if (plot_hline == T) {
      plotlist <- lapply(plotlist, function(x)
        x + geom_hline(yintercept = 0, linetype = "dashed"))
    }
    plots <- cowplot::plot_grid(plotlist = plotlist, ncol = n_col)
    return(plots)
  }
}
