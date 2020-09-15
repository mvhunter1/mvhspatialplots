#' @title nice_dim_plot
#' @description nicer looking version of the Seurat function DimPlot.
#' @param seurat_obj Seurat object.
#' @param group_by what to colour the points by, usually a column in the Seurat object metadata.
#' @param cols optional: colours to label the grouped points by.
#' @param pt_size point size.
#' @param label label groups with text on plot?
#' @param show_legend if T, will include the legend but won't label points w/text on plot.
#' @param reduction either "umap" or "pca".
#' @param dims_plot dimensions to plot.
#' @export
#' @return UMAP or PCA plot.

nice_dim_plot <- function(seurat_obj, group_by = "seurat_clusters", cols = NULL, pt_size = 1.3, label = T, show_legend = F, reduction = "umap", dims_plot = 1:2) {

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

  if (label == T) {

    plot <- Seurat::DimPlot(seurat_obj,
                            group.by = group_by,
                            cols = cols,
                            dims = dims_plot,
                            reduction = reduction,
                            pt.size = pt_size,
                            label = T,
                            label.size = 6) +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_text(size = 20),
            axis.line = element_line(size = 1)) +
      xlab(xlab) +
      ylab(ylab)

    if (show_legend == F) {
      plot <- plot + NoLegend()
      return(plot)
    } else {
      return(plot)
    }

  } else if (label == F) {

    plot <- Seurat::DimPlot(seurat_obj,
                            group.by = group_by,
                            cols = cols,
                            dims = dims_plot,
                            reduction = reduction,
                            pt.size = pt_size) +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_text(size = 20),
            axis.line = element_line(size = 1)) +
      xlab(xlab) +
      ylab(ylab)
    return(plot)
  }
}
