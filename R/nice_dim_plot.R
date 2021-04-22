#' @title nice_dim_plot
#' @description nicer looking version of the Seurat function DimPlot.
#' @param seurat_obj Seurat object.
#' @param group_by what to colour the points by, usually a column in the Seurat object metadata.
#' @param cols optional: colours to label the grouped points by.
#' @param pt_size point size.
#' @param label label groups with text on plot?
#' @param reduction either "umap" or "pca".
#' @param dims_plot dimensions to plot.
#' @export
#' @return UMAP or PCA plot.

nice_dim_plot <- function(seurat_obj, group_by = NULL, cols = NULL, pt_size = 1.3, label = T, reduction = "umap", dims_plot = 1:2) {

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

  if (length(group_by) == 1) {

    if (is.null(cols)) {
      n_cols <- seurat_obj[[]] %>% dplyr::select(all_of(group_by)) %>% unique() %>% nrow()
      if (n_cols <= 12) {
        cols <- pals::tol(n_cols)
      }
    }

    if (label == T) {
      plot <- Seurat::DimPlot(seurat_obj,
                              group.by = group_by,
                              pt.size = pt_size,
                              reduction = reduction,
                              cols = cols,
                              label = T) +
        NoLegend() +
        xlab(xlab) +
        ylab(ylab) +
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_text(size = 20),
              axis.line = element_line(size = 1))
    } else {
      plot <- Seurat::DimPlot(seurat_obj,
                              group.by = group_by,
                              pt.size = pt_size,
                              reduction = reduction,
                              cols = cols,
                              label = F) +
        xlab(xlab) +
        ylab(ylab) +
        theme(axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_text(size = 20),
              axis.line = element_line(size = 1))
    }
    return(plot)

  } else if (length(group_by) > 1) {

    if (label == T) {
      plots <- Seurat::DimPlot(seurat_obj,
                               group.by = group_by,
                               pt.size = pt_size,
                               combine = F,
                               reduction = reduction,
                               cols = cols,
                               label = T)
      plots <- lapply(plots, function(x) x +
                        NoLegend() +
                        xlab(xlab) +
                        ylab(ylab) +
                        theme(axis.ticks = element_blank(),
                              axis.text = element_blank(),
                              axis.title = element_text(size = 20),
                              axis.line = element_line(size = 1)))
    } else {
      plots <- Seurat::DimPlot(seurat_obj,
                               group.by = group_by,
                               pt.size = pt_size,
                               combine = F,
                               reduction = reduction,
                               cols = cols,
                               label = F)
      plots <- lapply(plots, function(x) x +
                        xlab(xlab) +
                        ylab(ylab) +
                        theme(axis.ticks = element_blank(),
                              axis.text = element_blank(),
                              axis.title = element_text(size = 20),
                              axis.line = element_line(size = 1)))
    }
    plots_use <- cowplot::plot_grid(plotlist = plots, ncol = n_col)
    return(plots_use)
  } else {
    stop("group by is missing")
  }
}




