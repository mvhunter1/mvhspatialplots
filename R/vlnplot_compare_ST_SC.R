#' @title vlnplot_compare_ST_SC
#' @description nicer looking version of the Seurat function VlnPlot, to compare expression of gene(s) in scRNA-seq and ST datasets.
#' @param SC_obj single cell Seurat object.
#' @param ST_obj ST Seurat object.
#' @param features gene(s) you want to plot expression of.
#' @param pt_size point size.
#' @param sort sort X-axis in descending order of mean expression?
#' @param dims_plot dimensions to plot.
#' @param plot_hline plot a dashed horizontal line at 0?
#' @export
#' @return violin plots.

vlnplot_compare_ST_SC <- function(SC_obj, ST_obj, features, pt_size = 0.3, sort = T, plot_hline = T) {

  require(Seurat)
  require(tidyverse)
  require(pals)
  require(cowplot)

  # create colourmap
  n_groups_SC <- Idents(SC_obj) %>% levels() %>% length()
  n_groups_ST <- Idents(ST_obj) %>% levels() %>% length()

  cols_SC <- pals::tol(n = n_groups_SC)
  cols_ST <- pals::tol(n = n_groups_ST)

  # plot
  if (length(features) == 1) {
    SC_plot <- Seurat::VlnPlot(SC_obj,
                               features = features,
                               pt.size = pt_size,
                               sort = sort) +
      Seurat::NoLegend() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 16, margin = margin(r = 10)),
            axis.text = element_text(size = 14)) +
      scale_fill_manual(values = cols_SC) +
      ylab("expression in single cell data")

    ST_plot <- Seurat::VlnPlot(ST_obj,
                               features = features,
                               pt.size = pt_size,
                               sort = sort) +
      Seurat::NoLegend() +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size = 16, margin = margin(r = 10)),
            axis.text = element_text(size = 14)) +
      scale_fill_manual(values = cols_ST) +
      ylab("expression in ST data")

    if (plot_hline) {
      SC_plot <- SC_plot + geom_hline(yintercept = 0, linetype = "dashed")
      ST_plot <- ST_plot + geom_hline(yintercept = 0, linetype = "dashed")
    }
    plots <- plot_grid(SC_plot, ST_plot, nrow = 1, axis = "lrbt", align = "hv")
    return(plots)

    } else {

      SC_plots <- Seurat::VlnPlot(SC_obj,
                                  features = features,
                                  pt.size = pt_size,
                                  sort = sort,
                                  combine = F)
      SC_plots <- lapply(SC_plots, function(x) {
        x + Seurat::NoLegend() +
          theme(axis.title.x = element_blank(),
                                       axis.title.y = element_text(size = 16, margin = margin(r = 10)),
                                       axis.text = element_text(size = 14)) +
          scale_fill_manual(values = cols_ST) +
          ylab("expression in single cell data") +
          geom_hline(yintercept = 0, linetype = "dashed")
        })

      ST_plots <- Seurat::VlnPlot(ST_obj,
                                  features = features,
                                  pt.size = pt_size,
                                  sort = sort,
                                  combine = F)
      ST_plots <- lapply(ST_plots, function(x) {
        x + Seurat::NoLegend() +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_text(size = 16, margin = margin(r = 10)),
                axis.text = element_text(size = 14)) +
          scale_fill_manual(values = cols_ST) +
          ylab("expression in ST data") +
          geom_hline(yintercept = 0, linetype = "dashed")
      })

      SC_plots <- plot_grid(plotlist = SC_plots, ncol = 1, align = "hv", axis = "lrbt")
      ST_plots <- plot_grid(plotlist = ST_plots, ncol = 1, align = "hv", axis = "lrbt")

      all_plots <- plot_grid(SC_plots, ST_plots, ncol = 2)
      return(all_plots)
  }

}
