#' @title cellcycle_bar
#' @description create bar graph of % of cells in each cell cycle phase across different Seurat object idents.
#' @param seurat_obj Seurat object.
#' @param group_by the ident from your Seurat object that you want to group the cells by.
#' @export
#' @return bar graph.
#'

cellcycle_bar <- function(seurat_obj, group_by) {
  col <- c("#F1ECEC", "#65DBBD", "#DB6583")
  obj_name <- deparse(substitute(seurat_obj))

  Idents(seurat_obj) <- group_by

  idents <- seurat_obj[[]] %>% pull(group_by) %>% unique()

  plot_data <- NULL
  for (ident in idents) {
    data <- subset(seurat_obj, idents = ident)
    count <- count(data[[]], Phase)
    count_per <- round(100*count$n/sum(count$n),1)
    plot_data[[ident]] <- data.frame(ident = ident,
                                     count_per = count_per,
                                     phase = c("G1", "G2M", "S"))
  }
  plot_data <- data.table::rbindlist(plot_data)

  # plot
  plot <- ggplot(plot_data, aes(x = ident, y = count_per, fill = phase)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col) +
    theme_minimal() +
    ylab("% of cells in phase") +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 24, color = "black"),
          legend.text = element_text(size = 20, color = "black"),
          legend.title = element_blank(),
          axis.title.y = element_text(size = 24, color = "black"),
          axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5))

  return(plot)
}
