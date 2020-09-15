#' @title barcode_plot
#' @description Create GSEA barcode plot from a chosen pathway.
#' @param pathway the chosen pathway to plot.
#' @param gsea_pathways your list of pathways used for GSEA.
#' @param gsea_vector the ranked vector you used for GSEA with your human gene names and logFC.
#' @param gsea_results your GSEA results dataframe.
#' @export
#' @return A nice looking barcode plot for your chosen pathway with gene ranks.

barcode_plot <- function(pathway, gsea_pathways, gsea_vector, gsea_results) {

  require(limma)

  pathway_plot <- pathway
  pathway_genes <- gsea_pathways[[pathway_plot]]
  nes_plot <- gsea_results[gsea_results$pathway == pathway_plot, ] %>% pull(NES)
  plot <- barcodeplot(statistics = gsea_vector, index = pathway_genes, quantiles = c(0, 0), main = paste0(pathway_plot, "\n", "NES = ", round(nes_plot, digits = 3)), xlab = "average log fold change")
  return(plot)

}
