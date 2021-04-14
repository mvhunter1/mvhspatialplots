#' @title convert_motif_targets_human
#' @description convert the genelist from HOMER, i.e genes enriched in a particular motif, to a ranked vector for GSEA, using the original FindMarkers results and adding human gene names. * Improved version!
#' @param motif_genes character vector with the list of genes enriched in your motif of interest.
#' @param expression_data FindMarkers results used to create original genelist for HOMER.
#' @param create_ranked_vector default = T, create ranked vector to use with fgsea?
#' @export
#' @return ranked vector or dataframe with converted human genes.


convert_motif_targets_human <- function(motif_genes, expression_data, create_ranked_vector = T) {

  if (!"gene" %in% colnames(expression_data)) {
    stop("You must have a column with the gene names called \"gene\" in your FindMarkers results to proceed with this function.")
  }

  # with newer versions of FindMarkers the logFC column is called "avg_log2FC" - add a check to change the column name if necessary
  if ("avg_log2FC" %in% colnames(expression_data)) {
    index <- match("avg_log2FC", colnames(expression_data))
    colnames(expression_data)[index] <- "avg_logFC"
  }

  motif_expression_data <- expression_data[expression_data$gene %in% motif_genes,]

  if (create_ranked_vector) {
    motif_expression_data_human <- motif_expression_data %>% convert_FindMarkers_human(., create_ranked_vector = T)
  } else {
    motif_expression_data_human <- motif_expression_data %>% convert_FindMarkers_human(., create_ranked_vector = F)
  }
  return(motif_expression_data_human)
}
