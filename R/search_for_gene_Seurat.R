#' @title search_for_gene_Seurat
#' @description small function to search for presence of a gene name or string in the rownames of a Seurat object.
#' @param seurat_obj Seurat object.
#' @param gene gene or portion of a gene name you want to search for
#' @param fixed if T, searches for exact matches only.
#' @export
#' @return matching gene names.
#'

search_for_gene_Seurat <- function(seurat_obj, gene, fixed = F) {

  search_result <- grep(pattern = gene, x = rownames(seurat_obj), fixed = fixed, ignore.case = T, value = T)
  return(search_result)

}
