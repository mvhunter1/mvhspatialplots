#' @title wilcox_test_expr
#' @description wilcoxon rank-sum test for multiple groups, to compare expression of a specific gene.
#' @param seurat_obj Seurat object.
#' @param gene gene to compare expression data across groups for.
#' @param group_by groups you want to compare, should be a column in metadata.
#' @export
#' @return P value matrix.

wilcox_test_expr <- function(seurat_obj, gene, group_by) {

  data <- GetAssayData(seurat_obj, slot = "data") %>% .[rownames(.) == gene,] %>% enframe() # extract expression data
  colnames(data) <- c("cell", "expr_data")

  # merge with metadata to extract cluster IDs
  metadata <- seurat_obj[[]] %>% rownames_to_column(var = "cell")
  all_data <- merge(x = data,
                    y = metadata,
                    by = "cell")

  # perform wilcoxon rank-sum test with bonferroni correction
  pairwise.wilcox.test(all_data$expr_data, all_data[,group_by], p.adjust.method = "bonferroni")

}
