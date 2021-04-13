#' @title create_HOMER_genelist
#' @description Create txt file of genes with ENSEMBL IDs to be used for HOMER.
#' @param genelist List of genes you want to perform HOMER on. Must include "gene" and "avg_logFC" columns.
#' @param save_as name that HOMER txt file will be saved as, don't have to include txt at end.
#' @param remove_ribosomal_genes default = T, remove ribosomal genes from FindMarkers results?
#' @export
#' @return doesn't return anything into R, just saves txt file in the working directory.


create_HOMER_genelist <- function(genelist, save_as, remove_ribosomal_genes = T) {

  # genelist must have a column named "gene" with the gene names otherwise this function won't work.
  if (!"gene" %in% colnames(genelist)) {
    stop("You must have a column with the gene names called \"gene\" in your genelist to proceed with this function.")
  }

  # will save genelist in the working directory
  message(paste0("HOMER genelist txt file will be saved in the working directory: ", getwd()))

  # load conversion table
  convert.table.Z11 <- read.delim("/Users/hunterm/Documents/R/from_Nate/GRCz11_to_HS.txt")
  # select just fish gene names and ensembl IDs
  ensembl.fish <- convert.table.Z11 %>% dplyr::select(Ensembl, Zebrafish_Symbol)
  colnames(ensembl.fish)[2] <- "gene"

  if (remove_ribosomal_genes) {
    genelist <- genelist %>% .[grep(pattern = "^rps|^rpl", x = .$gene, invert = T, ignore.case = T),]
  }

  homer_genelist <- genelist %>%
    filter(p_val_adj <= 0.05) %>%
    plyr::join(x = ensembl.fish,
               y = .,
               by = "gene",
               type = "inner") %>%
    arrange(-avg_logFC)

  file_name <- paste0(save_as, "_homer.txt")
  write.table(homer_genelist, file = file_name, sep = "\t", row.names = F, col.names = T, quote = F)

}
