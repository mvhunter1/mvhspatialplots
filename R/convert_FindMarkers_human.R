#' @title convert_FindMarkers_human
#' @description convert the results of FindMarkers (dataframe) to human genes and create ranked vector for GSEA. * Improved version.
#' @param genelist FindMarkers output dataframe. Needs to have a column with the gene names called "gene" to proceed with this function.
#' @param remove_ribosomal_genes default = T, remove ribosomal genes from FindMarkers results?
#' @param create_ranked_vector default = T, create ranked vector to use with fgsea?
#' @export
#' @return ranked vector or dataframe with converted human genes.

convert_FindMarkers_human <- function(genelist, remove_ribosomal_genes = T, create_ranked_vector = T) {

  require(data.table)
  require(tidyverse)

  fish.human.convert.Z11 <- read.delim("/Users/hunterm/Documents/R/from_Nate/GRCz11_to_HS.txt")
  fish.human.convert.Z11 <- fish.human.convert.Z11[fish.human.convert.Z11$DIOPT_Score > 6,]

  # check for presence of a "gene" column in df. If not, give an error message.
  if (!"gene" %in% colnames(genelist)) {
    stop("You must have a column with the gene names called \"gene\" in your FindMarkers results to proceed with this function.")
  }

  if (remove_ribosomal_genes) {
    genelist <- genelist[grep(pattern = "^rps|^rpl", x = genelist$gene, invert = T),]
  }

  # merge with conversion table
  merged <- merge(x = genelist,
                  y = fish.human.convert.Z11,
                  by.x = "gene",
                  by.y = "Zebrafish_Symbol")

  # check for duplicates
  dups <- merged$Human_Symbol[duplicated(merged$Human_Symbol)] %>% unique() %>% as.character()

  if (length(dups) >= 1) {
    message("Some duplicated genes found, removing...")
    # pull out duplicated genes
    dups_data <- NULL
    for (dup in dups) {
      data <- merged[merged$Human_Symbol == dup,] %>% arrange(-avg_logFC)
      dups_data[[dup]] <- data[1,] # take the expression data for the gene with the highest logFC
    }
    dups_data <- rbindlist(dups_data)

    # merge back with non duplicated genes
    merged_new <- rbind(dups_data,
                        merged[!merged$Human_Symbol %in% dups,])
    message("Duplicated genes removed, conversion complete!")
    if (create_ranked_vector) {
      vector <- merged_new %>% dplyr::select(Human_Symbol, avg_logFC) %>% deframe()
      return(vector)
    } else {
      return(merged_new)
    }
  } else {
    if (create_ranked_vector) {
      message("No duplicated genes found, returning converted list.")
      vector <- merged %>% dplyr::select(Human_Symbol, avg_logFC) %>% deframe()
      return(vector)
    } else {
      return(merged)
    }
  }
}
