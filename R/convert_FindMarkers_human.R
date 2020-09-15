#' @title convert_FindMarkers_human
#' @description convert the output of the Seurat function FindMarkers from fish genes to human.
#' @param fish_genelist FindMarkers results dataframe.
#' @param create_ranked_vector if T, will convert the output into a ranked vector with human genes as name and avg_logFC as values, to use for GSEA.
#' @export
#' @return if create_ranked_vector = T, a ranked vector with human genes as name and avg_logFC as values, to use for GSEA. if F, a dataframe.
#'
convert_FindMarkers_human <- function(fish_genelist, create_ranked_vector = T) {
  # convert output of FindMarkers to human genes
  # if create_ranked_vector = T: will return a ranked vector with human gene symbols and logFC, to be used for GSEA

  fish.human.convert.Z11 <- read.delim("/Users/hunterm/Documents/R/from_Nate/GRCz11_to_HS.txt")
  fish.human.convert.Z11 <- fish.human.convert.Z11[fish.human.convert.Z11$DIOPT_Score > 6, ]

  genelist <- fish_genelist

  if (colnames(genelist)[1] == "p_val") {
    genelist <- genelist %>% tibble::rownames_to_column(var = "fish_gene")
  } else {
    colnames(genelist)[1] <- "fish_gene"
  }

  merge <- merge(x = genelist,
                 y = fish.human.convert.Z11,
                 by.x = "fish_gene",
                 by.y = "Zebrafish_Symbol")

  h.genelist.c <- merge$Human_Symbol %>% as.character()
  h.genelist <- list(h.genelist.c)

  if (length(unique(duplicated(h.genelist.c))) == 1) {
    message("No duplicated genes, conversion complete")

    if (create_ranked_vector == T) {
      h.genelist.vector <- merge %>% dplyr::select(Human_Symbol, avg_logFC) %>% deframe()
      return(h.genelist.vector)
    } else {
      return(merge)
    }

  } else {
    message("Some genes are duplicated, removing...")

    gene.dups <- as.character(unique(merge$Human_Symbol[duplicated(merge$Human_Symbol)]))
    dups <- vector("list", length(gene.dups))

    for (ii in 1:length(gene.dups)) {
      genes <- merge[merge$Human_Symbol == gene.dups[ii], ] %>%
        dplyr::arrange(desc(avg_logFC))
      dups[[ii]] <- genes[1,]
    }

    dups <- data.table::rbindlist(dups)

    nondups <- merge[grep(paste(pattern = gene.dups,
                                collapse = "|"),
                          x = merge$Human_Symbol,
                          invert = TRUE), ]

    test1 <- unique(duplicated(nondups$Human_Symbol))
    test2 <- unique(duplicated(dups$fish_gene))

    if (test1 == FALSE & test2 == FALSE) {
      h.genelist <- rbind(dups, nondups)
      test3 <- unique(duplicated(h.genelist$fish_gene))
      test4 <- unique(duplicated(h.genelist$Human_Symbol))
      if (test3 == FALSE & test4 == FALSE) {
        message("Duplicated genes removed, conversion complete")

        if (create_ranked_vector == T) {
          h.genelist.vector <- h.genelist %>% dplyr::select(Human_Symbol, avg_logFC) %>% deframe()
          return(h.genelist.vector)
        } else {
          return(h.genelist)
        }

      } else {
        stop("Some duplicated genes weren't removed properly") }

    } else {
      stop("Some duplicated genes weren't removed properly")
    }
  }
}
