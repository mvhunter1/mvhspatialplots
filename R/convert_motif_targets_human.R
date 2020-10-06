#' @title convert_motif_targets_human
#' @description convert HOMER motif target genes to human and add corresponding data - intended for cases when you want to perform GSEA on the HOMER motif target genes.
#' @param motif_genes character vector of the motif target genes (fish).
#' @param expression_data dataframe of the corresponding expression data for your genelist (usually output of FindMarkers). NOTE: the function assumes the gene names are in column 2 of this dataframe - if you're getting an error, check that first.
#' @param create_ranked_vector if T, will create a ranked vector of the human orthologs and expression data - to be used for GSEA.
#' @export
#' @return if create_ranked_vector = T, will return a ranked vector to be used for GSEA. if F, a dataframe.
#'

convert_motif_targets_human <- function(motif_genes, expression_data, create_ranked_vector = T) {
  # convert HOMER target genes to human

  motif.genes <- motif_genes
  expression.data <- expression_data

  # colnames(expression.data)[7] <- "Zebrafish_Symbol"

  ###### This assumes the gene names are in column 2! If getting an error, check that first
  colnames(expression.data)[2] <- "Zebrafish_Symbol"

  # load conversion table
  convert.table.Z11 <- read.delim("/Users/hunterm/Documents/R/from_Nate/GRCz11_to_HS.txt")
  fish.human.convert.Z11 <- convert.table.Z11[convert.table.Z11$DIOPT_Score > 6, ]

  # add expression data to gene lists
  motif.genes.use <- subset(expression.data, Zebrafish_Symbol %in% motif.genes) %>%
    dplyr::select(Zebrafish_Symbol, avg_logFC)

  # merge with conversion table
  motif.genes.convert <- merge(x = motif.genes.use,
                               y = fish.human.convert.Z11,
                               by = 'Zebrafish_Symbol')

  # extract duplicated human genes
  motif.gene.dups <- as.character(unique(motif.genes.convert$Human_Symbol[duplicated(motif.genes.convert$Human_Symbol)]))

  if (length(motif.gene.dups >= 1)) {
    message("Some human genes are duplicated, removing...")

    # initialize lists to store duplicates in
    motif.dups <- vector("list", length(motif.gene.dups))

    # remove duplicates
    for (ii in 1:length(motif.gene.dups)) {
      genes <- motif.genes.convert[motif.genes.convert$Human_Symbol == motif.gene.dups[ii], ] %>% arrange(desc(avg_logFC))
      motif.dups[[ii]] <- genes[1,]
    }

    motif.dups <- data.table::rbindlist(motif.dups)

    ## take all other non duplicate genes
    motif.nondups <- motif.genes.convert[grep(paste(pattern = motif.gene.dups,
                                                    collapse = "|"),
                                              x = motif.genes.convert$Human_Symbol,
                                              invert = TRUE), ]

    # make sure no genes are duplicated
    nondups.fish <- unique(duplicated(motif.nondups$Zebrafish_Symbol))
    nondups.human <- unique(duplicated(motif.nondups$Human_Symbol))
    dups.fish <- unique(duplicated(motif.dups$Zebrafish_Symbol))
    dups.human <- unique(duplicated(motif.dups$Human_Symbol))
    if (nondups.fish == FALSE & nondups.human == FALSE & dups.fish == FALSE & dups.human == FALSE) {
      message("No duplicated genes remain, ok to proceed")
    } else {
      stop("Some duplicated genes remain")
    }

    # merge to create final gene list, without duplicated fish genes
    motif.genes.human <- rbind(motif.dups, motif.nondups)

    # check for duplicates one more time
    nondups.final.fish <- unique(duplicated(motif.genes.human$Zebrafish_Symbol))
    nondups.final.human <- unique(duplicated(motif.genes.human$Human_Symbol))
    if (nondups.final.fish == FALSE & nondups.final.human == FALSE) {
      rm(nondups.final.fish, nondups.final.human, nondups.fish, nondups.human, dups.fish, dups.human, motif.dups, motif.nondups, genes, fish.human.convert.Z11, convert.table.Z11, motif.genes.use, motif.genes.convert, motif.gene.dups, ii)
      message("Fish to human gene conversion finished")

      if (create_ranked_vector == T) {
        motif.genes.human.vector <- motif.genes.human %>%
          dplyr::select(Human_Symbol, avg_logFC) %>%
          deframe()
        return(motif.genes.human.vector)

      } else {
        return(motif.genes.human)
      }

    } else {
      stop("Some duplicated genes remain")
    }
  } else {
    message("No duplicated genes, returning gene list")

    if (create_ranked_vector == T) {
      motif.genes.human.vector <- motif.genes.convert %>%
        dplyr::select(Human_Symbol, avg_logFC) %>%
        deframe()
      return(motif.genes.human.vector)

    } else {
      return(motif.genes.convert)
    }

  }
}
