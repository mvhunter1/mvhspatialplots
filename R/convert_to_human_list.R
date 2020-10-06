#' @title convert_to_human_list
#' @description convert a fish genelist (character vector) to human, using DIOPT scores to remove duplicated genes. Intended to be used in cases where you don't have corresponding expression data for the genelist.
#' @param genelist list of fish genes as character vector.
#' @export
#' @return a character vector of the human orthologs of the given fish genes.
#'

convert_to_human_list <- function(genelist) {

  fish.human.convert.Z11 <- read.delim("/Users/hunterm/Documents/R/from_Nate/GRCz11_to_HS.txt")
  fish.human.convert.Z11 <- fish.human.convert.Z11[fish.human.convert.Z11$DIOPT_Score > 6, ]

  genelist <- genelist[!is.na(genelist)] %>% data.frame()

  if (nrow(genelist) < 1) { stop("Gene list to convert is empty") }

  colnames(genelist)[1] <- "fish_gene"

  merge <- merge(x = genelist,
                 y = fish.human.convert.Z11,
                 by.x = "fish_gene",
                 by.y = "Zebrafish_Symbol")

  h.genelist.c <- merge$Human_Symbol %>% as.character()
  h.genelist <- list(h.genelist.c)

  if (length(unique(duplicated(h.genelist.c))) == 1) {
    message("No duplicated genes, conversion complete")
    return(h.genelist.c)

    } else {
    message("Some genes are duplicated, removing...")

    gene.dups <- as.character(unique(merge$Human_Symbol[duplicated(merge$Human_Symbol)]))
    dups <- vector("list", length(gene.dups))

    for (ii in 1:length(gene.dups)) {
      genes <- merge[merge$Human_Symbol == gene.dups[ii], ] %>%
        arrange(desc(DIOPT_Score))
      dups[[ii]] <- genes[1,]
    }

    dups <- data.table::rbindlist(dups)

    nondups <- merge[grep(paste(pattern = gene.dups,
                                collapse = "|"),
                          x = merge$Human_Symbol,
                          invert = TRUE), ]

    test1 <- unique(duplicated(nondups$Human_Symbol))
    test2 <- unique(duplicated(dups$Human_Symbol))

    if (test1 == FALSE & test2 == FALSE) {
      h.genelist.c <- rbind(dups, nondups) %>% pull(Human_Symbol) %>% as.character()
      test3 <- unique(duplicated(h.genelist.c))

      if (test3 == FALSE) {
        h.genelist <- list(h.genelist.c)
        message("Duplicated genes removed, conversion complete")
        return(h.genelist.c)
      }
    } else {
      stop("Some duplicated genes weren't removed properly")
    }
  }
}
