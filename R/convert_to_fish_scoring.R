#' @title convert_to_fish_scoring
#' @description convert a human genelist (character vector) to fish, using DIOPT scores to remove duplicated genes. Intended to be used in cases where you don't have corresponding expression data for the genelist.
#' @param genelist list of human genes as character vector.
#' @export
#' @return a character vector of the fish orthologs of the given fish genes.
#'

convert_to_fish_scoring <- function(genelist) {
  # this will return the converted genelist as a character vector!!

  fish.human.convert.Z11 <- read.delim("/Users/hunterm/Documents/R/from_Nate/GRCz11_to_HS.txt")
  fish.human.convert.Z11 <- fish.human.convert.Z11[fish.human.convert.Z11$DIOPT_Score > 6, ]

  genelist <- genelist[!is.na(genelist)] %>% data.frame()

  if (nrow(genelist) < 1) { stop("Gene list to convert is empty") }

  colnames(genelist)[1] <- "human_gene"

  merge <- merge(x = genelist,
                 y = fish.human.convert.Z11,
                 by.x = "human_gene",
                 by.y = "Human_Symbol")

  f.genelist.c <- merge$Zebrafish_Symbol %>% as.character()
  f.genelist <- list(f.genelist.c)

  if (length(unique(duplicated(f.genelist.c))) == 1) {
    message("No duplicated genes, conversion complete")
    return(f.genelist.c)
  } else {
    message("Some genes are duplicated, removing...")
    gene.dups <- as.character(unique(merge$Zebrafish_Symbol[duplicated(merge$Zebrafish_Symbol)]))
    dups <- vector("list", length(gene.dups))
    for (ii in 1:length(gene.dups)) {
      genes <- merge[merge$Zebrafish_Symbol == gene.dups[ii], ] %>%
        arrange(desc(DIOPT_Score))
      dups[[ii]] <- genes[1,]
    }
    dups <- data.table::rbindlist(dups)
    nondups <- merge[grep(paste(pattern = gene.dups,
                                collapse = "|"),
                          x = merge$Zebrafish_Symbol,
                          invert = TRUE), ]

    test1 <- unique(duplicated(nondups$Zebrafish_Symbol))
    test2 <- unique(duplicated(dups$Zebrafish_Symbol))
    if (test1 == FALSE & test2 == FALSE) {
      f.genelist.c <- rbind(dups, nondups) %>% pull(Zebrafish_Symbol) %>% as.character()
      test3 <- unique(duplicated(f.genelist.c))
      if (test3 == FALSE) {
        f.genelist <- list(f.genelist.c)
        message("Duplicated genes removed, conversion complete")
        return(f.genelist.c)
      }
    } else {
      stop("Some duplicated genes weren't removed properly")
    }
  }
}
