#' @title gsea_waterfall_plot
#' @description double waterfall plot of the chosen top and bottom pathways by NES from GSEA. Optional labelling of pathways of interest.
#' @param gsea_results your GSEA results dataframe, should have one column called "pathway" and another called "NES", additional columns are fine too.
#' @param pathways_to_highlight a character vector with the exact names of the pathways you want highlighted on your waterfall plot in red circles.
#' @param labels_for_plot optional: a character vector with the exact names of the pathways you want labelled as text on your plot.
#' @param n_terms the number of terms to be plotted.
#' @export
#' @return double waterfall plot.

gsea_waterfall_plot <- function(gsea_results, pathways_to_highlight, labels_for_plot = NULL, n_terms = 250) {

  # find the top pathways by NES
  gsea.up <- gsea_results %>%
    dplyr::arrange(desc(NES)) %>%
    head(., n_terms) %>%
    dplyr::select(pathway, NES) %>%
    tibble::rownames_to_column(var = "index")

  # find the bottom pathways by NES
  gsea.down <- gsea_results %>%
    dplyr::arrange(desc(NES)) %>%
    tail(., n_terms) %>%
    dplyr::select(pathway, NES) %>%
    tibble::rownames_to_column(var = "index")

  # organize data frame for plotting
  gsea.up$index = 1:n_terms
  gsea.down$index = n_terms:1
  gsea.waterfall <- rbind(gsea.up, gsea.down)
  label_pathways <- intersect(pathways_to_highlight, gsea.waterfall$pathway)

  if (length(label_pathways) < 1) {
    stop("pathways_to_highlight are not in the top or bottom terms by NES")
  }

  labelled_pathways <- data.frame(pathway = label_pathways, label_pathway = "yes")
  unlabelled_pathways <- data.frame(pathway = gsea.waterfall$pathway[!gsea.waterfall$pathway %in% label_pathways],
                                    label_pathway = "no")
  pathways_all <- rbind(labelled_pathways, unlabelled_pathways)
  gsea.waterfall <- merge(gsea.waterfall,
                          pathways_all,
                          by = "pathway")

  gsea.waterfall$label_pathway <- factor(gsea.waterfall$label_pathway, levels = c("yes", "no"))

  # create plot
  plot <- ggplot(gsea.waterfall, aes(x = index, y = NES)) +
    geom_point(aes(color = label_pathway), size = 4) +
    scale_color_manual(values = c("#E01111", "#B2B2B220")) +
    geom_hline(yintercept = 0, size = 1) +
    theme_minimal() +
    # scale_y_continuous(limits = c(-3.1, 3.1), breaks = seq(-3,3,1)) + # change this to change the y axis limits
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.y = element_text(size = 16, color = "black"),
          panel.grid.major.x = element_blank()) +
    xlim(0, n_terms + 1) +
    ylab("normalized enrichment score")

  if (is.null(labels_for_plot)) {
    return(plot)
  } else {
    # optional: add text labels to the plot highlighting certain pathways
    plot <- plot + ggrepel::geom_text_repel(data = subset(gsea.waterfall, pathway %in% labels_for_plot),
                                            aes(label = pathway),
                                            box.padding = unit(1, "lines"),
                                            point.padding = unit(1, "lines"))
    return(plot)
  }
}
