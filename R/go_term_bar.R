#' @title go_term_bar
#' @description nice looking bar graph of top GO terms from GSEA by pval or NES.
#' @param gsea_results dataframe with GSEA results.
#' @param n_terms number of top GO terms to plot.
#' @param metric must be "NES" or "pval".
#' @param fill if T, will fill bars with gradient colours, if F, will fill bars with grey only.
#' @param change_labels to change pathway names to lowercase and remove "GO:" from the pathway names.
#' @param pval_fill if T, will colour bars based on pval.
#' @export
#' @return bar graph.

go_term_bar <- function(gsea_results, n_terms = 10, metric = "NES", fill = T, change_labels = T, pval_fill = T) {
  if (change_labels == T) {
    gsea_results$pathway <- gsea_results$pathway %>%
      gsub(pattern = "GO_", replacement = "", x = ., ignore.case = F) %>%
      gsub(pattern = "_", replacement = " ", x = .) %>%
      gsub(pattern = "plus", replacement = "+", x = .) %>%
      tolower()
  }
  if (fill == T) {
    cols_pval <- pals::viridis(n_terms)
    cols_NES <- pals::magma(n_terms)
  } else if (fill == F) {
    cols_pval <- "grey"
    cols_NES <- "grey"
  }
  if (metric == "pval") {
    data.bar <- gsea_results %>%
      dplyr::filter(NES > 0) %>%
      dplyr::arrange(-log10pval) %>%
      dplyr::slice(., 1:n_terms)
    if (pval_fill) {
      # updated to fix issue with ordering of Y axis
      data.bar$pathway <- factor(data.bar$pathway, levels = rev(data.bar$pathway))
      bar <- ggplot(data.bar, aes(x = log10pval, y = pathway, fill = log10pval)) +
        geom_bar(stat = "identity") +
        scale_fill_gradientn(colours = ocean.matter(n = 1000) %>% rev()) +
        theme_minimal() +
        labs(x = "-log10 p-value",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar)
    } else {
      data.bar$pathway <- factor(data.bar$pathway, levels = rev(data.bar$pathway))
      bar <- ggplot(data.bar, aes(x = log10pval, y = pathway)) +
        geom_bar(stat = "identity", fill = cols_pval) +
        theme_minimal() +
        labs(x = "-log10 p-value",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar)
    }
  } else if (metric == "NES") {
    data.bar <- gsea_results %>%
      dplyr::arrange(-NES) %>%
      dplyr::slice(., 1:n_terms)
    if (pval_fill) {
      data.bar$pathway <- factor(data.bar$pathway, levels = rev(data.bar$pathway))
      bar <- ggplot(data.bar, aes(x = NES, y = pathway, fill = pval)) +
        geom_bar(stat = "identity") +
        scale_fill_gradientn(colours = ocean.matter(n = 1000)) +
        theme_minimal() +
        labs(x = "normalized enrichment score",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar)
    } else {
      data.bar$pathway <- factor(data.bar$pathway, levels = rev(data.bar$pathway))
      bar <- ggplot(data.bar, aes(x = NES, y = pathway)) +
        geom_bar(stat = "identity", fill = cols_NES) +
        theme_minimal() +
        labs(x = "normalized enrichment score",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(bar)
    }
  } else {
    stop("metric must be NES or pval")
  }
}
