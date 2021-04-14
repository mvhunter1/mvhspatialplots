#' @title nice_scatter_plot_w_outliers
#' @description nicer looking scatter plot with outlier points labelled.
#' @param df dataframe with data to plot.
#' @param x_var name of X variable column in df as string.
#' @param y_var name of Y variable column in df as string.
#' @param variable_to_label column in df with the labels you want on the plot for outliers.
#' @param low_outlier_cutoff bottom residual cutoff to select outliers, default -0.5
#' @param high_outlier_cutoff top residual cutoff to select outliers, default 0.5
#' @param fit_line plot a fitted line on the plot?
#' @param calculate_rsq calculate Rsq?
#' @param calculate_R calculate R?
#' @param label_outliers label outliers on plot?
#' @export
#' @return scatter plot.
#'


nice_scatter_plot_w_outliers <- function(df, x_var, y_var, variable_to_label, low_outlier_cutoff = -0.5, high_outlier_cutoff = 0.5, fit_line = F, calculate_rsq = T, calculate_R = T, label_outliers = T) {

  require(ggrepel)
  require(tidyverse)

  if (typeof(df) == "list") {
    df <- df %>% data.frame()
  }

  # this function performs a linear regression by default.
  col1 <- which(colnames(df) == x_var)
  col2 <- which(colnames(df) == y_var)
  lm_res <- lm(df[, col1] ~ df[, col2])

  plot <- ggplot(df, aes_string(x = x_var, y = y_var)) +
    geom_point() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 13))

  if (fit_line) {
    plot <- plot + geom_smooth(method = "lm", color = "#B8122690",
                               size = 2)
  }

  if (calculate_rsq) {
    rsq <- summary(lm_res)$r.squared %>% round(digits = 4)
  }

  if (calculate_R) {
    cor <- cor(df %>% dplyr::select(x_var, y_var))
    r_val <- cor[1, 2] %>% signif(digits = 3)
  }

  if (calculate_R & calculate_rsq) {
    plot <- plot + labs(subtitle = paste("Rsq", "=", rsq,
                                         "R", "=", r_val))
  }
  else if (calculate_R == T & calculate_rsq == F) {
    plot <- plot + labs(subtitle = paste("R", "=", r_val))
  }
  else if (calculate_R == F & calculate_rsq == T) {
    plot <- plot + labs(subtitle = bquote(R^2 ~ "=" ~ .(rsq)))
  }

  if (label_outliers) {
    df$residuals <- residuals(lm_res)
    outliers <- df[df$residuals >= high_outlier_cutoff | df$residuals <= low_outlier_cutoff, ] %>% pull(variable_to_label)

    plot <- plot + geom_text_repel(data = (df[df$gene %in% outliers,]),
                                  aes_string(x = x_var, y = y_var, label = variable_to_label, size = 14),
                                  segment.alpha = 0.5,
                                  show.legend = F,
                                  point.padding = 0.25,
                                  min.segment.length = 0.05)

  }
  return(plot)
}
