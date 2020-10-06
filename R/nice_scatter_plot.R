#' @title nice_scatter_plot
#' @description create a nice looking scatter plot comparing 2 variables, with optional linear regression and calculation of R and Rsq.
#' @param df dataframe containing the data you want to plot.
#' @param x_var variable to be plotted on the x axis, as string.
#' @param y_var variable to be plotted on the y axis, as string.
#' @param fit_line if T, will fit a line using linear regression.
#' @param calculate_rsq if T, will calculate Rsq.
#' @param calculate_R if T, will calculate R.
#' @export
#' @return scatter plot. R and Rsq will be printed on plot if calculated.
#'

nice_scatter_plot <- function(df, x_var, y_var, fit_line = T, calculate_rsq = T, calculate_R = T) {

  require(tidyverse)

  plot <- ggplot(df, aes_string(x = x_var, y = y_var)) +
    geom_point() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 13))

  if (fit_line) {
    plot <- plot + geom_smooth(method = "lm", color = "#B8122690", size = 2)
  }

  if (calculate_rsq == T) {
    # source("/Users/hunterm/Dropbox/MH_ST/Miranda_R/functions/calculate_Rsq.R")
    rsq <- calculate_Rsq(df, x_var, y_var)
  }

  if (calculate_R) {
    cor <- cor(df %>% dplyr::select(x_var, y_var))
    r_val <- cor[1,2] %>% signif(digits = 3)
  }

  if (calculate_R & calculate_rsq) {
    plot <- plot + labs(subtitle = paste("Rsq", "=", rsq, "R", "=", r_val))
  } else if (calculate_R == T & calculate_rsq == F) {
    plot <- plot + labs(subtitle = paste("R", "=", r_val))
  } else if (calculate_R == F & calculate_rsq == T) {
    plot <- plot + labs(subtitle = bquote(R^2 ~ "=" ~ .(rsq)))
  }

  return(plot)

}
