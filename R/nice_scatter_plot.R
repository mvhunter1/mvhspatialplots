#' @title nice_scatter_plot
#' @description create nice looking scatter plot with optional linear regression/line fitting.
#' @param df dataframe with variables to plot.
#' @param x_var variable to be plotted on X axis inputted as string, should be a column in df.
#' @param y_var variable to be plotted on Y axis inputted as string, should be a column in df.
#' @param fit_line set to T to perform linear regression and fit a line.
#' @param calculate_rsq set to T to calculate Rsq of line fitting.
#' @export
#' @return scatter plot.

nice_scatter_plot <- function(df, x_var, y_var, fit_line = T, calculate_rsq = T) {
  plot <- ggplot(df, aes_string(x = x_var, y = y_var)) +
    geom_point() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 13))
  if (fit_line) {
    plot <- plot + geom_smooth(method = "lm", color = "#B8122690", size = 2)
  }
  if (calculate_rsq) {
    rsq <- mvhspatialplots::calculate_Rsq(df, x_var, y_var)
    plot <- plot + labs(subtitle = bquote(R^2 ~ "=" ~ .(rsq)))
  }
  return(plot)
}
