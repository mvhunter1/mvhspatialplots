#' @title calculate_Rsq
#' @description calculate Rsq value from linear regression.
#' @param dataframe dataframe with the values you're using for linear regression.
#' @param var_1 one of the variables you're using for linear regression inputted as string.
#' @param var_2 the other variable you're using for linear regression inputted as string.
#' @export
#' @return Rsq value.

calculate_Rsq <- function(dataframe, var_1, var_2) {
  if (typeof(dataframe) == "list") {
    dataframe <- dataframe %>% data.frame()
  }
  col1 <- which(colnames(dataframe) == var_1)
  col2 <- which(colnames(dataframe) == var_2)
  lm_res <- lm(dataframe[,col1] ~ dataframe[,col2])
  # summary.lm(lm_res)
  rsq <- summary(lm_res)$r.squared %>% round(digits = 4)
  return(rsq)
}
