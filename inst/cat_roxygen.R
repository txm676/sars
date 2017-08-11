

cat_roxygen <- function(model, funName, fileName){
  
  cat1 <- function(...){cat(..., file = fileName, append = T)}
  
  cat1(paste0("#' ","Fit the ", model$name," model", "\n"))
  cat1("\n")
  cat1(paste0("#' @description ","Fit the ", model$name," model to SAR data", "\n"))
  cat1(paste0("#' @usage ", funName,"(data, custstart = NULL, normtest = "'lillie'")", "\n"))
  
  
  
  
  
  cat("#@'export \n")
  
  cat("\n")
  
  
}#eo cat_roxygen



#' Fit the log-log version of the power model
#'
#' @description Fit the log-log version of the power model to SAR data and
#'   return parameter values, summary statistics and the fitted values. A check
#'   is made for any islands with zero species. If any zero species islands are
#'   found, a constant (default = 1) is added to each species richness value to
#'   enable log transformation. Natural logarithms are used.
#' @usage lin_pow(dat, con = 1)
#' @param x A dataset in the form of a dataframe with two columns:
#'   the first with island/site areas, and the second with the species richness of each
#'   island/site.
#' @param con The constant to add to the species richness values in cases where
#'   one of the islands has zero species
#' @return A list of class "sars" with four elements. The first element is an
#'   object of class 'summary.lm'. This is the summary of the linear model fit
#'   using the \link[stats]{lm} function and the user's data. The second element
#'   is a numeric vector of the model's fitted values, and the third and fourth
#'   contain the island areas and observed richness values, respectively.
#'
#'   The \code{\link{summary.sars}} function returns a more useful summary of the
#'   model fit results, and the \code{\link{plot.}} plots the model.
#' @examples
#' data(galap)
#' fit <- lin_pow(galap, con = 1)
#' summary(fit)
#' plot(fit)
#' @export






