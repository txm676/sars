#' Create a Collection of SAR Model Fits

#' @description Creates a fit collection of SAR model fits, which can then be
#'   plotted using \code{\link{plot.sars}}.
#' @usage fit_collection(..., fits = list(...))
#' @param ... A set of one or more SAR model fits (all of class 'sars').
#' @param fits Internal argument that creates a list of the model fits.
#' @return A list of class 'sars' with n elements, corresponding to the
#'   n individual SAR model fits.
#' @examples
#' data(galap)
#' fit <- sar_linear(galap)
#' fit2 <- sar_power(galap)
#' fitC <- fit_collection(fit, fit2)
#' plot(fitC)
#' @export


fit_collection <- function(..., fits = list(...)){
  if (all(is.na(fits))) stop("all fits where NA :(, cannot create the 'fit_collection'")
  if (!all(vapply(fits[-1], FUN = function(x,y)base::identical(x$data,y$data), y = fits[[1]], FUN.VALUE = logical(1)))) stop("All models were not fitted on the same data set !")
  mod_names <- vapply(fits, FUN = function(x){base::gsub(pattern = " ", replacement = "_",  x = x$model$name)}, FUN.VALUE = character(1))
  names(fits) <- mod_names
  class(fits) <- "sars"
  attr(fits, "type") <- "fit_collection"
  return(fits)
}#eo fit_collection
