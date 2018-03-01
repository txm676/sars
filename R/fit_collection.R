######Fit Collection########

#' @export

fit_collection <- function(..., fits = list(...)){
  if (all(is.na(fits))) stop("all fits where NA :(, cannot create the 'fit_collection'")
  if (!all(sapply(fits[-1], FUN = function(x,y)base::identical(x$data,y$data), y = fits[[1]]))) stop("All models were not fitted on the same data set !")
  mod_names <- vapply(fits, FUN = function(x){base::gsub(pattern = " ", replacement = "_",  x = x$model$name)}, FUN.VALUE = character(1))
  names(fits) <- mod_names
  class(fits) <- "sars"
  attr(fits, "type") <- "fit_collection"
  return(fits)
}#eo fit_collection
