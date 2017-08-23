######Fit Collection########

#' @export

fit_collection <- function(..., fits = list(...)){
  
  if (!all(sapply(fits[-1], FUN = function(x,y)identical(x$data,y$data), y = fits[[1]]))) stop("All models were not fitted on the same data set !")
  
  class(fits) <- "sars"
  attr(fits, "type") <- "fit_collection"
  invisible(fits)
}#eo fit_collection
