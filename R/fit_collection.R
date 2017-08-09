######Fit Collection########

#' @export

fit_collection <- function(..., fits = list(...)){
  #res <- list(...)
  class(fits) <- "sars"
  attr(fits, "type") <- "fit_collection"
  invisible(fits)
}
