######Fit Collection########

#' @export

fit_collection <- function(...){
  res <- list(...)
  class(res) <- "mmsar2"
  attr(res, "type") <- "fitcollection"
  invisible(res)
}
