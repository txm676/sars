######Fit Collection########

#' @export

fit_collection <- function(...){
  res <- list(...)
  class(res) <- "sars"
  attr(res, "type") <- "fitcollection"
  invisible(res)
}
