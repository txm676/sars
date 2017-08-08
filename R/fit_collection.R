######Fit Collection########

#' @export

fit_collection <- function(...){
  res <- list(...)
  class(res) <- "mmsar2_fitcollection"
  invisible(res)
}
