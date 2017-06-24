#' Fit the log-log version of the power model
#'
#' @description Fit the log-log version of the power model to SAR data and
#'   return parameter values, summary statistics and the fitted values. A check
#'   is made for any islands with zero species. If any zero species islands are
#'   found, a constant (default = 1) is added to each species richness value to
#'   enable log transformation. Natural logarithms are used.
#' @usage lin_pow(dat, a = 1, s = 2, con = 1)
#' @param x A dataset in the form of list. The first element of the list is
#'   the name of the dataset (e.g. "galap"), and the second element is a
#'   dataframe with a minimum of two columns: one with island/site areas, and
#'   one with the species richness of each island/site.
#' @param a The column number of the area values. The default is 1 (i.e. column 1)
#' @param s The column number of the species richness values. The default is 2 (i.e. column 2)
#' @param con The constant to add to the species richness values in cases where one of the islands has zero species
#' @return A list with three elements. The first element is an object of class
#'   'summary.lm'. This is the summary of the linear model fit using the
#'   \link[stats]{lm} function and the user's data. The second element is a
#'   numeric vector of the model's fitted values, and the third contains the island
#'   areas.
#'
#'   The \code{\link{summary.}} method provides a more useful summary. The \code{\link{plot.}} plots the model.
#' @examples
#' data(galap)
#' fit <- lin_pow(galap, a = 1, s = 2, con = 1)
#' @export
#' @importFrom stats lm

lin_pow <- function(x, a = 1, s = 2, con = 1) {

  dat <- data.frame(x$data[, a],x$data[, s])
  names(dat) = c("A", "S")
  dat.rssoptim <- list(Dataset = x$Dataset, data = dat)
  tri <- order(dat.rssoptim$data$A)
  dat.rssoptim$data <- dat.rssoptim$data[tri,]
  dataRes <- dat.rssoptim$data
  logDat <- log(dataRes$A)
  if (any(dat.rssoptim$data$S == 0)){
      log.data = data.frame(A = log(dat.rssoptim$data$A), S = log(dat.rssoptim$data$S + con))
  } else {
       log.data = log(dat.rssoptim$data)
  }
  linearPower.fit = tryCatch((lm(S ~ A, data = log.data)), error = function(e){NA})
  if (any(is.na(linearPower.fit))){
    error("ERROR in LINEAR POWER fit\n")
  }
  fv <- linearPower.fit$fitted.values
  linearPower.fit <- summary(linearPower.fit)
  res <- list(Model = linearPower.fit, Fitted = fv, Area = dat[,1])
  class(res) <- "mmSAR2"
  attr(res, "Type") <- "lin_pow"
  attr(res, "Dataset") <- x$Dataset
  return(res)
}







