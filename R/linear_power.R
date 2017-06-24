#' Fit the log-log version of the power model
#'
#' @description Fit the log-log version of the power model to SAR data and
#'   return parameter values, summary statistics and the fitted values. A check
#'   is made for any islands with zero species. If any zero species islands are
#'   found, a constant (default = 1) is added to each species richness value to
#'   enable log transformation.Natural logarithms are used.
#' @usage lin_pow(dat, a = 1, s = 2, con = 1)
#' @param dat A dataset in the form of list. The first element of the list is
#'   the name of the dataset (e.g. "galap"), and the second element is a
#'   dataframe with a minimum of two columns: one with island/site areas, and
#'   one with the species richness of each island/site.
#' @param a The column number of the area values. The default is 1 (i.e. column 1)
#' @param s The column number of the species richness values. The default is 2 (i.e. column 2)
#' @param con The constant to add to the species richness values in cases where one of the islands has zero species
#' @return A list of class 'mmSAR' with a 'Type' attribute
#'   'lin_pow'. The first element in the list is a vector of overlap values for
#'   the "in nodes" and the second element is a vector of overlap values for the
#'   "out nodes".
#'
#'   The \code{\link{summary.}} methods provides more useful summary
#'   statistics.
#' @examples
#' data(galap)
#' y <-  boreal[1:300,] #subset 300 rows for speed
#' d <- sample(nrow(y), 200, replace = FALSE) #create a random pot_net
#' pot_net <- y[d,] #by randomly sampling 200 rows from boreal
#' x <- NOSM_POT_dir(y, pot_net, perc = 1, sl = 1)
#' summary(x)
#' @export

getLinearPower <- function(dat, name = "mmsAR_dataset", a = 1, s = 2, con = 1) {

  dat <- data.frame(dat$data[, a],dat$data[, s])
  names(dat) = c("A", "S")

  dat.rssoptim <- list(nam = name, data = dat)

  tri <- order(dat.rssoptim$data$A)

  dat.rssoptim$data <- dat.rssoptim$data[tri,]

  dataRes <- dat.rssoptim$data

  logDat <- log(dataRes$A)

  if(any(dat.rssoptim$data$S == 0)){
    log.data = data.frame(A = log(dat.rssoptim$data$A), S = log(dat.rssoptim$data$S + con))
  } else {
    log.data = log(dat.rssoptim$data)
  }

  linearPower.fit = tryCatch(summary(lm(S ~ A, data = log.data)), error = function(e){NA})

  if (is.na(linearPower.fit)){
    error("ERROR in LINEAR POWER fit\n")
  }

  c = linearPower.fit$coefficients[1,1]
  z = linearPower.fit$coefficients[2,1]
  z.sig = linearPower.fit$coefficients[2,4]
  r2 = linearPower.fit$r.squared

  res <- c(c, z, z.sig, r2)
  names(res) <- c("c", "z", "z.sig", "r2")

  res <- c(res.dat,res)

  return(res)
}
