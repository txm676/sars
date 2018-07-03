#' Fit the log-log version of the power model
#'
#' @description Fit the log-log version of the power model to SAR data and
#'   return parameter values, summary statistics and the fitted values.
#' @usage lin_pow(dat, con = 1, compare = FALSE, normaTest =  "lillie", homoTest = "cor.fitted")
#' @param dat A dataset in the form of a dataframe with two columns: the first
#'   with island/site areas, and the second with the species richness of each
#'   island/site.
#' @param con The constant to add to the species richness values in cases where
#'   one of the islands has zero species.
#' @param compare Fit the standard (non-linear) power model and return the
#'   z-value for comparison (default: \code{compare = FALSE}).
#' @param normaTest The test used to test the normality of the residuals of the
#'   model. Can be any of "lillie" (Lilliefors Kolmogorov-Smirnov test; the
#'   default), "shapiro" (Shapiro-Wilk test of normality), "kolmo"
#'   (Kolmogorov-Smirnov test), or "none" (no residuals normality test is
#'   undertaken).
#' @param homoTest The test used to check for homogeneity of the residuals of
#'   the model. Can be any of "cor.fitted" (a correlation of the residuals with
#'   the model fitted values; the default), "cor.area" (a correlation of the
#'   residuals with the area values), or "none" (no residuals homogeneity test
#'   is undertaken).
#' @details A check is made for any islands with zero species. If any zero
#'   species islands are found, a constant (default: \code{con = 1}) is added to each
#'   species richness value to enable log transformation. Natural logarithms are
#'   used.
#'   
#'   The \code{compare} argument can be used to compare the c and z values
#'   calculated using the log-log power model with that calculated using the
#'   non-linear power model. Note that the log-log function returns logc.
#' @return A list of class "sars" with up to six elements. The first element is
#'   an object of class 'summary.lm'. This is the summary of the linear model
#'   fit using the \link[stats]{lm} function and the user's data. The second
#'   element is a numeric vector of the model's fitted values, and the third
#'   contains the log-transformed observed data. The remaining elements depend
#'   on the function arguments selected and can include the results of the
#'   non-linear power model fit, and of the residuals normality and
#'   heterogeneity tests.
#'
#'   The \code{\link{summary.sars}} function returns a more useful summary of
#'   the model fit results, and the \code{\link{plot.sars}} plots the model.
#' @examples
#' data(galap)
#' fit <- lin_pow(galap, con = 1)
#' summary(fit)
#' plot(fit)
#' @export


lin_pow <- function(data, con = 1, compare = F, normaTest =  "lillie", homoTest = "cor.fitted") {

  if (!(is.matrix(data) || is.data.frame(data))) stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop("NAs present in data")
  
  data <- data[order(data[,1]),]
  colnames(data) <- c("A", "S")

  if (any(data$S == 0)){
      log.data = data.frame(A = log(data$A), S = log(data$S + con))
  } else {
       log.data = data.frame(A = log(data$A), S = log(data$S))
  }
  linearPower.fit = stats::lm(S ~ A, data = log.data)

  fv <- linearPower.fit$fitted.values
  linearPower.fit <- stats::summary.lm(linearPower.fit)
  resid <- linearPower.fit$residuals
  res <- list(Model = linearPower.fit, calculated = fv, data = log.data)
  
  if (compare == T){
    pow <- sar_power(data)
    res$power <- pow
  }
  
  #normality and homogeneity tests
  normaTest <- match.arg(normaTest, c("none", "shapiro", "kolmo", "lillie"))
  homoTest <- match.arg(homoTest, c("none","cor.area","cor.fitted"))
  #normality of residuals
  if (normaTest == "shapiro") {
    normaTest <- list("test" = "shapiro", tryCatch(shapiro.test(resid), error = function(e)NA))
  } else if (normaTest == "lillie"){ 
    normaTest <- list("test" = "lillie", tryCatch(nortest::lillie.test(resid), error = function(e)NA))
  } else if (normaTest == "kolmo"){ 
    normaTest <- list("test" = "kolmo", tryCatch(ks.test(resid, "pnorm"), error = function(e)NA))
  } else{
    normaTest <- "none"
  }
  #Homogeneity of variance
  if (homoTest == "cor.area"){
    homoTest  <- list("test" = "cor.area", tryCatch(cor.test(resid,data$A), error = function(e)list(estimate=NA,p.value=NA)))
  } else if (homoTest == "cor.fitted"){
    homoTest  <- list("test" = "cor.fitted", tryCatch(cor.test(resid,as.vector(res$calculated)), 
                                                      error = function(e)list(estimate=NA,p.value=NA)))
  } else {
    homoTest = "none"
  }
  res$normaTest <- normaTest
  res$homoTest <- homoTest
  class(res) <- "sars"
  attr(res, "type") <- "lin_pow"
  return(res)
}







