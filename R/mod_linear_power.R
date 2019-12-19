#' Fit the log-log version of the power model
#'
#' @description Fit the log-log version of the power model to SAR data and
#'   return parameter values, summary statistics and the fitted values.
#' @usage lin_pow(data, con = 1, logT = log, compare = FALSE, normaTest =
#'   "lillie", homoTest = "cor.fitted")
#' @param data A dataset in the form of a dataframe with two columns: the first
#'   with island/site areas, and the second with the species richness of each
#'   island/site.
#' @param con The constant to add to the species richness values in cases where
#'   one of the islands has zero species.
#' @param logT The log-transformation to apply to the area and richness values.
#'   Can be any of \code{log}(default), \code{log2} or \code{log10}.
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
#'   species islands are found, a constant (default: \code{con = 1}) is added
#'   to each species richness value to enable log transformation. Natural
#'   logarithms are used as default, but log2 and log10 can be used instead
#'   using the \code{logT} argument.
#'
#'   The \code{compare} argument can be used to compare the c and z values
#'   calculated using the log-log power model with that calculated using the
#'   non-linear power model. Note that the log-log function returns log(c).
#' @importFrom stats lm summary.lm shapiro.test ks.test cor.test
#' @importFrom nortest lillie.test
#' @importFrom utils capture.output
#' @return A list of class "sars" with up to seven elements. The first element
#'   is an object of class 'summary.lm'. This is the summary of the linear model
#'   fit using the \link[stats]{lm} function and the user's data. The second
#'   element is a numeric vector of the model's fitted values, and the third
#'   contains the log-transformed observed data. The remaining elements depend
#'   on the function arguments selected and can include the results of the
#'   non-linear power model fit, the log-transformation function used (i.e.
#'   \code{logT}) and the results of the residuals normality and heterogeneity
#'   tests.
#'
#'   The \code{\link{summary.sars}} function returns a more useful summary of
#'   the model fit results, and the \code{\link{plot.sars}} plots the model.
#' @examples
#' data(galap)
#' fit <- lin_pow(galap, con = 1)
#' summary(fit)
#' plot(fit)
#' @export


lin_pow <- function(data, con = 1, logT = log,
                    compare = FALSE, normaTest =  "lillie", 
                    homoTest = "cor.fitted") {

  if (!(is.matrix(data) | is.data.frame(data)))
    stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop("NAs present in data")
  if (!is.primitive(logT)) stop("logT should be a (primitive) function,
                                specifically: log, log2 or log10")

  data <- data[order(data[,1]),]
  colnames(data) <- c("A", "S")
  
  xr <- range(data$S)/ mean(data$S)
  if (isTRUE(all.equal(xr[1], xr[2]))) {
    if (data$S[1] == 0){
      if (compare){
        warning("All richness values are zero: parameter estimates of",
               " non-linear models should be interpreted with caution")
      } else{
        warning("All richness values identical")
      }
    } else{
      warning("All richness values identical")
    }
  }
  
  if (any(data$S == 0)){
      log.data <- data.frame(A = logT(data$A), S = logT(data$S + con))
      linearPower.fit <- lm(logT(S + con) ~ logT(A), data = data)
      } else {
       log.data <- data.frame(A = logT(data$A), S = logT(data$S))
       linearPower.fit <- lm(logT(S) ~ logT(A), data = data)
  }

  fv <- linearPower.fit$fitted.values
  linearPower.fit <- summary.lm(linearPower.fit)
  resid <- linearPower.fit$residuals
  res <- list(Model = linearPower.fit, calculated = fv, data = log.data)

  if (compare == TRUE){
    pow <- sar_power(data)
    res$power <- pow
  }

  #normality and homogeneity tests
  normaTest <- match.arg(normaTest, c("none", "shapiro", "kolmo", "lillie"))
  homoTest <- match.arg(homoTest, c("none","cor.area","cor.fitted"))
  #normality of residuals
  if (normaTest == "shapiro") {
    normaTest <- list("test" = "shapiro", tryCatch(shapiro.test(resid),
                                                   error = function(e)NA))
  } else if (normaTest == "lillie"){
    normaTest <- list("test" = "lillie", tryCatch(lillie.test(resid),
                                                  error = function(e)NA))
  } else if (normaTest == "kolmo"){
    normaTest <- list("test" = "kolmo", tryCatch(ks.test(resid, "pnorm"),
                                                 error = function(e)NA))
  } else{
    normaTest <- "none"
  }
  #Homogeneity of variance
  if (homoTest == "cor.area"){
    homoTest  <- list("test" = "cor.area", tryCatch(cor.test(resid,data$A),
                            error = function(e)list(estimate=NA,p.value=NA)))
  } else if (homoTest == "cor.fitted"){
    homoTest  <- list("test" = "cor.fitted",
                  tryCatch(cor.test(resid,as.vector(res$calculated)),
                  error = function(e)list(estimate=NA,p.value=NA)))
  } else {
    homoTest <- "none"
  }
  res$normaTest <- normaTest
  res$homoTest <- homoTest
  lt <- capture.output(logT)
  res$logT <- switch(lt,
               "function (x, base = exp(1))  .Primitive(\"log\")" = "log()",
               "function (x)  .Primitive(\"log10\")" = "log10()",
               "function (x)  .Primitive(\"log2\")" = "log2()")
  class(res) <- "sars"
  attr(res, "type") <- "lin_pow"
  return(res)
}
