# POWER MODEL (ARRHENIUS 1921)

#' @export

sar_power <- function(data=galap, custstart = NULL, normtest = "lillie"){

  if (!(is.matrix(data) || is.data.frame(data))) stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop("NAs present in data")
  normtest <- match.arg(normtest, c("none", "shapiro", "kolmo", "lillie"))

  data <- data[order(data[,1]),]
  colnames(data) <- c("A", "S")

  model <- list(
    name = c("Power"),
    formula = expression(S == c * A ^ z),
    exp = expression(c * A ^ z),
    shape="convex",
    asymp=function(pars)FALSE,
    parLim = c("R", "R"),
    init = function(data){
      if (any(data$S == 0)){
        log.data = data.frame(A = log(data$A), S = log(data$S + .5))
      } else {
        log.data = log(data)
      }
      res = lm(S ~ A, log.data)$coefficients
      res = c(exp(res[1]), res[2])
      names(res) = c("c", "z")
      return(res)
    }
  )

  if (is.null(custstart)){
    model$custstart <- function(data) c(5, .25, 0)
  } else {
    model$custstart <- custstart
  }

  model <- compmod(model)

  fit <- rssoptim(model, data, custstart, normtest, algo = "Nelder-Mead")
  obs <- obs_shape(fit)
  fit$observed_shape <- obs$fitShape
  fit$asymptote <- obs$asymp

  class(fit) <- "sars"
  attr(fit, "type") <- "fit"
  return(fit)
}
