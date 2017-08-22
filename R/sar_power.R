#' Fit the Asymptotic regression model

#' @description Fit the Asymptotic regression model to SAR data
#' @usage sar_power(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_power(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_power <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
# POWER MODEL (ARRHENIUS 1921)
model <- list(
    name = c("Power"),
    formula = expression(S == c * A ^ z),
    exp = expression(c * A ^ z),
    shape="convex",
    asymp=function(pars)FALSE,
    parLim = c("R", "R"),
    custStart=function(data)c(5,0.25),
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
model <- compmod(model) 
fit <- rssoptim(model = model, data = data, custstart = start, normtest = 'lillie', algo = 'Nelder-Mead') 
obs <- obs_shape(fit) 
fit$observed_shape <- obs$fitShape 
fit$asymptote <- obs$asymp 
class(fit) <- 'sars' 
attr(fit, 'type') <- 'fit' 
return(fit) 
}#end of sar_power
