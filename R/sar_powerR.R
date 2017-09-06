#' Fit the PowerR model

#' @description Fit the PowerR model to SAR data
#' @usage sar_powerR(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_powerR(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_powerR <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
# POWER MODEL BIS (ROSENSWEIG 1995)
model <- list(
  name = c("PowerR"),
  formula = expression(S == f + c*A^z),
  exp = expression(f + c*A^z),
  shape = "convex",
  asymp = function(pars)FALSE,
  custStart = function(data)c(5,.25,0),
  #limits for parameters
  parLim = c("R","R","R"),
  #initials values function
  init = function(data) {
    if (any(data$S == 0)) {
      log.data <- data.frame(A = log(data$A), S = log(data$S + .5))
    }else{
      log.data <- log(data)
    }
    res <- stats::lm(S~A,log.data)$coefficients
    res <- c(exp(res[1]), res[2])
    res <- c(stats::lm(S~A,data)$coefficients[1], res)
    names(res) <- c("f","c","z")
    return(res)
    }
)

model <- compmod(model) 
fit <- rssoptim(model = model, data = data, start = start, algo = 'Nelder-Mead') 
obs <- obs_shape(fit) 
fit$observed_shape <- obs$fitShape 
fit$asymptote <- obs$asymp 
class(fit) <- 'sars' 
attr(fit, 'type') <- 'fit' 
return(fit) 
}#end of sar_powerR
