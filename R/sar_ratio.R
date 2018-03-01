#' Fit the Rational function model

#' @description Fit the Rational function model to SAR data
#' @usage sar_ratio(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in optim
#'   to run a grid search.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_ratio(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_ratio <- function(data = galap, start = NULL, grid_start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#RATIONAL FUNCITON ratkowski (1990)
model <- list(
  name=c("Rational function"),
  formula=expression(S == (c + z*A)/(1+d*A)),
  exp=expression((c + z*A)/(1+d*A)),
  shape="convex",
  asymp=function(pars)pars["z"]/pars["d"],
  parLim = c("R","Rplus","unif"),
  custStart=function(data)c(1,1,0.000001),
  #initial values function
  init=function(data){c(0,0,.5)}
)

model <- compmod(model) 
fit <- get_fit(model = model, data = data, start = start, grid_start = grid_start, algo = 'Nelder-Mead', verb = TRUE) 
if(is.na(fit$value)){ 
  return(NA) 
}else{ 
  obs <- obs_shape(fit) 
  fit$observed_shape <- obs$fitShape 
  fit$asymptote <- obs$asymp 
  class(fit) <- 'sars' 
  attr(fit, 'type') <- 'fit' 
  return(fit) 
} 
}#end of sar_ratio
