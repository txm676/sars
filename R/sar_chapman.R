#' Fit the Chapman Richards model

#' @description Fit the Chapman Richards model to SAR data
#' @usage sar_chapman(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in inital optimization
#'   to run a grid search.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_chapman(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_chapman <- function(data = galap, start = NULL, grid_start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (base::anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#Chapman–Richards 3 S = a [1 − exp(−bA)]c Flather (1996)
model <- list(
  name=c("Chapman Richards"),
  exp=expression(d * (1 - exp(-z*A)^c )),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","R","R"),
  #initials values function
  init=function(data){d=max(data$S); Z=( -log( (-data$S/(max(data$S)+1))+1))/data$A; z = mean(Z); c(d,z,1)}
)

model <- compmod(model) 
fit <- get_fit(model = model, data = data, start = start, grid_start = grid_start, algo = 'Nelder-Mead', verb = TRUE) 
if(is.na(fit$value)){ 
  return(list(value = NA)) 
}else{ 
  obs <- obs_shape(fit) 
  fit$observed_shape <- obs$fitShape 
  fit$asymptote <- obs$asymp 
  class(fit) <- 'sars' 
  attr(fit, 'type') <- 'fit' 
  return(fit) 
} 
}#end of sar_chapman
