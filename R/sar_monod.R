#' Fit the Monod model

#' @description Fit the Monod model to SAR data
#' @usage sar_monod(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in optim
#'   to run a grid search.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_monod(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_monod <- function(data = galap, start = NULL, grid_start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
# MONOD CURVE (MONOD 1950, Willimas et al. 2009 formula)
model <- list(
  name=c("Monod"),
  formula=expression(S==d/(1+c*A^(-1))),
  exp=expression(d/(1+c*A^(-1))),
  shape="convex",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","Rplus"),
  custStart=function(data)c(stats::quantile(data$A,c(0.25)),max(data$S)),
  #initials values function
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    d=as.double(max(data$A)+max(data$S)/4)
    c=data[[1]]*(d/data$S - 1)
    c(d,stats::quantile(c,c(0.25)))
  }
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
}#end of sar_monod
