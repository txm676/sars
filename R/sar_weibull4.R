#' Fit the Cumulative Weibull 4 par. model

#' @description Fit the Cumulative Weibull 4 par. model to SAR data
#' @usage sar_weibull4(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in optim
#'   to run a grid search.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_weibull4(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_weibull4 <- function(data = galap, start = NULL, grid_start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#CUMULATIVE WEIBULL DISTRIBUTION with 4 parameters
model <- list(
  name=c("Cumulative Weibull 4 par."),
  formula=expression(S == d * (1 - exp(-c*A^z))^f ),
  exp=expression(d * (1 - exp(-c*A^z))^f ),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","Rplus","Rplus","Rplus"),
  custStart=function(data)c(10,.01,max(data$S),.01),
  #initial values function
  init=function(data){
    data = data[data$S!=0,]
    #c calculation (asymptote)
    d=max(data$S)+max(data$S)/4
    #z calculation
    Z=log(-log((d-data$S)/d)) # = log z + flogX
    Z[][Z == Inf]=NA
    c=exp(min(Z))
    dat=data.frame("A"=log(data[[1]]),"S"=Z)
    c=exp(stats::lm(S~A,dat)$coefficients[[1]])
    #f calculation
    z=stats::lm(S~A,dat)$coefficients[[2]]
    c(d,c,z,1)
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
}#end of sar_weibull4
