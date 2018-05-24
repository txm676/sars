#' Fit the Asymptotic regression model

#' @description Fit the Asymptotic regression model to SAR data
#' @usage sar_asymp(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in inital optimization
#'   to run a grid search.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_asymp(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_asymp <- function(data = galap, start = NULL, grid_start = NULL, normaTest =  "lillie",
              homoTest = "cor.fitted"){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (base::anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#Asymptotic Regression
model <- list(
  name=c("Asymptotic regression"),
  exp=expression(d - c*z^A),
  shape="convex",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R"),
  #initial values function
  init=function(data){#Ratkowsky 1983 p178
    #d determination (asymptote)
    d=max(data$S)+max(data$S)/4
    #Intermediate variable calculation
    Z=log(d-data$S)
    #we have also Z=log(c)+Xlog(z) -> linear regression
    dat=data.frame("a"=data$A,"Z"=Z)
    zf=stats::lm(Z~a,dat)$coefficients
    c(d,exp(zf[1]),exp(zf[2]))
  }
)

model <- compmod(model) 
fit <- get_fit(model = model, data = data, start = start, grid_start = grid_start, algo = 'Nelder-Mead', 
       normaTest =  normaTest, homoTest = homoTest, verb = TRUE) 
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
}#end of sar_asymp
