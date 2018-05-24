#' Fit the Heleg(Logistic) model

#' @description Fit the Heleg(Logistic) model to SAR data
#' @usage sar_heleg(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in inital optimization
#'   to run a grid search.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_heleg(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_heleg <- function(data = galap, start = NULL, grid_start = NULL, normaTest =  "lillie",
              homoTest = "cor.fitted"){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (base::anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#LOGISTIC FUNCTION (HE & LEGENDRE 1996)
model <- list(
  name=c("Heleg(Logistic)"),
  formula=expression(S == c/(f + A^(-z))),
  exp=expression(c/(f + A^(-z))),
  shape="sigmoid",
  asymp=function(pars)pars["c"]/pars["f"],
  parLim = c("Rplus","Rplus","Rplus"),
  custStart=function(data)c(max(data$S),10,.01),
  #initial values function
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    #c calculation (asymptote)
    c=max(data$S)+1
    #Intermediate variable calculation
    #long=length(data[[2]])
    Z=log((c/data$S) -1)
    #We have the z and f init values by linear regression of Z on data[[1]]
    dat=data.frame("a"=data$A,"Z"=Z)
    zf=stats::lm(Z~a,dat)$coefficients
    c(max(data$S)*exp(-zf[1]),exp(-zf[1]),-zf[2])
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
}#end of sar_heleg
