#' Fit the Asymptotic regression model

#' @description Fit the Asymptotic regression model to SAR data
#' @usage sar_asymp(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_asymp(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_asymp <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#Asymptotic Regression
model <- list(
  name=c("Asymptotic regression"),
  formula=expression(S == d - c*z^A),
  exp=expression(d - c*z^A),
  shape="convex",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R"),
  custStart=function(data)c(1,1,max(data$S)*2),
  #initial values function
  init=function(data){#Ratkowsky 1983 p178
    #d determination (asymptote)
    d=max(data$S)+max(data$S)/4
    #Intermediate variable calculation
    Z=log(d-data$S)
    #we have also Z=log(c)+Xlog(z) -> linear regression
    dat=data.frame("a"=data$A,"Z"=Z)
    zf=lm(Z~a,dat)$coefficients
    c(d,exp(zf[1]),exp(zf[2]))
  }
)

model <- compmod(model) 
fit <- rssoptim(model = model, data = data, custstart = start, algo = 'Nelder-Mead') 
obs <- obs_shape(fit) 
fit$observed_shape <- obs$fitShape 
fit$asymptote <- obs$asymp 
class(fit) <- 'sars' 
attr(fit, 'type') <- 'fit' 
return(fit) 
}#end of sar_asymp
