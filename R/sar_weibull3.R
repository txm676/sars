#' Fit the Cumulative Weibull 3 par. model

#' @description Fit the Cumulative Weibull 3 par. model to SAR data
#' @usage sar_weibull3(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_weibull3(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_weibull3 <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#CUMULATIVE WEIBULL DISTRIBUTION with 3 parameters
model <- list(
  name=c("Cumulative Weibull 3 par."),
  formula=expression(S == d(1 - exp(-c*A^z)) ),
  exp=expression(d*(1 - exp(-c*A^z)) ),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","Rplus","Rplus"),
  custStart=function(data)c(10,.01,max(data$S)),
  init=function(data){
    data = data[data$S!=0,]
    #c calculation (asymptote)
    d=max(data$S)+max(data$S)/4
    #z calculation
    Z=log(-log((d-data$S)/d)) # = log z + flogX
    Z[][Z == Inf]=NA
    c=exp(min(Z))
    dat=data.frame("A"=log(data$A),"S"=Z)
    c=exp(stats::lm(S~A,dat)$coefficients[[1]])
    #f calculation
    z=stats::lm(S~A,dat)$coefficients[[2]]
    c(d,c,z)
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
}#end of sar_weibull3
