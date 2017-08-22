#' Fit the Asymptotic regression model

#' @description Fit the Asymptotic regression model to SAR data
#' @usage sar_ratio(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_ratio(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_ratio <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#RATIONAL FUNCITON ratkowski (1990)
model <- list(
  name=c("Rational function"),
  formula=expression(S == over( (c + z*A) , (1+d*A) ) ),
  exp=expression((c + z*A)/(1+d*A)),
  shape="convex",
  asymp=function(pars)pars["z"]/pars["d"],
  parLim = c("R","Rplus","unif"),
  custStart=function(data)c(1,1,0.000001),
  #initial values function
  init=function(data){c(0,0,.5)}
)

model <- compmod(model) 
fit <- rssoptim(model = model, data = data, custstart = start, normtest = 'lillie', algo = 'Nelder-Mead') 
obs <- obs_shape(fit) 
fit$observed_shape <- obs$fitShape 
fit$asymptote <- obs$asymp 
class(fit) <- 'sars' 
attr(fit, 'type') <- 'fit' 
return(fit) 
}#end of sar_ratio
