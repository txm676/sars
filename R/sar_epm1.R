#' Fit the Asymptotic regression model

#' @description Fit the Asymptotic regression model to SAR data
#' @usage sar_epm1(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_epm1(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_epm1 <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
# EXTENDED POWER MODEL 1 (TJORVE 2009)
model <- list(
  name=c("Extended Power model 1"),
  formula=expression(s==c*a^(z*a^-d)),
  exp=expression(c*A^(z*A^-d)),
  shape="sigmoid",
  asymp=function(pars)FALSE,
  parLim = c("Rplus","R","R"),
  custStart=function(data)c(5,.25,.15),
  #initials values function
  init=function(data){if(any(data$S==0)){log.data=data.frame(S=log(data$a),S=log(data$S+.5))}else{log.data=log(data)};res=lm(S~A,log.data)$coefficients;res=c(exp(res[1]),res[2],.15);names(res)=model$parNames;return(res)}
)


model <- compmod(model) 
fit <- rssoptim(model = model, data = data, custstart = start, algo = 'Nelder-Mead') 
obs <- obs_shape(fit) 
fit$observed_shape <- obs$fitShape 
fit$asymptote <- obs$asymp 
class(fit) <- 'sars' 
attr(fit, 'type') <- 'fit' 
return(fit) 
}#end of sar_epm1
