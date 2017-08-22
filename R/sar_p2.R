#' Fit the Persistence function 2 model

#' @description Fit the Persistence function 2 model to SAR data
#' @usage sar_p2(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_p2(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_p2 <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
# PERSISTENCE FUNCTION 2 (TJORVE 2009)
model <- list(
  name=c("Persistence function 2"),
  formula=expression(S == c*A^z * exp(-d/A)),
  exp=expression(c*A^z * exp(-d/A)),
  mod=s~c*a^z * exp(-d/a),
  shape="sigmoid",
  asymp=function(pars)FALSE,
  custStart=function(data)c(5,.25,.15),
  #limits for parameters
  parLim = c("Rplus","Rplus","R"),
  #initials values function
  init=function(data){if(any(data$S==0)){log.data=data.frame(A=log(data$A),S=log(data$S+.5))}else{log.data=log(data)};res=lm(S~A,log.data)$coefficients;c(res[1],res[2],.1)}
)

model <- compmod(model) 
fit <- rssoptim(model = model, data = data, custstart = start, algo = 'Nelder-Mead') 
obs <- obs_shape(fit) 
fit$observed_shape <- obs$fitShape 
fit$asymptote <- obs$asymp 
class(fit) <- 'sars' 
attr(fit, 'type') <- 'fit' 
return(fit) 
}#end of sar_p2
