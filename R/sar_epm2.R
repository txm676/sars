#' Fit the Extended Power model 2 model

#' @description Fit the Extended Power model 2 model to SAR data
#' @usage sar_epm2(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_epm2(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_epm2 <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
# EXTENDED POWER MODEL 2 (TJORVE 2009)
model <- list(
  name=c("Extended Power model 2"),
  formula=expression(S==c*A^(z-(d/A))),
  exp=expression(c*A^(z-(d/A))),
  shape="sigmoid",
  asymp=function(pars)FALSE,
  custStart=function(data)c(5,.25,1),
  #limits for parameters
  parLim = c("Rplus","unif","R"),
  #initials values function
  init=function(data){return(c(0,0,0))}
)

model <- compmod(model) 
fit <- rssoptim(model = model, data = data, custstart = start, algo = 'Nelder-Mead') 
obs <- obs_shape(fit) 
fit$observed_shape <- obs$fitShape 
fit$asymptote <- obs$asymp 
class(fit) <- 'sars' 
attr(fit, 'type') <- 'fit' 
return(fit) 
}#end of sar_epm2
