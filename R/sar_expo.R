#' Fit the Asymptotic regression model

#' @description Fit the Asymptotic regression model to SAR data
#' @usage sar_expo(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_expo(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_expo <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
# EXPONENTIAL MODEL (GLEASON 1922)
model <- list(
    name=c("Exponential"),
    formula=expression(S==c+z*log(A)),
    exp=expression(c+z*log(A)),
    shape="convex",
    asymp=function(pars)FALSE,
    #limits for parameters
    parLim = c("R","R"),
    custStart=function(data)c(5,.25),
    #initials values function
    init=function(data){
      semilog.data = data.frame(log(data$A),data$S)
      names(semilog.data)=c("A","S")
      par=lm(S~A,semilog.data)$coefficients
      names(par)=c("c","z")
      par
    }
)

model <- compmod(model) 
fit <- rssoptim(model = model, data = data, custstart = start, normtest = 'lillie', algo = 'Nelder-Mead') 
obs <- obs_shape(fit) 
fit$observed_shape <- obs$fitShape 
fit$asymptote <- obs$asymp 
class(fit) <- 'sars' 
attr(fit, 'type') <- 'fit' 
return(fit) 
}#end of sar_expo
