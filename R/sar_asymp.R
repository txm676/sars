#' Fit the Asymptotic regression model

#' @description Fit the Asymptotic regression model to SAR data
#' @usage sar_asymp(data, start = NULL, grid_start = NULL, normaTest =  'lillie',
#'   homoTest = 'cor.fitted')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param start 
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in inital optimization
#'   to run a grid search.
#' @param normaTest The test used to test the normality of the residuals of the
#'   model. Can be any of 'lillie' (Lilliefors Kolmogorov-Smirnov test; the
#'   default), 'shapiro' (Shapiro-Wilk test of normality), 'kolmo'
#'   (Kolmogorov-Smirnov test), or 'none' (no residuals normality test is undertaken).
#' @param homoTest The test used to check for homogeneity of the residuals of
#'   the model. Can be any of 'cor.fitted' (a correlation of the residuals with
#'   the model fitted values; the default), 'cor.area' (a correlation of the
#'   residuals with the area values), or 'none' (no residuals homogeneity test is undertaken).
#' @return A list of class 'sars' with the following components: 
#'   \itemize{
#'     \item{par} { The model parameters}
#'     \item{value} { Residual sum of squares}
#'     \item{counts} { **}
#'     \item{convergence} { Numeric code indicating model convergence (0 = converged)}
#'     \item{message} { Any message from the model fit algorithm}
#'     \item{hessian} { ***}
#'     \item{verge} { Logical code indicating model convergence}
#'     \item{startValues} { The start values for the model parameters used in the optimisation}
#'     \item{data} { Observed data}
#'     \item{model} { A list of model information (e.g. the model name and formula)}
#'     \item{calculated} {  The fitted values of the model}
#'     \item{residuals} { The model residuals}
#'     \item{AIC} { The AIC value of the model}
#'     \item{AICc} { The AICc value of the model}
#'     \item{BIC} { The BIC value of the model}
#'     \item{R2} { The R2 value of the model}
#'     \item{R2a} { The adjusted R2 value of the model}
#'     \item{sigConf} { The model coefficients table}
#'     \item{normaTest} { The results of the residuals normality test}
#'     \item{homoTest} { The results of the residuals homogeneity test}
#'     \item{observed_shape} { The observed shape of the model fit}
#'     \item{asymptote} { A logical value indicating whether the observed fit is asymptotic}}

#'   The \code{\link{summary.sars}} function returns a more useful summary of
#'   the model fit results, and the \code{\link{plot.sars}} plots the model fit.
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
