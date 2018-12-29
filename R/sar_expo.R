#' Fit the Exponential model

#' @description Fit the Exponential model to SAR data.
#' @usage sar_expo(data, start = NULL, grid_start = NULL, normaTest =  'lillie',
              
#'   homoTest = 'cor.fitted')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param start NULL or custom parameter start values for the optimisation algorithm.
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in initial optimization
#'   to run a grid search.
#' @param normaTest The test used to test the normality of the residuals of the
#'   model. Can be any of 'lillie' (Lilliefors Kolmogorov-Smirnov test; the
#'   default), 'shapiro' (Shapiro-Wilk test of normality), 'kolmo'
#'   (Kolmogorov-Smirnov test), or 'none' (no residuals normality test is undertaken).
#' @param homoTest The test used to check for homogeneity of the residuals of
#'   the model. Can be any of 'cor.fitted' (a correlation of the residuals with
#'   the model fitted values; the default), 'cor.area' (a correlation of the
#'   residuals with the area values), or 'none' (no residuals homogeneity test is undertaken).
#' @details The model is fitted using non-linear regression. The model parameters are estimated
#'   by minimizing the residual sum of squares with an unconstrained Nelder-Mead optimization algorithm
#'   and the \code{\link{optim}} function. To avoid numerical problems and speed up the convergence process,
#'   the starting values used to run the optimization algorithm are carefully chosen, or custom values can be provided
#'   using the argument \code{start}. The fitting process also determines the observed shape of the model fit,
#'   and whether or not the observed fit is asymptotic (see Triantis et al. 2012 for further details).

#'   Model validation is undertaken by assessing the normality (\code{normaTest}) and homogeneity (\code{homoTest})
#'   of the residuals and a warning is provided in \code{\link{summary.sars}} if either test is failed.

#'   A selection of information criteria (e.g. AIC, BIC) are returned and can be used to compare models
#'   (see also \code{\link{fit_collection}} and \code{\link{sar_multi}})
#' @return A list of class 'sars' with the following components: 
#'   \itemize{
#'     \item{par} { The model parameters}
#'     \item{value} { Residual sum of squares}
#'     \item{counts} {  The number of iterations for the convergence of the fitting algorithm}
#'     \item{convergence} { Numeric code indicating model convergence (0 = converged)}
#'     \item{message} { Any message from the model fit algorithm}
#'     \item{hessian} { A symmetric matrix giving an estimate of the Hessian at the solution found}
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
#' @references Triantis, K.A., Guilhaumon, F. & Whittaker, R.J. (2012) The island species-area
#'   relationship: biology and statistics. Journal of Biogeography, 39, 215-231.
#' @examples
#' data(galap)
#' fit <- sar_expo(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_expo <- function(data, start = NULL, grid_start = NULL, 
normaTest =  "lillie", homoTest = "cor.fitted"){
if (!(is.matrix(data) | is.data.frame(data)))  
stop('data must be a matrix or dataframe')
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
      par=stats::lm(S~A,semilog.data)$coefficients
      names(par)=c("c","z")
      par
    }
)

model <- compmod(model)
fit <- get_fit(model = model, data = data, start = start,  
grid_start = grid_start, algo = 'Nelder-Mead', 
       
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
}#end of sar_expo
