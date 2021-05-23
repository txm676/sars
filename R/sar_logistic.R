#' Fit the Logistic(Standard) model

#' @description Fit the Logistic(Standard) model to SAR data.
#' @usage sar_logistic(data, start = NULL, grid_start = 'partial',
#'   grid_n = NULL, normaTest = 'none',
#'   homoTest = 'none', homoCor = 'spearman', verb = TRUE)
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param start NULL or custom parameter start values for the optimisation algorithm.
#' @param grid_start Should a grid search procedure be implemented to test multiple starting parameter values. Can be one of 'none', 'partial' or 'exhaustive' The default is set to 'partial'.
#' @param grid_n If \code{grid_start = exhaustive}, the number of points sampled in the starting parameter space.
#' @param normaTest The test used to test the normality of the residuals of the
#'   model. Can be any of 'lillie' (Lilliefors test
#', 'shapiro' (Shapiro-Wilk test of normality), 'kolmo'
#'   (Kolmogorov-Smirnov test), or 'none' (no residuals normality test is undertaken; the default).
#' @param homoTest The test used to check for homogeneity of the residuals of
#'   the model. Can be any of 'cor.fitted' (a correlation of the residuals with
#'   the model fitted values), 'cor.area' (a correlation of the
#'   residuals with the area values), or 'none' (no residuals homogeneity test is undertaken; the default).
#' @param homoCor The correlation test to be used when \code{homoTest !='none'}. Can be any of 'spearman'
#'   (the default), 'pearson', or 'kendall'.
#' @param verb Whether or not to print certain warnings (default = TRUE)
#' @details The model is fitted using non-linear regression. The model parameters are estimated
#'   by minimizing the residual sum of squares with an unconstrained Nelder-Mead optimization algorithm
#'   and the \code{\link{optim}} function. To avoid numerical problems and speed up the convergence process,
#'   the starting values used to run the optimization algorithm are carefully chosen. However, if this does
#' not work, custom values can be provided (using the \code{start} argument), or a more comprehensive search
#'   can be undertaken using the \code{grid_start} argument. See the vignette for more information.
#'   The fitting process also determines the observed shape of the model fit,
#'   and whether or not the observed fit is asymptotic (see Triantis et al. 2012 for further details).

#'   Model validation can be undertaken by assessing the normality (\code{normaTest}) and homogeneity (\code{homoTest})
#'   of the residuals and a warning is provided in \code{\link{summary.sars}} if either test is chosen and fails.

#'   A selection of information criteria (e.g. AIC, BIC) are returned and can be used to compare models
#'   (see also \code{\link{sar_average}}).
#'   
#'   As grid_start has a random component, when \code{grid_start != 'none'} in your model fitting, you can
#'    get slightly different results each time you fit a model
#'   
#'    The parameter confidence intervals returned in sigConf are just simple confidence intervals, calculated as 2 * standard error.
#' @importFrom stats lm quantile
#' @return A list of class 'sars' with the following components: 
#'   \itemize{
#'     \item{par} { The model parameters}
#'     \item{value} { Residual sum of squares}
#'     \item{counts} {  The number of iterations for the convergence of the fitting algorithm}
#'     \item{convergence} { Numeric code returned from optim indicating model convergence (0 = converged)}
#'     \item{message} { Any message from the model fit algorithm}
#'     \item{hessian} { A symmetric matrix giving an estimate of the Hessian at the solution found}
#'     \item{verge} { Logical code indicating that optim model convergence value is zero}
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
#'     \item{asymptote} { A logical value indicating whether the observed fit is asymptotic}
#'     \item{neg_check} { A logical value indicating whether negative fitted values have been returned}}

#'   The \code{\link{summary.sars}} function returns a more useful summary of
#'   the model fit results, and the \code{\link{plot.sars}} plots the model fit.
#' @references Triantis, K.A., Guilhaumon, F. & Whittaker, R.J. (2012) The island species-area
#'   relationship: biology and statistics. Journal of Biogeography, 39, 215-231.
#' @examples
#' data(galap)
#' fit <- sar_logistic(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_logistic <- function(data, start = NULL, 
grid_start = "partial", grid_n = NULL, 
normaTest =  "none", homoTest = 
"none", homoCor = "spearman", verb = TRUE){
if (!(is.matrix(data) | is.data.frame(data)))  
stop('data must be a matrix or dataframe')
if (is.matrix(data)) data <- as.data.frame(data)
if (anyNA(data)) stop('NAs present in data')
normaTest <- match.arg(normaTest, c('none', 'shapiro', 'kolmo',
'lillie')) 
homoTest <- match.arg(homoTest, c('none', 'cor.area',
'cor.fitted')) 
if (homoTest != 'none'){
homoCor <- match.arg(homoCor, c('spearman', 'pearson',
'kendall')) 
}
if (!(grid_start %in% c('none', 'partial', 'exhaustive'))){
stop('grid_start should be one of none, partial or exhaustive')
}
if (grid_start == 'exhaustive'){
  if (!is.numeric(grid_n))
  stop('grid_n should be numeric if grid_start == exhaustive')
  }
if (!is.logical(verb)){
stop('verb should be logical')
}
data <- data[order(data[,1]),]
colnames(data) <- c('A','S')
#check for all equal richness values (particuarly zeros)
xr <- range(data$S)/mean(data$S)
if (isTRUE(all.equal(xr[1], xr[2]))) {
  if (data$S[1] == 0){
   warning('All richness values are zero: parameter estimates of',
           ' non-linear models should be interpreted with caution')
     } else{
       warning('All richness values identical')
     }}
#Standard Logistic function
#Formula given in Tjorve (2003)
#d confirmed as asymptote in Tjorve (2003)
#Starting parameter method taken from:
#https://bscheng.com/2014/05/07/modeling-logistic-growth-data-in-r/
#Shape given as sigmoid, but it can have observed shapes of linear, convex up
#and convex down. 
#Pars all set to Rplus, similar to heleg and mmf. d has to be as is the 
#asymptote, and have tested the other two - z should be Rplus but with c it is
#trickier. For some datasets, setting c to R results in a better fit, but then 
#this also means for other datasets it converges on a worse fit. Having it as
#Rplus seems to allow it to better fit the "typical" sigmoid shape so have left
#it as that for now.
model <- list(
  name=c("Logistic(Standard)"),
  formula=expression(S==d/(1 + exp(-z*A + c))),
  exp=expression(d/(1 + exp(-z*A + c))),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","Rplus","Rplus"),
  init=function(data){
    if (any(data$S == 0)){
      p <- (data$S + 0.1)/(max(data$S) + 1)
    } else{
      p <- data$S/(max(data$S) + 1)
    }
    logitR <- log(p / (1 - p)) #logit function (same as car::logit)
    cc <- coef(lm(logitR ~ data$A))
    d1 <- max(data$S) + 1
    #use abs, because in the formulation they use in that website, the sign
    #of the c parameter is reversed (seemingly negative most of the time)
    c(d1, cc[2], abs(cc[1]))
  }
)
model <- compmod(model)
fit <- get_fit(model = model, data = data, start = start,  
grid_start = grid_start, grid_n = grid_n, algo = 'Nelder-Mead', 
       normaTest =  normaTest, homoTest = homoTest, 
       homoCor = homoCor, verb = verb)
if(is.na(fit$value)){
  return(list(value = NA))
}else{ 
  obs <- obs_shape(fit, verb = verb)
  fit$observed_shape <- obs$fitShape
  fit$asymptote <- obs$asymp
  fit$neg_check <- any(fit$calculated < 0)
  class(fit) <- 'sars'
  attr(fit, 'type') <- 'fit'
  return(fit)
}
}#end of sar_logistic
