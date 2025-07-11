#' Fit the Chapman Richards model

#' @description Fit the Chapman Richards model to SAR data.
#' @usage sar_chapman(data, start = NULL, grid_start = 'partial',
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
#'     \item \strong{par}  The model parameters
#'     \item \strong{value}  Residual sum of squares
#'     \item \strong{counts}   The number of iterations for the convergence of the fitting algorithm
#'     \item \strong{convergence}  Numeric code returned from optim indicating model convergence (0 = converged)
#'     \item \strong{message}  Any message from the model fit algorithm
#'     \item \strong{hessian}  A symmetric matrix giving an estimate of the Hessian at the solution found
#'     \item \strong{verge}  Logical code indicating that optim model convergence value is zero
#'     \item \strong{startValues}  The start values for the model parameters used in the optimisation
#'     \item \strong{data}  Observed data
#'     \item \strong{model}  A list of model information (e.g. the model name and formula)
#'     \item \strong{calculated}   The fitted values of the model
#'     \item \strong{residuals}  The model residuals
#'     \item \strong{AIC}  The AIC value of the model
#'     \item \strong{AICc}  The AICc value of the model
#'     \item \strong{BIC}  The BIC value of the model
#'     \item \strong{R2}  The R2 value of the model
#'     \item \strong{R2a}  The adjusted R2 value of the model
#'     \item \strong{sigConf}  The model coefficients table
#'     \item \strong{normaTest}  The results of the residuals normality test
#'     \item \strong{homoTest}  The results of the residuals homogeneity test
#'     \item \strong{observed_shape}  The observed shape of the model fit
#'     \item \strong{asymptote}  A logical value indicating whether the observed fit is asymptotic
#'     \item \strong{neg_check}  A logical value indicating whether negative fitted values have been returned}

#'   The \code{\link{summary.sars}} function returns a more useful summary of
#'   the model fit results, and the \code{\link{plot.sars}} plots the model fit.
#' @references Triantis, K.A., Guilhaumon, F. & Whittaker, R.J. (2012) The island species-area
#'   relationship: biology and statistics. Journal of Biogeography, 39, 215-231.
#' @examples
#' data(galap)
#' fit <- sar_chapman(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_chapman <- function(data, start = NULL, 
grid_start = "partial", grid_n = NULL, 
normaTest =  "none", homoTest = 
"none", homoCor = "spearman", verb = TRUE){
if (!(is.matrix(data) | is.data.frame(data)))  
stop('data must be a matrix or dataframe')
data <- as.data.frame(data)
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
#Chapman–Richards 3 S = a [1 − exp(−bA)]c Flather (1996)
#note there was an error in our original formula (corrected
#Nov 2020)
model <- list(
  name = c("Chapman Richards"),
  formula = expression(S == d * (1 - exp(-z*A))^c),
  exp = expression(d * (1 - exp(-z*A))^c),
  shape = "sigmoid",
  asymp = function(pars)pars["d"],
  #limits for parameters
  parLim  =  c("Rplus","R","R"),
  #initials values function
  init = function(data){
    d=max(data$S) 
    Z=(-log((-data$S/(max(data$S)+1))+1))/data$A
    z = mean(Z)
    c(d,z,0.5)}
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
}#end of sar_chapman
