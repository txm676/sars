#' Fit the Monod model

#' @description Fit the Monod model to SAR data.
#' @usage sar_monod(data, start = NULL, grid_start = FALSE, grid_n = NULL, normaTest = 'none',
#'   homoTest = 'none', homoCor = 'spearman')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param start NULL or custom parameter start values for the optimisation algorithm.
#' @param grid_start Logical argument specifying whether a grid search procedure should be implemented to test multiple starting parameter values. The default is set to FALSE, but for certain models (e.g. Gompertz, Chapman Richards), we advice using it to ensure an optimal fit.
#' @param grid_n If \code{grid_start = TRUE}, the number of points sampled in the model parameter space.
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
#'   (see also \code{\link{sar_average}})
#' @importFrom stats lm quantile
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
#'     \item{asymptote} { A logical value indicating whether the observed fit is asymptotic}
#'     \item{neg_check} { A logical value indicating whether negative fitted values have been returned}}

#'   The \code{\link{summary.sars}} function returns a more useful summary of
#'   the model fit results, and the \code{\link{plot.sars}} plots the model fit.
#' @references Triantis, K.A., Guilhaumon, F. & Whittaker, R.J. (2012) The island species-area
#'   relationship: biology and statistics. Journal of Biogeography, 39, 215-231.
#' @examples
#' data(galap)
#' fit <- sar_monod(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_monod <- function(data, start = NULL, grid_start = FALSE, 
grid_n = NULL, normaTest =  "none", homoTest = "none", homoCor = "spearman"){
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
if (!is.logical(grid_start)) stop('grid_start should be logical')
if (grid_start){
  if (!is.numeric(grid_n))
  stop('grid_n should be numeric if grid_start == TRUE')
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
# MONOD CURVE (MONOD 1950, Willimas et al. 2009 formula)
#Have checked and this formula is equivalent to that in 
#Tjorve and generates the same output etc.
model <- list(
  name=c("Monod"),
  formula=expression(S==d/(1+c*A^(-1))),
  exp=expression(d/(1+c*A^(-1))),
  shape="convex",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","Rplus"),
  custStart=function(data)c(stats::quantile(data$A,c(0.25)),max(data$S)),
  #initials values function
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    d=as.double(max(data$A)+max(data$S)/4)
    c=data[[1]]*(d/data$S - 1)
    c(d,stats::quantile(c,c(0.25)))
  }
)

model <- compmod(model)
fit <- get_fit(model = model, data = data, start = start,  
grid_start = grid_start, grid_n = grid_n, algo = 'Nelder-Mead', 
       normaTest =  normaTest, homoTest = homoTest, 
       homoCor = homoCor)
if(is.na(fit$value)){
  return(list(value = NA))
}else{ 
  obs <- obs_shape(fit)
  fit$observed_shape <- obs$fitShape
  fit$asymptote <- obs$asymp
  fit$neg_check <- any(fit$calculated < 0)
  class(fit) <- 'sars'
  attr(fit, 'type') <- 'fit'
  return(fit)
}
}#end of sar_monod
