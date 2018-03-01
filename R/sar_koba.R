#' Fit the Kobayashi model

#' @description Fit the Kobayashi model to SAR data
#' @usage sar_koba(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in optim
#'   to run a grid search.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_koba(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_koba <- function(data = galap, start = NULL, grid_start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
# Kobayashi logarithmic (KOBAYASHI 1975), convex upward, no asymptote
model <- list(
  name = c("Kobayashi"),
  formula = expression(S == c*log(1+ A/z)),
  exp = expression(c*log(1 + A/z)),
  shape = "convex",
  asymp = function(pars)FALSE,
  #limits for parameters
  parLim = c("R","Rplus"),
  custStart = function(data) c(5,0.1),
  #initials values function
  init = function(data){
    c(max(data$S), max(data$S))
  }
)



################## test init

#from mKobayashi 1975

# S = c if A = (e - 1) * z 
# with e = natural base log ~ 2.718

# p269
# c is the number of species occuring in it's caracteristic area (e - 1) * z





model <- compmod(model) 
fit <- get_fit(model = model, data = data, start = start, grid_start = grid_start, algo = 'Nelder-Mead', verb = TRUE) 
if(is.na(fit$value)){ 
  return(NA) 
}else{ 
  obs <- obs_shape(fit) 
  fit$observed_shape <- obs$fitShape 
  fit$asymptote <- obs$asymp 
  class(fit) <- 'sars' 
  attr(fit, 'type') <- 'fit' 
  return(fit) 
} 
}#end of sar_koba
