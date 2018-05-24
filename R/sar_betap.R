#' Fit the Beta-P cumulative model

#' @description Fit the Beta-P cumulative model to SAR data
#' @usage sar_betap(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param grid_start NULL or the number of points sampled in the model parameter space
#'   or FALSE to prevent any grid start after a fail in inital optimization
#'   to run a grid search.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_betap(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_betap <- function(data = galap, start = NULL, grid_start = NULL, normaTest =  "lillie",
              homoTest = "cor.fitted"){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (base::anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#Beta-P function (cumulative)
model <- list(
  name=c("Beta-P cumulative"),
  exp=expression(d*(1-(1+(A/c)^z)^-f)),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R","R"),
  #initial values function
  init=function(data){c(max(data$S)+1,.5,.5,.5)}
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
}#end of sar_betap
