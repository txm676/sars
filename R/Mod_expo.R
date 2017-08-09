# EXPONENTIAL MODEL (GLEASON 1922)

#' @export

sar_expo <- function(data, custstart = NULL, normtest = "lillie"){
  
  if (!(is.matrix(data) || is.data.frame(data))) stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop("NAs present in data")
  normtest <- match.arg(normtest, c("none", "shapiro", "kolmo", "lillie"))
  
  data <- data[order(data[,1]),]
  colnames(data) <- c("A", "S")
  

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
  
  fit <- rssoptim(model, data, custstart, normtest, algo = "Nelder-Mead")
  
  fit$model <- model
  
  class(fit) <- "sars"
  attr(fit, "type") <- "fit"
  return(fit)
}

