#' Fit the linear model

#' @description Fit the linear model to SAR data
#' @usage sar_linear(data,)
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_linear(galap)
#' summary(fit)
#' plot(fit)
#' @export 

sar_linear <- function(data = galap){
  if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
  if (is.matrix(data)) data <- as.data.frame(data) 
  if (base::anyNA(data)) stop('NAs present in data') 
  data <- data[order(data[,1]),] 
  colnames(data) <- c('A','S') 
  #standard linear regression
  mod <- lm(S ~ A, data = data)
  fit <- vector("list", length = 13)
  fit$par <- mod$coefficients
  names(fit$par) <- c("c", "m")
  res <- as.vector(mod$residuals)
  fit$value <- sum(res^2)#residual sum of squares
  fit$verge <- TRUE#should always converge
  fit$data <- data
  fit$model$name <- "Linear model"
  fit$model$formula <- "S == c + m*A"
  fit$model$parNames <- c("c", "m")
  fit$calculated <- as.vector(mod$fitted.values)
  fit$residuals <- res
  #AIC, R2 etc
  n <- dim(data)[1]
  value <- sum(mod$residuals^2)
  #A common mistake when calculating the number of parameters is failure to include error.
  #This is because it is not normally thought of as a parameter as, strictly speaking, 
  #you're not really predicting it. As a results the number of parameters in a standard linear 
  #equation (y=mx+c) is 3 (mx, c, and error) rather than 2.
  P <- 3
  fit$AIC <- n * log(value / n) + 2 * P
  fit$AICc <- (n * log(value / n)) + (2*P*(n / (n - P - 1)))
  fit$BIC <- (n *log(value / n)) + (log(n) * P)
  #R2 (Kvaleth, 1985, Am. Statistician)
  fit$R2 <-  1 - ( (value) /  sum((data$S - mean(data$S))^2) )
  #R2a (He & Legendre 1996, p724)
  fit$R2a <-  1 - ( ((n-1)*(value)) /  ((n-P)*sum((data$S - mean(data$S))^2)) )
  fit$sigConf <- cbind(summary(mod)$coefficients, confint(mod))#confidence intervals calculated using in built function
  fit$observed_shape <- "linear" 
  fit$asymptote <- FALSE
  class(fit) <- 'sars' 
  attr(fit, 'type') <- 'fit' 
  return(fit) 
}#end of sar_linear

#nor = shapiro.test(resi)
#pear = cor.test(resi,dat$A)

#normaTest <- list(shapiro=nor)
#homoTest <- list(cor.area=pear)
#fit <- list(value=value,verge=69,normaTest=normaTest,homoTest=homoTest,calculated=mod$fitted.values,residuals=resi,AICc=AICc,AIC=AIC,BIC=BIC,R2=R2,R2a=R2a)

