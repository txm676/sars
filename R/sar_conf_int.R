#' Generate a non parametric confidence interval for a SAR

#' @description Generate a non parametric confidence interval for a SAR
#' @usage sar_conf_int(obj, n_boot = 99)
#' @param obj An object of class 'sar' or 'multi.sar'.
#' @param n_boot the number of bootstrap sample used to construct the confidence interval.
#' @return 
#' @examples
#' #multi_sar fit
#' data(galap)
#' fit <- sar_multi(data = galap)
#' conf_fit <- sar_conf_int(fit)
#' #simple sar fit
#' data(galap)
#' fit <- sar_power(data = galap)
#' conf_fit <- sar_conf_int(fit)

#' @export


#n = number of iterations

sar_conf_int <- function(fit, n, normaTest = "lillie",
                         homoTest = "cor.fitted",
                         neg_check = TRUE,
                         alpha_normtest = 0.05,
                         alpha_homotest = 0.05){
  
  if (!"multi" %in% class(fit)) stop ("class of 'fit' should be 'multi'")
  if (length(fit$details$mod_names) < 2) stop ("less than two models in the sar multi object")
  
  #model names for matching
  
  x1 <- c("Power", "PowerR", "Extended_Power_model_1", "Extended_Power_model_2", "Persistence_function_1", 
          "Persistence_function_2", "Exponential", "Kobayashi", "MMF", "Monod", "Negative_exponential", 
          "Chapman_Richards", "Cumulative_Weibull_3_par.", "Asymptotic_regression", "Rational_function", 
          "Gompertz", "Cumulative_Weibull_4_par.", "Beta-P_cumulative", "Heleg(Logistic)", "Linear_model")
  
  x2 <- c("sar_power(", "sar_powerR(", "sar_epm1(", "sar_epm2(", "sar_p1(", "sar_p2(", 
          "sar_expo(", "sar_koba(", "sar_mmf(", "sar_monod(", "sar_negexpo(", 
          "sar_chapman(", "sar_weibull3(", "sar_asymp(", "sar_ratio(", "sar_gompertz(", 
          "sar_weibull4(", "sar_betap(", "sar_heleg(", "sar_linear(")
  
  #observed data
  dat <- fit$details$fits[[1]]$data
  
  #weights and model names
  wei <- fit$details$weights_ics
  nams <- as.vector(names(wei))
  
  #pick model name based on model weight
  sm <- sample(nams, n, prob = wei)
  
  #select the expression of the selected model
  wn <- which(x1 %in% sm)
  w2 <- x2[wn]
  
  #fit the best model to observed data; extract fitted values and residuals
  me <- eval(parse(text = paste(w2, "dat)", sep = "")))
  meF <- me$calculated
  meR <- me$residuals
  
  
  #based on code from mmSAR in Rforge
  jacob <- matrix(nrow = nrow(me$data), ncol = 2)
  
  for (k in 1:nrow(jacob)) {
    jacob[k, ] <- numDeriv::jacobian(me$model$rss.fun, me$par, data = me$data[k, ], opt = FALSE)
  }
  
  jacobbis <- t(jacob) %*% jacob
  s <- svd(jacobbis)
  jacobbismun <- s$v %*% (diag(1 / s$d)) %*% (t(s$u))
  hatMat <- jacob %*% jacobbismun %*% t(jacob)
  matList <- list(jacob = jacob, hatMat = hatMat)
  
  #Residuals transformation from Davidson and Hinkley, 1997 "Bootstrap methods and their applications" p 259 eq (6.9)
  diagHatMat <- diag(hatMat)
  transResiduals <- meR - mean(meR)
  transResiduals <- transResiduals / sqrt(1 - diagHatMat)#checked and this gives same values as mmSAR
  
  #sample residuals with replacement until n matches data
  sr <- sample(transResiduals, nrow(me$data), replace = T)
  
  #add modified residuals to fitted values
  mf <- meF + sr

  
}