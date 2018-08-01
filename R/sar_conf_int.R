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
  
  wei <- fit$details$weights_ics
  nams <- as.vector(names(wei))
  
  #pick model name based on model weight
  sm <- sample(nams, n, prob = wei)
  
  
  
  
  
  
}