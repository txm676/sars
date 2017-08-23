#' Summarising the results of the model fitting functions
#'
#' @description S3 method for class 'mmSAR2'. summary.mmSAR2 creates summary
#'   statistics for objects of class mmSAR2. The exact summary statistics
#'   computed depends on the 'Type' attribute (e.g. 'lin_pow') of the mmSAR2
#'   object (see below). The summary method generates more useful information
#'   for the user than the standard model fitting functions. Another S3 method
#'   (print.summary.mmSAR2; not documented) is used to print the output.
#' @param object An object of class 'mmSAR2'.
#' @param ...	further arguments passed to or from other methods.
#' @return The summary function returns an object of class "summary.mmSAR2". A
#'   print function is used to obtain and print a summary of the model fit
#'   results.
#'
#'   An object of class "summary.mmSAR2" is a list containing the following
#'   components:
#'   \itemize{
#'     \item{"Summary"}{   A vector of model fit values (e.g. the c, z,
#'     z-significance, and R2 value for the linear power model fit).}
#'     \item{"df"}{   A dataframe containing the island areas and the models' fitted values.}
#'    }
#' @examples
#' data(galap)
#' fit <- lin_pow(galap, a = 1, s = 2, con = 1)
#' summary(fit)
#' @export


summary.sars <- function(object){

  if (attributes(object)$type == "lin_pow"){
    object2 <- object$Model
    logc <- object2$coefficients[1, 1]
    z <- object2$coefficients[2, 1]
    z.sig <- object2$coefficients[2, 4]
    r2 <- object2$r.squared
    md_res <- round(c(logc, z, z.sig, r2), 2) 
    names(md_res) <- c("logc", "z", "z.sig", "r2")
    fit_df <- round(data.frame(Area = object$Area, Fitted = object$Fitted), 2) 
    res <- list(Summary = md_res, df = fit_df)
    if(length(object) == 4){
      cp <- object[[4]]$par[1]
      zp <- object[[4]]$par[2]
      res$power <- round(c("logc" = cp, "z" = zp), 2)
    }
  }
  
  if (attributes(object)$type == "fit"){
    name <- object$model$name
    resid <- object$residuals
    pars <- object$par
    parN <- object$model$parNames
    ic <- object$AIC
    ic2 <- object$AICc
    bi <- object$BIC
    R2 <- object$R2
    R2a <- object$R2a
    shape <- object$observed_shape
    asymp <- object$asymptote
    res <- list("Model" = name, "residuals" = round(resid, 1), "Parameters" = round(pars, 2), 
                "parNames" = parN, "AIC" = round(ic, 2), "AICc" = round(ic2, 2), "BIC" = round(bi, 2),
                "R2" = round(R2, 2), "R2a" = round(R2a, 2), "observed_shape" = shape, "asymptote" = asymp)
  }
  
  class(res) <- "summary.sars"
  attr(res, "type") <- attr(object, "type")
  return(res)
}



