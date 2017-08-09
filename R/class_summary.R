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
#' @importFrom dplyr %>%
#' @export


summary.sars <- function(object){

  if (attributes(object)$type == "linpow"){
    object2 <- object$Model
    logc <- object2$coefficients[1, 1]
    z <- object2$coefficients[2, 1]
    z.sig = object2$coefficients[2, 4]
    r2 <- object2$r.squared
    md_res <- c(logc, z, z.sig, r2) %>% round(2)
    names(md_res) <- c("logc", "z", "z.sig", "r2")
    fit_df <- round(data.frame(Area = object$Area, Fitted = object$Fitted), 2) 
    res <- list(Summary = md_res, df = fit_df)
    if(length(object) == 4){
      cp <- object[[4]]$par[1]
      zp <- object[[4]]$par[2]
      res$power <- round(c("logc" = cp, "z" = zp), 2)
    }
  }
  
  
  
  
  
  
  
  
  class(res) <- "summary.sars"
  attr(res, "type") <- attributes(object)$type
  return(res)
}



