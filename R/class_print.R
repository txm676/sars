
#' @export



print.summary.sars <- function(object){
  
  if (attributes(object)$type == "lin_pow"){
    cat("Model = ","Log-log power", "\n", "Logc =", object$Summary[1], "\n", 
        "c =", exp(object$Summary[1]), "\n",
        "z =", object$Summary[2], "\n", "z.sig =", object$Summary[3], 
    "\n", "R2 =", object$Summary[4], "\n")
    cat("\n")
    if (length(object) == 3){
      cat("Power (non-linear) parameters:", "\n",
          "c =", object$power[1], "\n",
          "z =", object$power[2])
    }
  }
  
  if (attributes(object)$type == "fit"){
    cat("\n", "Model: ","\n", object$Model, "\n", sep = "")
    cat("\n", "Call: ","\n", as.character(object$formula), "\n", sep = "")
    cat("\n", "Did the model converge: ", object$convergence , "\n", sep = "")
    cat("\n", "Residuals: ", "\n", sep = "")
    print(stats::quantile(object$residuals))
    cat("\n", "Parameters: ", "\n", sep = "")
    mm <- object$Parameters
    rownames(mm) <- object$parNames
    stats::printCoefmat(mm)
    cat("\n", "R-squared: ", object$R2 , ", Adjusted R-squared: ", object$R2a, "\n", sep = "")
    cat("AIC: ", object$AIC , ", AICc: ", object$AICc, ", BIC: ", object$BIC, "\n", sep = "")
    cat("Observed shape: ", object$observed_shape, ", Asymptote: ", object$asymptote, "\n", "\n", sep = "")
  }

  
}
  

#' @export
#' 

print.sars <- function(object){
  
  if (attributes(object)$type == "fit"){ 
    cat("\n", "Model: ","\n", object$model$name, "\n", sep = "")
    cat("\n", "Call: ","\n", as.character(object$model$formula), "\n", sep = "")
    cat("\n", "Coefficients: ", "\n", sep = "")
    print(object$par)
    cat("\n")
  }
  
  if (attributes(object)$type == "fit_collection"){ 
    cat("\n", "This is a fit collection", "\n", sep = "")
    cat("\n", length(object), " models contained in the fit collection: ","\n", sep = "")
    cat(unlist(lapply(ff, function(x) x$model$name)), "\n", "\n")

  }
  
  
  
  
  
  
  
  
  
  
  
}