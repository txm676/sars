#' @importFrom stats printCoefmat
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
    cat("\n", "Model:","\n", object$Model, "\n", sep = "")
    mm <- matrix(object$Parameters, nrow = 1, ncol = length(object$par))
    colnames(mm) <- object$parNames
    rownames(mm) <- "Estimates"
    cat("\n", "Parameters: ", "\n", sep = "")
    printCoefmat(mm)
    cat("\n", "R-squared: ", object$R2 , ", Adjusted R-squared: ", object$R2a, "\n", sep = "")
    cat("AIC: ", object$AIC , ", AICc: ", object$AICc, ", BIC: ", object$BIC, "\n", sep = "")
    cat("Observed shape: ", object$observed_shape, ", Asymptote: ", object$asymptote, "\n", "\n", sep = "")
  }
  
  
}
  