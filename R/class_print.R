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
}
  