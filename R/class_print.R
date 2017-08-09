

print.summary.sars <- function(object){
  
  if (attributes(object)$type == "linpow"){
    cat("Model = ", "Log-log power", "\n", "c =", x$Summary[1], "\n", 
        "z =", x$Summary[2], "\n")
    cat("z.sig =", x$Summary[3], "\n", "R2 =", x$Summary[4], "\n")
    cat("\n")
    if (length(object == 3)){
      cat("Power (non-linear) parameters:", "\n",
          "c =", object$power[1], "z =", object$power[1])
    }
  }
}
  