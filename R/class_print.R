#' @inherit  summary.sars
#' @export

print.summary.sars <- function(x, ...){
  if (attributes(x)$Type == "lin_pow") cat(" Data = ",attributes(x)$Dataset, "\n",
                                           "Model = ", "log-log power", "\n")
  cat("c =", x$Summary[1], "\n")
  cat("z =", x$Summary[2], "\n")
  cat("z.sig =", x$Summary[3], "\n")
  cat("R2 =", x$Summary[4], "\n")
  cat("\n")
  #cat("Fitted values:", "\n", x$df)
}


