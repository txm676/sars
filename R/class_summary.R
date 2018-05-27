#' Summarising the results of the model fitting functions
#'
#' @description S3 method for class 'sars'. summary.sars creates summary
#'   statistics for objects of class 'sars'. The exact summary statistics
#'   computed depends on the 'Type' attribute (e.g. 'lin_pow') of the 'sars'
#'   object. The summary method generates more useful information
#'   for the user than the standard model fitting functions. Another S3 method
#'   (\code{print.summary.sars}; not documented) is used to print the output.
#' @param object An object of class 'sars'.
#' @param ...	further arguments passed to or from other methods.
#' @return The summary function returns an object of class "summary.sars". A
#'   print function is used to obtain and print a summary of the model fit
#'   results.
#' @examples
#' data(galap)
#' fit <- lin_pow(galap, con = 1, compare = T)
#' summary(fit)
#' @export

summary.sars <- function(object){
  if (attributes(object)$type == "lin_pow"){
    rownames(object$Model$coefficients) <- c("LogC", "z")
    fit_df <- round(data.frame(Area = object$data$A, Fitted = object$calculated), 2) 
    res <- list("Model" = object$Model, df = fit_df,
                "Normality_test" = object$normaTest, "Homogeneity_test" = object$homoTest)
    if ("power" %in% names(object)){
      cp <- object$power$par[1]
      zp <- object$power$par[2]
      res$power <- round(c("logc" = cp, "z" = zp), 2)
    }
  }
  
  if (attributes(object)$type == "fit"){
    name <- object$model$name
    resid <- object$residuals
    pars_tab <- object$sigConf
    parN <- object$model$parNames
    formula <- object$model$formula
    ic <- object$AIC
    ic2 <- object$AICc
    bi <- object$BIC
    R2 <- object$R2
    R2a <- object$R2a
    shape <- object$observed_shape
    asymp <- object$asymptote
    conv <- object$verge
    neg <- ifelse(any(object$calculated < 0), 1, 0)
    res <- list("Model" = name, "residuals" = round(resid, 1), "Parameters" = pars_tab, 
                "parNames" = parN, "formula" = formula, "AIC" = round(ic, 2), "AICc" = round(ic2, 2), "BIC" = round(bi, 2),
                "R2" = round(R2, 2), "R2a" = round(R2a, 2), "observed_shape" = shape, "asymptote" = asymp, 
                "convergence" = conv, "Normality_test" = object$normaTest, "Homogeneity_test" = object$homoTest, 
                "Negative_values" = neg)
  }
  
  
  if (attributes(object)$type == "fit_collection"){ 
    return(cat("\n", "No summary method for fit_collection", "\n", sep = ""))
    
  }
  
  if (attributes(object)$type == "multi"){
    Mods <- as.vector(object$details$mod_names)
    nf <- as.vector(object$details$no_fit)
    cri <- object$details$ic
    ranks <- object$details$weights_ics
    df <- data.frame("Model" = names(ranks), "Weight" = as.vector(ranks))
    df$IC <- vapply(object$details$fits, function(x){x[[cri]]}, FUN.VALUE = numeric(1))
    colnames(df)[3] <- paste(cri)
    df$R2 <- vapply(object$details$fits, function(x){x$R2}, FUN.VALUE = numeric(1))
    df$R2a <- vapply(object$details$fits, function(x){x$R2a}, FUN.VALUE = numeric(1))
    df$Shape <- vapply(object$details$fits, function(x){x$observed_shape}, FUN.VALUE = character(1))
    df$Asymptote <- vapply(object$details$fits, function(x){x$asymptote}, FUN.VALUE = logical(1))
    df <- df[order(-df$Weight),]
    df[, 2:5] <- round(df[, 2:5], 3)
    rownames(df) <- NULL
    res <- list("Models" = Mods, "Criterion" = cri, "Model_table" = df, "no_fit" = nf)
  }
  
  class(res) <- "summary.sars"
  attr(res, "type") <- attr(object, "type")
  return(res)
}


