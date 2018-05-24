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
    fit_df <- round(data.frame(Area = object$data$A, Fitted = object$calculated), 2) 
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
    #normality
    if (object$normaTest$test == "shapiro"|| object$normaTest$test == "lillie" ||
        object$normaTest$test == "kolmo"){
      normP <- object$normaTest[[2]]$p.value
    } else{
      normP <- "No normality test undertaken"
    }
    #homogeneity
    if (object$homoTest$test == "cor.area" || object$homoTest$test == "cor.fitted"){
      homoP <- object$homoTest[[2]]$p.value
    } else{
      homoP <- "No homogeneity test undertaken"
    }
    res <- list("Model" = name, "residuals" = round(resid, 1), "Parameters" = pars_tab, 
                "parNames" = parN, "formula" = formula, "AIC" = round(ic, 2), "AICc" = round(ic2, 2), "BIC" = round(bi, 2),
                "R2" = round(R2, 2), "R2a" = round(R2a, 2), "observed_shape" = shape, "asymptote" = asymp, 
                "convergence" = conv, "Normality_test_P" = normP, "Homogeneity_test_P" = homoP)
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


