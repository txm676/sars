#' Summarising the results of the model fitting functions
#'
#' @description S3 method for class 'sars'. \code{summary.sars} creates
#'   summary statistics for objects of class 'sars'. The exact summary
#'   statistics computed depends on the 'Type' attribute (e.g. 'multi') of
#'   the 'sars' object. The summary method generates more useful information
#'   for the user than the standard model fitting functions. Another S3
#'   method (\code{print.summary.sars}; not documented) is used to print the
#'   output.
#' @param object An object of class 'sars'.
#' @param \dots Further arguments.
#' @return The \code{summary.sars} function returns an object of class
#'   "summary.sars". A print function is used to obtain and print a summary
#'   of the model fit results.
#'
#'   For a 'sars' object of Type 'fit', a list with 16 elements is returned
#'   that contains useful information from the model fit, including the model
#'   parameter table (with t-values, p-values and confidence intervals),
#'   model fit statistics (e.g. R2, AIC), the observed shape of the model and
#'   whether or not the fit is asymptotic, and the results of any additional
#'   model checks undertaken (e.g. normality of the residuals).
#'
#'   For a 'sars' object of Type 'multi', a list with 4 elements is returned:
#'   (i) a vector of the names of the models that were successfully fitted
#'   and passed any additional checks, (ii) a character string containing the
#'   name of the criterion used to rank models, (iii) a data frame of the
#'   ranked models, and (iv) a vector of the names of any models that were
#'   not fitted or did not pass any additional checks. In regards to (iii;
#'   \code{Model_table}), the dataframe contains the fit summaries for each
#'   successfully fitted model (including the value of the model criterion
#'   used to compare models, the R2 and adjusted R2, and the observed shape
#'   of the fit); the models are ranked in decreasing order of information
#'   criterion weight.
#'
#'   For a 'sars' object of Type 'lin_pow', a list with 5 elements is
#'   returned: (i) the model fit output from the \code{\link{lm}} function,
#'   (ii) the fitted values of the model, (iii) the observed data, and (iv
#'   and v) the results of the residuals normality and heterogeneity tests.
#'   If the argument \code{compare = TRUE} is used in \code{\link{lin_pow}},
#'   a sixth element is returned that contains the parameter values from the
#'   non-linear power model.
#' @examples
#' data(galap)
#' #fit a multimodel SAR and get the model table
#' mf <- sar_average(data = galap)
#' summary(mf)
#' summary(mf)$Model_table
#' #Get a summary of the fit of the linear power model
#' fit <- lin_pow(galap, con = 1, compare = TRUE)
#' summary(fit)
#' @export

summary.sars <- function(object, ...){
  if (attributes(object)$type == "lin_pow"){
    rownames(object$Model$coefficients) <- c("LogC", "z")
    fit_df <- round(data.frame(Area = object$data$A,
                               Fitted = object$calculated), 2)
    res <- list("Model" = object$Model, df = fit_df,
                "normaTest" = object$normaTest, "homoTest" = object$homoTest)
    if ("power" %in% names(object)){
      cp <- object$power$par[1]
      zp <- object$power$par[2]
      res$power <- round(c("logc" = cp, "z" = zp), 2)
    }
  }
  
  if (attributes(object)$type == "pred"){
    stop("No summary method for a 'sars' object of type 'pred'\n", sep = "")
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
    negCheck <- ifelse(any(object$calculated < 0), 1, 0)
    res <- list("Model" = name, "residuals" = round(resid, 1),
                "Parameters" = pars_tab,
                "parNames" = parN, "formula" = formula, "AIC" = round(ic, 2),
                "AICc" = round(ic2, 2), "BIC" = round(bi, 2),
                "R2" = round(R2, 2), "R2a" = round(R2a, 2),
                "observed_shape" = shape, "asymptote" = asymp,
                "convergence" = conv, "normaTest" = object$normaTest,
                "homoTest" = object$homoTest,
                "Negative_values" = negCheck)
  }


  if (attributes(object)$type == "fit_collection"){
    return(cat("\nNo summary method for fit_collection\n", sep = ""))
  }

  if (attributes(object)$type == "multi"){
    Mods <- as.vector(object$details$mod_names)
    nf <- as.vector(object$details$no_fit)
    cri <- object$details$ic
    ranks <- object$details$weights_ics
    df <- data.frame("Model" = names(ranks),
                     "Weight" = as.vector(ranks))
    df$IC <- vapply(object$details$fits,
                    function(x){x[[cri]]}, FUN.VALUE = numeric(1))
    colnames(df)[3] <- paste(cri)
    df$R2 <- vapply(object$details$fits,
                    function(x){x$R2}, FUN.VALUE = numeric(1))
    df$R2a <- vapply(object$details$fits,
                     function(x){x$R2a}, FUN.VALUE = numeric(1))
    df$Shape <- vapply(object$details$fits,
                    function(x){x$observed_shape}, FUN.VALUE = character(1))
    df$Asymptote <- vapply(object$details$fits,
                           function(x){x$asymptote}, FUN.VALUE = logical(1))
    df <- df[order(-df$Weight),]
    df[, 2:5] <- round(df[, 2:5], 3)
    rownames(df) <- NULL
    res <- list("Models" = Mods, "Criterion" = cri, "Model_table" = df,
                "no_fit" = nf)
  }

  class(res) <- "summary.sars"
  attr(res, "type") <- attr(object, "type")
  res
}
