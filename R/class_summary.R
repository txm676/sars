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
#' @param order The information criterion used to order the model summary table
#'   for objects of Type 'threshold' (Default = "BIC").
#' @param \dots Further arguments.
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
#'   For a 'sars' object of Type 'multi', a list with 5 elements is returned:
#'   (i) a vector of the names of the models that were successfully fitted and
#'   passed any additional checks, (ii) a character string containing the name
#'   of the criterion used to rank models, (iii) a data frame of the ranked
#'   models, (iv) a vector of the names of any models that were not fitted or
#'   did not pass any additional checks, and (v) a logical vector specifying
#'   whether the \code{\link{optim}} convergence code for each model that passed
#'   all the checks is zero. In regards to (iii; \code{Model_table}), the
#'   dataframe contains the fit summaries for each successfully fitted model
#'   (including the value of the model criterion used to compare models, the R2
#'   and adjusted R2, and the observed shape of the fit); the models are ranked
#'   in decreasing order of information criterion weight.
#'
#'   For a 'sars' object of Type 'lin_pow', a list with up to 7 elements is
#'   returned: (i) the model fit output from the \code{\link{lm}} function, (ii)
#'   the fitted values of the model, (iii) the observed data, (iv and v) the
#'   results of the residuals normality and heterogeneity tests, and (vi) the
#'   log-transformation function used. If the argument \code{compare = TRUE} is
#'   used in \code{\link{lin_pow}}, a 7th element is returned that contains the
#'   parameter values from the non-linear power model.
#'   
#'   For a 'sars' object of Type 'threshold', a list with three elements is
#'   returned: (i) the information criterion used to order the ranked model
#'   summary table ('order'), (ii) a model summary table (models are ranked
#'   using the 'order' argument), and (iii) details of any axes
#'   log-transformations undertaken. Note that in the model summary table, if
#'   log-area is used as the predictor, the threshold values will be on the log
#'   scale used. Thus it may be preferable to back-transform them (e.g. using
#'   \code{exp(th)} if natural logarithms are used) so that they are on the
#'   scale of untransformed area. Th1 and Th2 in the table are the threshold
#'   value(s), and seg1, seg2, seg3 provide the number of datapoints within each
#'   segment (for the threshold models); one-threshold models have two
#'   segements, and two-threshold models have three segments.
#' @examples
#' data(galap)
#' #fit a multimodel SAR and get the model table
#' mf <- sar_average(data = galap)
#' summary(mf)
#' summary(mf)$Model_table
#' #Get a summary of the fit of the linear power model
#' fit <- lin_pow(galap, con = 1, compare = TRUE)
#' summary(fit)
#' @importFrom stats logLik
#' @export

summary.sars <- function(object, order = "BIC", ...){
  if (attributes(object)$type == "lin_pow"){
    rownames(object$Model$coefficients) <- c("LogC", "z")
    fit_df <- round(data.frame(Area = object$data$A,
                               Fitted = object$calculated), 2)
    res <- list("Model" = object$Model, df = fit_df,
                "normaTest" = object$normaTest, "homoTest" = object$homoTest, 
                logT = object$logT)
    if ("power" %in% names(object)){
      cp <- object$power$par[1]
      zp <- object$power$par[2]
      res$power <- round(c("logc" = cp, "z" = zp), 2)
    }
  }
  
  if (attributes(object)$type == "pred"){
    return(cat("\nNo summary method for 'pred' object of class 'sars'\n", 
               sep = ""))
  }
  
  if (attributes(object)$type == "threshold_ci"){
    return(cat("\nNo summary method for 'threshold_ci' object of class 'sars'\n", 
               sep = ""))
  }
  if (attributes(object)$type == "threshold_coef"){
    return(cat("\nNo summary method for 'threshold_coef' object of class 'sars'\n", 
               sep = ""))
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
    res <- list("Model" = name, "residuals" = resid,
                "Parameters" = pars_tab,
                "parNames" = parN, "formula" = formula, "AIC" = round(ic, 2),
                "AICc" = round(ic2, 2), "BIC" = round(bi, 2),
                "R2" = round(R2, 2), "R2a" = round(R2a, 2),
                "observed_shape" = shape, "asymptote" = asymp,
                "convergence" = conv, "normaTest" = object$normaTest,
                "homoTest" = object$homoTest,
                "Negative_values" = object$neg_check)
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
    if(!all(object$details$ics == df$IC)){
      stop("error with ICs - contact package author")
    }
    colnames(df)[3] <- paste(cri)
    df$R2 <- vapply(object$details$fits,
                    function(x){x$R2}, FUN.VALUE = numeric(1))
    df$R2a <- vapply(object$details$fits,
                     function(x){x$R2a}, FUN.VALUE = numeric(1))
    df$Shape <- vapply(object$details$fits,
                    function(x){x$observed_shape}, FUN.VALUE = character(1))
    #a warning produces a long shape value sometimes: "observed shape 
    #algorithm failed: observed shape ..". So change these cases to sigmoid
    #shape_check <- vapply(df$Shape, FUN = function(x){grepl("failed", x)},
    #                      FUN.VALUE = logical(1))
    #if (any(shape_check)){
      #wsc <- which(shape_check)
      #df$Shape[wsc] <- "sigmoid"
   # }
    df$Asymptote <- vapply(object$details$fits,
                           function(x){x$asymptote}, FUN.VALUE = logical(1))
    df <- df[order(-df$Weight),]
    df[, 2:5] <- round(df[, 2:5], 3)
    rownames(df) <- NULL
    res <- list("Models" = Mods, "Criterion" = cri, "Model_table" = df,
                "no_fit" = nf, "Convergence" = object$details$convergence)
  }
  if (attributes(object)$type == "threshold"){
    mods <- object[[1]]
    names <- object[[2]]
    th <- object[[3]]
    #no. parameter for each model
    k <- c("ContOne" = 5, "ZslopeOne" = 4, "DiscOne" = 6, "ContTwo" = 7, 
           "ZslopeTwo" = 6, "DiscTwo" = 9, "Linear" = 3, "Intercept" = 2)
    #get IC values
    ICs <- mapply(function(x, y){
      val <- logLik(x)
      w <- which(names(k) == y)
      P <- k[w]
      n <- length(x$residuals)
      lAIC <- (2 * P) - (2 * val)
      #if denominator of AICc is 0 or negative, return Inf
      den <- n - P - 1
      if (den < 1){
        lAICc <- Inf
      } else {
        lAICc <- -2 * val + 2 * P * (n / (n - P - 1))
      }
      lBIC <- (-2 * val) + (P * log(n))
      lR2 <- summary(x)$r.squared
      lR2a <- summary(x)$adj.r.squared
      c(val, P, lAIC, lAICc, lBIC, lR2, lR2a)
    }, x = mods, y = names)
    rownames(ICs) <- c("LL", "Pars", "AIC", "AICc", "BIC", "R2", "R2a")
    colnames(ICs) <- names
    ICs <- round(t(ICs), 2)
    ICs <- as.data.frame(ICs)
    if (any(is.infinite(ICs$AICc))){
      warning("AICc not calculated for some models due to small sample size")
    }
    
    #get thresholds
    tdf <- matrix(NA, nrow = length(names), ncol = 2)
    for (i in 1:length(th)){
      if (!is.na(th[[i]][1])){
        #if (!is.primitive(object[[5]][[2]])){
        # if (object[[5]][[2]] == "none") {
        thV <- th[[i]]
        # } else {
        #  stop("something wrong with the chosen log function: A")
        #  }
        #} else if (identical(object[[5]][[2]], log10)){
        #   thV <- 10^th[[i]]
        # } else if (identical(object[[5]][[2]], log)){
        #   thV <- exp(th[[i]])
        # } else if (identical(object[[5]][[2]], log2)){
        #   thV <- 2^th[[i]]
        # } else{
        #   stop("something wrong with the chosen log function: B")
        # }
        if (length(th[[i]]) == 1){
          tdf[i, 1] <- round(thV, 3)
        } else{
          tdf[i,] <- round(thV, 3)
        }
      }
    }#eo for
    tdf <- as.data.frame(tdf)
    colnames(tdf) <- c("Th1", "Th2")
    
    #Work out the number of datapoints in each segment
    a <- object[[4]]$A
    nbi <- matrix(NA, nrow = length(names), ncol = 3)
    for (i in 1:length(th)){
      if (!is.na(th[[i]][1])){
        if (length(th[[i]]) == 1){
          nbi[i, 1:2] <- c(length(a[a <= th[[i]]]), length(a[a > th[[i]]]))
          if (sum(nbi[i, 1:2]) != length(a)){
            stop("issue with calculating no. points in each segment")
          }
        } else {
          nbi[i, 1:3] <- c(length(a[a <= th[[i]][1]]),
                           length(a[a > th[[i]][1] & a <= th[[i]][2]]),
                           length(a[a > th[[i]][2]]))
          if (sum(nbi[i, 1:3]) != length(a)){
            stop("issue with calculating no. points in each segment")
          }
        }
      }
    }
    nbi <- as.data.frame(nbi)
    colnames(nbi) <- c("seg1", "seg2", "seg3")
    #combine variables and order by 'order' IC
    mt <- cbind(ICs, tdf, nbi)
    mt <- mt[order(mt[order]),]
    res <- list("order" = order, "Model_table" = mt, "Axes transformation" = object[[5]][[1]])
  }
  class(res) <- "summary.sars"
  attr(res, "type") <- attr(object, "type")
  res
}
