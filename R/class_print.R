#' @importFrom stats quantile printCoefmat AIC
#' @export


print.summary.sars <- function(x, ...){
  object <- x
  if (attributes(object)$type == "lin_pow"){
    cat("Model = Log-log power\n")
    cat("\nLog-transformation function used: ", object$logT,"\n", sep = "")
    # rownames(object$Model$coefficients) <- c("LogC", "z")
    print(object$Model)
    if (object$normaTest$test == "shapiro"| object$normaTest$test == "lillie" |
        object$normaTest$test == "kolmo"){
      normP <- object$normaTest[[2]]$p.value
    } else{
      normP <- "No normality test undertaken"
    }
    #homogeneity
    if (object$homoTest$test == "cor.area" |
        object$homoTest$test == "cor.fitted"){
      homoP <- object$homoTest[[2]]$p.value
    } else{
      homoP <- "No homogeneity test undertaken"
    }
    if (is.numeric(normP) & normP < 0.05 ){
      cat("\nWarning: The normality test selected indicated the model",
          " residuals are not normally distributed (i.e. P < 0.05)\n",
          sep = "")
    }
    if (is.numeric(homoP) & homoP < 0.05){
      tr <- ifelse(object$homoTest$test == "cor.area", "area values",
                   "fitted values")
      cat("\n", paste("Warning: The homogeneity test selected indicated a",
            " signifiant correlation between the residuals and the",
            tr, "(i.e. P < 0.05)"), "\n", sep = "")
    }
    #non-linear power comparison
    if ("power" %in% names(object)){
      cat("Power (non-linear) parameters:\n",
          "c =", object$power[1], "\n",
          "z =", object$power[2], "\n")
    }
    }#eo if lin_pow

  if (attributes(object)$type == "fit"){
    cat("\nModel:\n", object$Model, "\n", sep = "")
    cat("\nCall: ","\n", as.character(object$formula), "\n", sep = "")
    cat("\nDid the model converge: ", object$convergence , "\n", sep = "")
    cat("\nResiduals:\n", sep = "")
    print(quantile(object$residuals))
    cat("\nParameters:\n", sep = "")
    mm <- object$Parameters
    #if singular gradient at parameter estimates there are no pars to print
    if (length(mm) == 1){
      cat("\nsingular gradient at parameter estimates: no parameters",
          " significance and conf. interval\n")
    } else{
       rownames(mm) <- object$parNames
       printCoefmat(mm)
    }
    cat("\nR-squared: ", object$R2 , ", Adjusted R-squared: ",
        object$R2a, "\n", sep = "")
    cat("AIC: ", object$AIC , ", AICc: ", object$AICc, ", BIC: ",
        object$BIC, "\n", sep = "")
    cat("Observed shape: ", object$observed_shape, ", Asymptote: ",
        object$asymptote, "\n", "\n", sep = "")
    #normality
    if (object$normaTest$test == "shapiro" |
        object$normaTest$test == "lillie" |
        object$normaTest$test == "kolmo"){
      normP <- object$normaTest[[2]]$p.value
    } else{
      normP <- "No normality test undertaken"
    }
    #homogeneity
    if (object$homoTest$test == "cor.area" |
        object$homoTest$test == "cor.fitted"){
      homoP <- object$homoTest[[2]]$p.value
    } else{
      homoP <- "No homogeneity test undertaken"
    }

    if (is.na(normP)){
      cat("\nWarning: The normality test returned NA, check results \n",
          sep = "")
    } else {
      if (is.numeric(normP) & normP < 0.05 ){
        cat("\nWarning: The normality test selected indicated the model",
            " residuals are not normally distributed (i.e. P < 0.05)\n",
            sep = "")
      }
    }
    if (is.na(homoP)){
      cat("\nWarning: The homogeneity test returned NA, check results \n",
          sep = "")
    } else {
      if (is.numeric(homoP) & homoP < 0.05){
        tr <- ifelse(object$homoTest$test == "cor.area", "area values",
                     "fitted values")
        cat("\n", paste("Warning: The homogeneity test selected indicated a",
                        " signifiant correlation between the residuals and the",
                        tr, "(i.e. P < 0.05)"), "\n", sep = "")
      }
    }
    #negative values check
    if (object$Negative_values){
      cat("\nWarning: The fitted values of the model contain negative",
          " values (i.e. negative species richness values)\n", sep = "")
    }
  }#eo if fit

  if (attributes(object)$type == "multi"){
    cat("\nSar_average object summary:\n", sep = "")
    cat("\n", paste(length(object$Models), " models successfully fitted"),
        "\n", sep = "")
    if (length(object$no_fit) > 1) {
      cat("\n", paste("The following models could not be fitted or were",
                      " removed due to model checks:")
                                       , "\n", sep = "")
      cat(paste(object$no_fit, collapse = ", "), "\n")
      } else if (object$no_fit == 0){
      cat("\n", paste("All models were fitted successfully"),
          "\n", sep = "")
      } else if (length(object$no_fit) == 1){
      cat("\nThe following model could not be fitted or was removed due",
          " to model checks:\n", sep = "")
      cat(object$no_fit, "\n")
    }
    cat("\n", paste("Ranked models based on", object$Criterion, " weights:"),
        "\n\n" ,sep = "")
    print(object$Model_table)
  }
  if (attributes(object)$type == "threshold"){
    cat("\nSar_threshold object summary:\n", sep = "")
    cat("\n", paste(length(object[[2]]$LL), "models fitted"),
        "\n", sep = "")
    lT <- object[[3]]
    lT2 <- switch(lT,
                  "none" = "untransformed axes",
                  "area" = "area log-transformed",
                  "both" = "area and richness log-transformed")
    cat("\n", paste("Models fitted with", lT2),
        "\n", sep = "")
    cat("\n", paste("Models ranked using", object[[1]]),
        "\n\n" ,sep = "")
    print(object[[2]])
  }
  
  if (attributes(object)$type == "habitat"){
    cat("\nsar_habitat model fit summary. Models ranked by AICc:\n\n", sep = "")
    print(object$Model_table)
  }
}


#' @export


print.sars <- function(x, ...){
  object <- x
  if (attributes(object)$type == "lin_pow"){
    cat("Model = Log-log power\n")
    cat("\nCall:\nlogS = logc + z.logA\n", sep = "")
    cat("\nLog-transformation function used: ", object$logT,"\n", sep = "")
    cat("\nCoefficients:\n", sep = "")
    logc <- object$Model$coefficients[1, 1]
    names(logc) <- "logc"
    z <- object$Model$coefficients[2, 1]
    names(z) <- "z"
    print(c(logc, z))
    cat("\n")
  }

  if (attributes(object)$type == "fit"){
    cat("\nModel:\n", object$model$name, "\n", sep = "")
    cat("\nCall:\n", as.character(object$model$formula), "\n", sep = "")
    cat("\nCoefficients:\n", sep = "")
    print(object$par)
    cat("\n")
  }

  if (attributes(object)$type == "fit_collection"){
    cat("\nThis is a fit collection\n", sep = "")
    cat("\n", length(object),
        " models contained in the fit collection\n", sep = "")
    cat( "\n", paste(unlist(lapply(object, function(x) x$model$name)),
                     collapse = ", "), "\n\n")
  }

  if (attributes(object)$type == "multi"){
    cat("\nThis is a sar_average fit object:\n", sep = "")
    cat("\n", paste(length(object$details$mod_names),
                    "models successfully fitted"), "\n", sep = "")
    if (length(object$details$no_fit) > 1) {
      cat("\n", paste(length(object$details$no_fit),
      "models were unable to be fitted or were removed due to model checks"),
          "\n", sep = "")
    } else if (object$details$no_fit != 0){
      cat("\n", paste(length(object$details$no_fit),
        "model was unable to be fitted or was removed due to model checks"),
          "\n", sep = "")
    }
     cat("\n", paste(object$details$ic, "used to rank models"),
         "\n", sep = "")
  }
  
  if (attributes(object)$type == "pred"){
    cat("\nThis is a sar_pred object:\n\n", sep = "")
    print.data.frame(object)
  }
  if (attributes(object)$type == "threshold"){
    object2 <- object
    class(object2) <- "list"
    print(object2)
  }
  if (attributes(object)$type == "threshold_ci"){
    m <- object$Method
    cat("\nThreshold confidence interval summary\n", sep = "")
    cat("\nMethod: ", m, "\n", sep = "")
    #cat("\n", "Confidence intervals generated for:", names(object[[1]]),
    #   "\n", sep = " ")
    if (m == "boot"){
      for (i in 1:length(object[[2]])){
        CI <- round(object[[2]][[i]], 2)
        cat("\nModel: ", names(object[[1]][i]),"\n", sep = "")
        cat("Confidence interval of the breakpoint:", CI[1],"-", CI[2], "\n")
      }
    } else{
      for (i in 1:(length(object) - 1)){
        CI <- round(object[[i]], 2)
        cat("\nModel: ", names(object)[i],"\n", sep = "")
        cat("Confidence interval of the breakpoint:", CI[1],"-", CI[2], "\n")
      }
    }
  }
  if (attributes(object)$type == "threshold_coef"){
    cat("\nSlopes and intercepts of fitted breakpoint models: \n\n", sep = "")
    print.data.frame(object)
  }
  
  if (attributes(object)$type == "habitat"){
    object2 <- object
    cat("\nThe individual model fit objects from sar_habitat,",
    " use summary() for the model comparison results\n\n", sep = "")
    class(object2) <- "list"
    print(object2)
  }
  
}


##internal function to calculate R2 for GDM models (same as in optim)
#' @importFrom stats fitted resid
R2_gdm <-  function(fit){
  r1 <- resid(fit)
  f1 <- fitted(fit)
  s1 <- r1 + f1 #have to do this as nls object doesn't have data in it
  rss <- sum(r1 ^ 2)
  res <- 1 - ((rss) /  sum((s1 - mean(s1))^2))
  round(res, 3)
}


#' @export
#' @importFrom AICcmodavg AICc
#' @importFrom stats AIC
print.gdm <- function(x, ...){
  
  object <- x
  if (attributes(object)$Type %in% c("loga", "linear", "power_area",
                                     "power_area_time")){
    mod <- match.arg(attributes(object)$Type, c("logarithmic",
                                                "linear", "power_area",
                                                "power_area_time"))
    if (!attributes(object)$mod_sel){
      cat("\n",paste("GDM fit using the", mod, "SAR model", sep = " "),
          "\n\n")
      object2 <- object
      class(object2) <- "nls" #need to do this as can't export :::print.nls
      print(object2)
    } else {
      cat("\n",paste("GDM fit using the", mod, "SAR model", sep = " "),
          "\n")
      cat("\nGDM model summary:\n\n")

      object2 <- object[[1]]
      print(object2)
      cat("\nAll model summaries:\n\n")
      obNL <- object[1:3]
      df <- data.frame("RSE" = vapply(obNL, function(x) summary(x)$sigma,
                                      numeric(1)),
                       "R2" = vapply(obNL, R2_gdm, numeric(1)),
                       "AIC" = vapply(obNL, AIC, numeric(1)),
                       "AICc" = vapply(obNL, AICc, numeric(1)))
      df <- rbind(df, c(summary(object[[4L]])$sigma, R2_gdm(object[[4]]),
                        AIC(object[[4]]), AICc(object[[4]])))
      df$Delta.AIC <-  df$AIC - min(df$AIC)
      df$Delta.AICc <- df$AICc - min(df$AICc)
      rownames(df) <- c("GDM", "A + T", "A", "Intercept")
      df <- df[order(df$Delta.AIC),]
      print(df)
    }
  }

    if (attributes(object)$Type == "allMods"){
      cat("\nGDM model comparison:\n\n")

      if (!attributes(object)$mod_sel){
        df <- data.frame("RSE" = vapply(object, function(x) summary(x)$sigma,
                                        numeric(1)),
                        "R2" = vapply(object, R2_gdm, numeric(1)),
                       "AIC" = vapply(object, AIC, numeric(1)),
                       "AICc" = vapply(object, AICc, numeric(1)))
      } else {
        df <- data.frame("RSE" = vapply(object,
                                        function(x) summary(x[[1]])$sigma,
                                        numeric(1)),
                         "R2" = vapply(object, function(x) R2_gdm(x[[1]]),
                                        numeric(1)),
                         "AIC" = vapply(object, function(x) AIC(x[[1]]),
                                        numeric(1)),
                         "AICc" = vapply(object, function(x) AICc(x[[1]]),
                                         numeric(1)))
      }
      df$Delta.AIC <-  df$AIC - min(df$AIC)
      df$Delta.AICc <-  df$AICc - min(df$AICc)
      rownames(df) <- c("Logarithmic", "Linear", "Power_area", "Power_area_time")
      df <- df[order(df$Delta.AIC),]
      print(df)
    }
  
  if (attributes(object)$Type == "ATT2"){
    cat("\n",paste("GDM fit using the original ATT2 model", sep = " "),
        "\n\n")
    if (!attributes(object)$mod_sel){
    object2 <- object
    class(object2) <- "lm"
    print(object2)
    } else {
      object2 <- object[[1]]
      class(object2) <- "lm"
      print(object2)
    }
    if (attributes(object)$mod_sel){
    cat("\nAll model summaries:\n\n")
    df <- data.frame("RSE" = vapply(object, function(x) summary(x)$sigma,
                                    numeric(1)),
                     "R2" = vapply(object, R2_gdm, numeric(1)),
                     "AIC" = vapply(object, AIC, numeric(1)),
                     "AICc" = vapply(object, AICc, numeric(1)))
    df$Delta.AIC <-  df$AIC - min(df$AIC)
    df$Delta.AICc <- df$AICc - min(df$AICc)
    rownames(df) <- c("GDM", "A + T", "A", "Intercept")
    df <- df[order(df$Delta.AIC),]
    print(df)
    }
  }
}