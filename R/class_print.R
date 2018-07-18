
#' @export


print.summary.sars <- function(object){
  
  if (attributes(object)$type == "lin_pow"){
    cat("Model = ","Log-log power", "\n")
   # rownames(object$Model$coefficients) <- c("LogC", "z")
    print(object$Model)   
    if (object$Normality_test$test == "shapiro"|| object$Normality_test$test == "lillie" ||
        object$Normality_test$test == "kolmo"){
      normP <- object$Normality_test[[2]]$p.value
    } else{
      normP <- "No normality test undertaken"
    }
    #homogeneity
    if (object$Homogeneity_test$test == "cor.area" || object$Homogeneity_test$test == "cor.fitted"){
      homoP <- object$Homogeneity_test[[2]]$p.value
    } else{
      homoP <- "No homogeneity test undertaken"
    }
    if (is.numeric(normP) && normP < 0.05 ){
      cat("\n", "Warning: The normality test selected indicated the model residuals are
          not normally distributed (i.e. P < 0.05)", "\n", sep = "")
    }
    if (is.numeric(homoP) && homoP < 0.05){
      tr <- ifelse(object$Homogeneity_test$test == "cor.area", "area values", "fitted values")
      cat("\n", paste("Warning: The homogeneity test selected indicated a signficant correlation
                      between the residuals and the",tr, "(i.e. P < 0.05)"), "\n", sep = "")
    }
    #non-linear power comparison
    if ("power" %in% names(object)){
      cat("Power (non-linear) parameters:", "\n", 
          "c =", object$power[1], "\n",
          "z =", object$power[2])
    }
  }#eo if lin_pow
  
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
    #normality
    if (object$Normality_test$test == "shapiro"|| object$Normality_test$test == "lillie" ||
        object$Normality_test$test == "kolmo"){
      normP <- object$Normality_test[[2]]$p.value
    } else{
      normP <- "No normality test undertaken"
    }
    #homogeneity
    if (object$Homogeneity_test$test == "cor.area" || object$Homogeneity_test$test == "cor.fitted"){
      homoP <- object$Homogeneity_test[[2]]$p.value
    } else{
      homoP <- "No homogeneity test undertaken"
    }

    if (is.numeric(normP) && normP < 0.05 ){
      cat("\n", "Warning: The normality test selected indicated the model residuals are
          not normally distributed (i.e. P < 0.05)", "\n", sep = "")
    }
    if (is.numeric(homoP) && homoP < 0.05){
      tr <- ifelse(object$Homogeneity_test$test == "cor.area", "area values", "fitted values")
      cat("\n", paste("Warning: The homogeneity test selected indicated a signficant correlation
          between the residuals and the",tr, "(i.e. P < 0.05)"), "\n", sep = "")
    }
    #negative values check
    if (object$Negative_values == 1){
      cat("\n", "Warning: The fitted values of the model contain negative values (i.e. negative 
          species richness values)", "\n", sep = "")
      }
  }#eo if fit

  if (attributes(object)$type == "multi"){ 
    cat("\n", "Sar_multi object summary: ", "\n", sep = "")
    cat("\n", paste(length(object$Models), " models successfully fitted"), "\n", sep = "")
    if (length(object$no_fit) > 1) {
      cat("\n", paste("The following models could not be fitted or were removed due to model checks:")
                                       , "\n", sep = "")
      cat(paste(object$no_fit, collapse = ", "), "\n")
      } else if (object$no_fit == 0){
      cat("\n", paste("All models were fitted successfully"), 
          "\n", sep = "")
      } else if (length(object$no_fit) == 1){
      cat("\n", "The following model could not be fitted or was removed due to model checks:", "\n", sep = "")
      cat(object$no_fit, "\n")
    }
    cat("\n", paste("Ranked models based on", object$Criterion, " weights:"), "\n", "\n" ,sep = "")
    print(object$Model_table)
  }
}


#' @export
#' 

print.sars <- function(object){
  
  if (attributes(object)$type == "lin_pow"){
    cat("Model = ","Log-log power", "\n")
    cat("\n", "Call: ","\n", "logS = logc + z.logA", "\n", sep = "")
    cat("\n", "Coefficients: ", "\n", sep = "")
    logc <- object$Model$coefficients[1, 1]
    names(logc) <- "logc"
    z <- object$Model$coefficients[2, 1]
    names(z) <- "z"
    print(c(logc, z))
    cat("\n")
  }
  
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
    cat( "\n", unlist(lapply(object, function(x) x$model$name)), "\n", "\n")
  }
  
  if (attributes(object)$type == "multi"){ 
    cat("\n", "This is a sar_multi fit object:", "\n", sep = "")
    cat("\n", paste(length(object$details$mod_names), "models successfully fitted"), "\n", sep = "")
    if (length(object$details$no_fit) > 1) {
      cat("\n", paste(length(object$details$no_fit), "models were unable to be fitted or were removed due to model checks"), 
          "\n", sep = "")
    } else if (object$details$no_fit != 0){
      cat("\n", paste(length(object$details$no_fit), "model was unable to be fitted or was removed due to model checks"),
          "\n", sep = "")
    }
     cat("\n", paste(object$details$ic, "used to rank models"), "\n", sep = "")
  }
}

#' @export

print.gdm <- function(object){
  if (attributes(object)$Type == "lin_pow"){
    if (!attributes(object)$mod_sel){
      cat("\n","GDM fit using the linear (log-log) power SAR model", "\n")
      stats:::print.lm(object)
    } else {
      cat("\n","GDM fit using the linear (log-log) power SAR model", "\n")
      cat("\n","GDM model summary:", "\n")
      stats:::print.lm(object[[1]])
      cat("\n","All model summaries:", "\n", "\n")
      
      df <- data.frame("R2" = vapply(object, function(x) summary(x)$r.squared, numeric(1)),
                       "Adj.R2" = vapply(object, function(x) summary(x)$adj.r.squared, numeric(1)),
                       "AIC" = vapply(object, AIC, numeric(1)))
      df$Delta.AIC <- df$AIC - min(df$AIC)
      rownames(df) <- c("GDM", "A + T", "A", "Intercept")
      df <- df[order(df$Delta.AIC),]
      print(df)
    }
  }
  
  if (attributes(object)$Type %in% c("expo", "linear", "power")){
    mod <- match.arg(attributes(object)$Type, c("exponential", "linear", "power"))
    if (!attributes(object)$mod_sel){
      cat("\n",paste("GDM fit using the", mod, "SAR model", sep = " "),
          "\n", "\n")
      stats:::print.nls(object)
    } else {
      cat("\n",paste("GDM fit using the", mod, "SAR model", sep = " "),
          "\n")
      cat("\n","GDM model summary:", "\n", "\n")
      stats:::print.nls(object[[1]])
      cat("\n","All model summaries:", "\n", "\n")
      
      df <- data.frame("RSE" = vapply(object, function(x) summary(x)$sigma, numeric(1)),
                       "AIC" = vapply(object, AIC, numeric(1)))
      df$Delta.AIC <-  df$AIC - min(df$AIC)
      rownames(df) <- c("GDM", "A + T", "A")
      df <- df[order(df$Delta.AIC),]
      print(df)
    }
  }
    
    if (attributes(object)$Type == "allMods"){
      cat("\n","GDM model comparison:", "\n", "\n")
      
      if (!attributes(object)$mod_sel){
        df <- data.frame("RSE" = vapply(object, function(x) summary(x)$sigma, numeric(1)),
                       "AIC" = vapply(object, AIC, numeric(1)))
      } else {
        df <- data.frame("RSE" = vapply(object, function(x) summary(x[[1]])$sigma, numeric(1)),
                         "AIC" = vapply(object, function(x) AIC(x[[1]]), numeric(1)))
      }
      df$Delta.AIC <-  df$AIC - min(df$AIC)
      rownames(df) <- c("Exponential", "Linear")
      df <- df[order(df$Delta.AIC),]
      print(df)
    }
  }
  
  
  

      
  
  
