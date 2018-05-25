###multi model sars


#' @note Occasionally a model fit will converge and pass the model fitting
#'   checks (e.g. residual normality) but the resulting fit is nonsensical (e.g.
#'   a horizontal line with intercept at zero). Thus, it can be useful to plot
#'   the resultant 'multi' object to check the individual model fits. To re-run
#'   the \code{sar_multi} function without a particular model, simply remove it
#'   from the \code{obj} argument.

#' @export
#' 
#' 
#' 
#' 
#' 

 #state in documentation that multicurve removes  s na RSS models. 


#models that fail normality tests removed and so may have also failed homogeneity test but not checked (as removed);
#same with negative values check.


#https://help.github.com/articles/caching-your-github-password-in-git/


sar_multi <- function(data = galap,
                       obj = c("power", "powerR","epm1","epm2","p1","p2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg", "linear"),
                       keep_details = TRUE,
                       crit = "Info",
                       normaTest = "lillie",
                       homoTest = "cor.fitted",
                       neg_check = TRUE,
                       alpha_normtest = 0.05,
                       alpha_homotest = 0.05,
                       verb = TRUE){
  
  if (!(is.character(obj))  || (class(obj) == "sars") ) stop("obj must be of class character or sars")
  
  if (is.character(obj) & is.null(data)) stop("if obj is character then data should be provided")
  
  if (is.character(obj)) {
    if (any(!(obj %in% c("linear","power","powerR","epm1","epm2","p1","p2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg")))) stop("provided model names do not match with model functions")
  }
  
  if (length(obj) < 2) stop("more than 1 fit is required to construct a sar_multi")
  
  normaTest <- match.arg(normaTest, c("none", "shapiro", "kolmo", "lillie"))
  homoTest <- match.arg(homoTest, c("none","cor.area","cor.fitted"))
  
  #if (verb) cat_line(cli::rule(left = paste0(crayon::cyan(cli::symbol$bullet),crayon::bold(" multi_sars")),right="multi-model SAR"))
  if (verb) sars:::cat_line(cli::rule(left = crayon::bold(" multi_sars"),right="multi-model SAR"))
  #if (verb) cat_line(crayon::magenta(cli::symbol$arrow_right)," Data set is: ")
  #if (verb) cat_line(cli::rule(left = paste0(crayon::magenta(cli::symbol$bullet))))
  #if (verb) bullet("O | S : model", bullet = blue_arrow())
  
  #if not yet fitted, fit the models to the datasquare_small_filled
  if (is.character(obj)) {
   
    mods <- paste0("sar_",obj)
    names(mods) <- obj
    
    fits <- suppressWarnings(lapply(obj, function(x){
      
      f <- eval(parse(text = paste0(mods[x],"(data", ", normaTest = ", paste0("'", normaTest, "'"), ", homoTest = ", paste0("'", homoTest, "'"), ")")))
      
      if (verb) {
        if(is.na(f$value)) {
          sars:::cat_line( paste0(crayon::red(cli::symbol$arrow_right)," ",crayon::col_align(x,max(nchar(obj)))," : ", crayon::red(cli::symbol$cross)))
        }else{
          
          if (!is.matrix(f$sigConf)){
            sars::: cat_line( paste0(crayon::yellow(cli::symbol$arrow_right)," ",crayon::col_align(x,max(nchar(obj)))," : Warning: could not compute parameters statistics"))
          }else{
            sars:::cat_line( paste0(crayon::cyan(cli::symbol$arrow_right)," ",crayon::col_align(x,max(nchar(obj)))," : ",crayon::green(cli::symbol$tick)))
          }
        }
      }
      
      f
      
    }))#eo suppressWarnings(lapply)
  
    
    #remove models with no parameter estimates
  #  sigC <- vapply(fits, function(x) any(is.na(x$sigConf)), FUN.VALUE = logical(1))
    #if(all(sigC)) stop("No model could be fitted, aborting multi_sars\n")
  
   # if (any(sigC)){
    #  warning("Could not compute parameter statistics for one or more models and these have been excluded from the multi SAR", call. = FALSE)
    #  badNames2 <- vapply(fits[sigC], FUN = function(x){x$model$name}, FUN.VALUE = character(1))
    #  if (badMods != 0) {
    #    badMods <- c(badMods, badNames2)
    #  } else{
     #   badMods <- badNames2
      #}
      #fits <- fits[!sigC]
   #}
  
 
    
  }else{
    if (attributes(obj)$type == "fit_collection"){
      fits <- obj
      if(all(is.na(fits))){
        stop("The fit collection had no fits, aborting multi_sars\n")
      }
    }else{
      stop("an object of class 'sars' is passed but not of type 'fit_collection':(")
    }
  }#eo if else is.character(obj)
  
  #####BAD MODEL CHECKS#######################
  
  #NA CHECKS
  f_nas <- unlist(lapply(fits,function(b)b$value))
  
  if(all(is.na(f_nas))){
    stop("No model could be fitted, aborting multi_sars\n")
  }
  
  badMods <- vector(length = 0, mode = "character")
  
  if(any(is.na(f_nas))){
    badNames <- is.na(f_nas)
    warning(paste(length(which(is.na(f_nas))), "models could not be fitted and have been excluded from the multi SAR"), call. = FALSE)
    badMods <- obj[badNames] #extract the bad model names from the obj vector (not from fits, as no model name if NA)
    fits <- fits[!is.na(f_nas)]
  }
  
  #if checks for normality and / or homoscedasticity enabled, then check and remove bad fits from fits

  if (normaTest != "none") {
    np <- vapply(fits, function(x) x$normaTest[[2]]$p.value, FUN.VALUE = numeric(1))
    whp <- np < alpha_normtest
    if (any(whp)) {
      warning(paste(length(which(np < alpha_normtest)), "models failed the residuals normality test and
              have been excluded from the multi SAR"), call. = FALSE)
      #get model names
      mn <- vapply(fits, function(x) x$model$name, FUN.VALUE = character(1))#get all names in fit collection
      badMods <- c(badMods, mn[whp])#select the model names for models with p < 0.05
      fits <- fits[!whp]#then remove these models from the fit collection
    }
  }
  if (homoTest != "none") {
    hp <- vapply(fits, function(x) x$homoTest[[2]]$p.value, FUN.VALUE = numeric(1))
    whh <- hp < alpha_homotest
    if (any(whh)) {
      warning(paste(length(which(hp < alpha_homotest)),"models failed the residuals homogeneity test and have been excluded from the multi SAR"), call. = FALSE)
      mn <- vapply(fits, function(x) x$model$name, FUN.VALUE = character(1))#get all names in fit collection
      badMods <- c(badMods, mn[whh])#select the model names for models with p < 0.05
      fits <- fits[!whh]
    }
  }
  #negative values
  if (neg_check){
    nc <- vapply(fits, function(x) any(x$calculated < 0), FUN.VALUE = logical(1))
    if (any(nc)) {
      warning(paste(length(which(nc)), "models have negative fitted values and have been excluded from the multi SAR"), call. = FALSE)
      mn <- vapply(fits, function(x) x$model$name, FUN.VALUE = character(1))#get all names in fit collection
      badMods <- c(badMods, mn[nc])#select the model names for models with p < 0.05
      fits <- fits[!nc]
    }
  }

  if (length(badMods) == 0) badMods <- 0
  if (length(fits) < 2) stop("Fewer than two models could be fitted and / or passed the model checks")
  
  
  fits <- fit_collection(fits = fits)
  
####################################
  
  #setting variables
  nPoints <- length(fits[[1]]$data$A)
  nMods <- length(fits)

  modNames <- vapply(fits, FUN = function(x){x$model$name}, FUN.VALUE = character(1))
  if(is.character(obj)){  keep_details <- TRUE }

  
  #choosing an IC criterion (AIC or AICc or BIC)
  IC <- switch(crit,
         Info= if ( (nPoints / 3) < 40 ) { "AICc" } else { "AIC" },
         Bayes= "BIC"
  )
  
  #get ICs
  ICs <- vapply(X = fits, FUN = function(x){x[[IC]]}, FUN.VALUE = double(1))
  
  #get delta ICs
  delta_ICs <- ICs - min(ICs)
  
  #get akaike weights
  akaikesum <- sum(exp( -0.5*(delta_ICs)))
  weights_ICs <- exp(-0.5*delta_ICs) / akaikesum
  
  #ERROR: produce weight averaged diversity measures
 # mmi <- vapply(fits,FUN=function(x){x$calculated},FUN.VALUE=double(nPoints))
 # mmi <- apply((mmi * weights_ICs), 1 , sum)
  
  #produce weight averaged diversity measures
  mmi <- vapply(fits,FUN=function(x){x$calculated},FUN.VALUE=double(nPoints))
  wm <- matrix(nrow = nPoints, ncol = length(fits))
  for (i in seq_along(weights_ICs)) {wm[ ,i] <- mmi[ ,i] * weights_ICs[i]}
  mmi <- apply(wm, 1 , sum)

  res <- mmi
  
  if(keep_details){
    
    details <- list(
      mod_names = modNames,
      fits = fits,
      crit = crit,
      ic = IC,
      norm_test = normaTest,
      homo_test = homoTest,
      alpha_norm_test = alpha_normtest,
      alpha_homo_test = alpha_homotest,
      ics = ICs,
      delta_ics = delta_ICs,
      weights_ics = weights_ICs,
      n_points = nPoints,
      n_mods = nMods,
      no_fit = as.vector(badMods)
    )
    
    res <- list(mmi = mmi, details = details)
  }#eo if keep_details 
  
  class(res) <- c("multi", "sars")
  attr(res, "type") <- "multi"
  
  #if (verb) cat_line(cli::rule(left = crayon::cyan(cli::symbol$bullet)))
  if (verb) sars:::cat_line(cli::rule())
  
  invisible(res)
  
}#end of multi_sars

