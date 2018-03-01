###multi model sars

#' @export

multi_sars <- function(obj = paste0("sar_",c("power", "powerR","epm1","epm2","p1","p2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg")),
                       data = galap,
                       keep_details = TRUE,
                       crit = "Info",
                       normtest = "lillie",
                       homotest = "cor.fitted",
                       alpha_normtest = 0.05,
                       alpha_homotest = alpha_normtest,
                       verb = TRUE){
  
  if (!(is.character(obj))  || (class(obj) == "sars") ) stop("obj must be of class character or sars")
  
  if (is.character(obj) & is.null(data)) stop("if obj is character then data should be provided")
  
  if (is.character(obj)) {
    if (any(!(obj %in% paste0("sar_",c("linear","power","powerR","epm1","epm2","p1","p2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg"))))) stop("provided model names do not match with model functions")
  }
  
  if (length(obj) < 2) stop("more than 1 fit is required to construct a multi_sar")
  
  normtest <- match.arg(normtest, c("none", "shapiro", "kolmo", "lillie"))
  homotest <- match.arg(homotest, c("none","cor.area","cor.fitted"))
  
  #if not yet fitted, fit the models to the data
  if (is.character(obj)) {
   
    fits <- lapply(obj, function(x){
      
      if (verb) cat("-- fitting model: ", x, "\n")
      
      eval(parse(text = paste0(x,"(data)")))
    })
    
    if(all(is.na(fits))){
      stop("No model could be fitted, aborting multi_sars\n")
    }
    
    if(any(is.na(fits))){
      warning("One or more models could not be fitted and have been excluded from the multi SAR\n")
      fits <- fits[!is.na(fits)]
    }
    
    fits <- fit_collection(fits = fits)
    
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
  
  
  #if checks for normality and / or homoscedasticity enabled, then check and remove bad fits from fits
  #if length(fits) < 2 -> stop
  
  #setting variables
  nPoints <- length(fits[[1]]$data$A)
  nMods <- length(fits)
  modNames <- vapply(fits, FUN = function(x){x$model$name}, FUN.VALUE = character(1))
  if(is.character(obj)){  keep_fits <- TRUE }
  
  #choosing an IC criterion (AIC or AICc or BIC)
  IC <- switch(crit,
         Info= if ( (nPoints / 3) < 40 ) { "AICc" } else { "AIC" },
         Bayes= "BIC"
  )
  
  #get ICs
  ICs <- vapply(X = fits, FUN = function(x){x[[IC]]}, FUN.VALUE = double(1))
  
  #get delta ICs
  delta_ICs <- ICs <- min(ICs)
  
  #get akaike weights
  akaikesum <- sum(exp( -0.5*(delta_ICs)))
  weights_ICs <- exp(-0.5*delta_ICs) / akaikesum
  
  #produce weight averaged diversity measures
  mmi <- vapply(fits,FUN=function(x){x$calculated},FUN.VALUE=double(nPoints))
  mmi <- apply((mmi * weights_ICs), 1 , sum)
  
  res <- mmi
  
  if(keep_details){
    
    details <- list(
      mod_names = modNames,
      fits = fits,
      crit = crit,
      ic = IC,
      norm_test = normtest,
      homo_test = homotest,
      alpha_norm_test = alpha_normtest,
      alpha_homo_test = alpha_homotest,
      ics = ICs,
      delta_ics = delta_ICs,
      weights_ics = weights_ICs,
      n_points = nPoints,
      n_mods = nMods
    )
    
    res <- list(mmi = mmi, details = details)
  }#eo if keep_details 
  
  class(res) <- "sars"
  attr(res, "type") <- "multi_sars"
  
  invisible(res)
  
}#end of multi_sars

