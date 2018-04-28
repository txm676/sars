###multi model sars

#' @export

multi_sars <- function(data = galap,
                       obj = paste0("sar_",c("power", "powerR","epm1","epm2","p1","p2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg")),
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
  
  if (verb) cat_line(cli::rule(left = paste0(crayon::cyan(cli::symbol$bullet)," multi_sars: multi-model SAR")))
  #if (verb) bullet("O | S : model", bullet = blue_arrow())
  
  #if not yet fitted, fit the models to the datasquare_small_filled
  if (is.character(obj)) {
   
    fits <- suppressWarnings(lapply(obj, function(x){
      
      f <- eval(parse(text = paste0(x,"(data)")))
      
      if (verb) {
        if(is.na(f$value)) {
          cat_line(crayon::cyan(cli::symbol$arrow_right)," ",x," : ",crayon::red(cli::symbol$cross))
        }else{
          
          if (!is.matrix(f$sigConf)){
            cat_line(crayon::cyan(cli::symbol$arrow_right)," ",x," : ",crayon::yellow(cli::symbol$circle_filled)," | warning: could not compute parameters statistics")
          }else{
            cat_line(crayon::cyan(cli::symbol$arrow_right)," ",x," : ",crayon::green(cli::symbol$tick))
          }
        }
      }
      
      f
      
    }))#eo suppressWarnings(lapply)
    
    f_nas <- unlist(lapply(fits,function(b)b$value))
    
    if(all(is.na(f_nas))){
      stop("No model could be fitted, aborting multi_sars\n")
    }
    
    if(any(is.na(f_nas))){
      warning(warned()," One or more models could not be fitted and have been excluded from the multi SAR", call. = FALSE)
      fits <- fits[!is.na(f_nas)]
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
  
  class(res) <- c("multi.sars", "sars")
  attr(res, "type") <- "multi"
  
  if (verb) cat_line(cli::rule(left = crayon::cyan(cli::symbol$bullet)))
  
  invisible(res)
  
}#end of multi_sars

