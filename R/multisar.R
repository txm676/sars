###empty multisar for now

#' @export

multi_sar <- function(obj,data=NULL,keep_fits=FALSE,crit="Info",normtest="lillie",homotest="cor.fitted",alpha_normtest=0.05,alpha_homotest=alpha_normtest){
  
  if (!(is.character(obj))  || (attributes(m)$type != "fitcollection") ) stop("obj must be a character or fitcollection")
  
  if( is.character(obj) & is.null(data)) stop("if obj is character then data should be provided")
  
  if(is.character(obj)){
    if(any(!(obj %in% paste0("sar_",c("linear","power","power_R","epm1","epm2","P1","P2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg"))))) stop("all model names should be ok")
  }
  
  normtest <- match.arg(normtest, c("none", "shapiro", "kolmo", "lillie"))
  homotest <- match.arg(homotest, c("none","cor.area","cor.fitted"))
  
  #if not yet fitted, fit the models to the data
  if(is.character(obj)){
   
    fits <- lapply(obj,function(x){
      eval(parse(text=paste0(x,"(data)")))
    })
    
    fits <- fit_collection(fits)
    
  }
  
  #if checks for normality and / or homoscedasticity enabled, then check and remove bad fits from fits
  
  #setting variables
  nPoints <- length(fits[[1]]$data$A)
  nMods <- length(fits)
  modNames <- vapply(fits, FUN = function(x){x$model$name}, FUN.VALUE = character(1))
  if(is.character(obj)){  keep_fits <- TRUE }
  
  #choosing an IC criterion (AIC or AICc or BIC)
  IC <- switch(crit,
         Info= if ( (nPoints / 3) < 40 ) { "AICc" } else { "AIC"},
         Bayes= "BIC"
  )
  
  #get delta ICs
  delta_ICS <- ICs <- min(ICs)
  
  #get akaike weights
  akaikesum <- sum(exp( -0.5*(delta_ICs)))
  weights_ICs <- exp(-0.5*delta_ICs) / akaikesum
  
  #produce weight averaged diversity measures
  mmS <- vapply(fits,FUN=function(x){x$calculated},FUN.VALUE=double(nPoints))
  mmS <- apply((mmS * weights_ICs), 1 , sum)
  
  res <- mmS
  
  if(keep_fits) res$fits <- as.list(fits)
  
}#end of multisar