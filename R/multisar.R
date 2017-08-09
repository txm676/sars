###empty multisar for now

multi_sar <- function(obj,data=NULL){
  
  if (!(is.character(obj))  || (attributes(m)$type != "fitcollection") ) stop("obj must be a character or fitcollection")
  
  if( is.character(obj) & is.null(data)) stop("if obj is character then data should be provided")
  
  if(is.character(obj)){
    if(any(!(obj %in% paste0("sar_",c("linear","power","power_R","epm1","epm2","P1","P2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg"))))) stop("all model names should be ok")
  }
  
  #if not yet fitted, fit the models to the data
  if(is.character(obj)){
   
    fits <- lapply(obj,function(x){
      eval(parse(text=paste0(x,"(data)")))
    })
    
    
  }
  
  
  
  
}#end of multisar