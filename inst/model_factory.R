





model_factory <- function(model,overwrite=FALSE){
  
  #construct model file/function name
  modName <- tolower(model$name)
  funName <- paste0("sar_",modName)
  fileName <- paste0(funName,".R")
  
  #check if file already exists
  packRoot <- system.file(package="sars")
  packRDir <- file.path(packRoot,"R")
  fExists <- file.exists(file.path(packRDir,fileName))
  if(fExists & !overwrite) stop("model already exists. If you know what you do, please set overwrite=TRUE")
  
  
  
  
  #cat the roxygen2 comments
  cat_roxygen(model, funName)
  
  
  
  
}#eo model_factory