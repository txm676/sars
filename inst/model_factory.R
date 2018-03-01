########################
#function to cat rOxygen comments
########################

cat_roxygen <- function(model, funName, fileName){
  
  cat1 <- function(...){cat(..., file = fileName, append = TRUE)}
  #cat1 <- function(...){cat(...)}
  
  cat1(paste0("#' ","Fit the ", model$name," model", "\n"))
  cat1("\n")
  cat1(paste0("#' @description ","Fit the ", model$name," model to SAR data", "\n"))
  cat1(paste0("#' @usage ", funName,"(data, custstart = NULL, normtest = 'lillie')", "\n"))
  cat1(paste0("#' @param ", "data ", "A dataset in the form of a dataframe with two columns: ", "\n"))
  cat1(paste0("#'   the first with island/site areas, and the second with the species richness", "\n")) 
  cat1(paste0("#'   of each island/site.", "\n"))
  cat1(paste0("#' @param ", "grid_start ", "NULL or the number of points sampled in the model parameter space", "\n"))
  cat1(paste0("#'   to run a grid search.", "\n"))
  cat1(paste0("#' @return ", "\n"))
  
  cat1(paste0("#' @examples", "\n", "#' data(galap)", "\n", "#' fit <- ", funName,
              "(galap)", "\n", "#' summary(fit)","\n", "#' plot(fit)", "\n"))
  
  cat1(paste0("#' @export"))
  
}#eo cat_roxygen

########################
#model factory function
########################
model_factory <- function(f, overwrite = FALSE){
  
  #helper function
  cat1 <- function(...){cat(..., file = file.path("R",fileName), append = T)}
  
  #sourcing the model file
  #source(file.path(system.file(package="sars"),"non_lin_models",f))
  source(file.path(getwd(),"inst","non_lin_models",f))
  
  #construct R function name 
  funName <- paste0("sar_",substr(f,5,(nchar(f) - 2)))
  
  #construct model file/function name
  #modName <- tolower(model$name)
  #funName <- paste0("sar_",modName)
  fileName <- paste0(funName,".R")
  
  #check if file already exists
  #packRoot <- system.file(package="sars")
  # packRDir <- file.path(substr(packRoot,1,(nchar(packRoot)-5)),"R")
  # filePath <- file.path(packRDir,fileName)
  # fExists <- file.exists(filePath)
  # if(fExists & !overwrite) stop("a model function '", funName, "' already exists. If you know what you do, please set overwrite=TRUE")
  #if(!(fExists) | overwrite)
  
  file.create(file.path("R",fileName), overwrite = overwrite, showWarnings = FALSE)
  #cat the roxygen2 comments
  cat_roxygen(model, funName, file.path("R",fileName))
  cat1("\n")
  cat1("\n")
  
  #function definition
  cat1(paste0(funName," <- function(data = galap, start = NULL, grid_start = NULL){","\n"))
  
  #checks
  cat1("if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe')","\n")
  cat1("if (is.matrix(data)) data <- as.data.frame(data)","\n")
  cat1("if (anyNA(data)) stop('NAs present in data')","\n")
  #cat1("normtest <- match.arg(normtest, c('none', 'shapiro', 'kolmo', 'lillie'))","\n")
  
  #data ordering and column naming (assuming Area then Species Richness)
  cat1("data <- data[order(data[,1]),]","\n")
  cat1("colnames(data) <- c('A','S')","\n")
  
  #model definition (Appending it)
  #dump("model", file = filePath, append = TRUE) #THE DUMP IS NOT PASTING THE MODEL LIST
  file.append(file.path("R",fileName), file.path(getwd(),"inst","non_lin_models",f))
                
  
  cat1("\n")
  
  cat1("model <- compmod(model)","\n")
  
  cat1("fit <- get_fit(model = model, data = data, start = start, grid_start = grid_start, algo = 'Nelder-Mead', verb = TRUE)","\n")
  cat1("if(is.na(fit$value)){","\n")
  cat1("  return(NA)","\n")
  cat1("}else{","\n")
  cat1("  obs <- obs_shape(fit)","\n")
  cat1("  fit$observed_shape <- obs$fitShape","\n")
  cat1("  fit$asymptote <- obs$asymp","\n")
  cat1("  class(fit) <- 'sars'","\n")
  cat1("  attr(fit, 'type') <- 'fit'","\n")
  cat1("  return(fit)","\n")
  cat1("}","\n")
  
  #function end
  cat1(paste0("}#end of ",funName,"\n"))
  
}#eo model_factory

#using it

modFiles <- list.files(file.path(getwd(),"inst","non_lin_models"))

lapply(modFiles, model_factory, overwrite = TRUE)





