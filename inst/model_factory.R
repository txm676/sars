########################
#function to cat rOxygen comments
########################

cat_roxygen <- function(model, funName, fileName){
  
  cat1 <- function(...){cat(..., file = fileName, append = TRUE)}
  #cat1 <- function(...){cat(...)}
  
  cat1(paste0("#' ","Fit the ", model$name," model", "\n"))
  cat1("\n")
  cat1(paste0("#' @description ","Fit the ", model$name," model to SAR data.", "\n"))
  cat1(paste0("#' @usage ", funName,"(data, start = NULL, grid_start = NULL, normaTest =  'lillie',", "\n"))
  cat1(paste0("#'   homoTest = 'cor.fitted')", "\n"))
  cat1(paste0("#' @param ", "data ", "A dataset in the form of a dataframe with two columns: ", "\n"))
  cat1(paste0("#'   the first with island/site areas, and the second with the species richness", "\n")) 
  cat1(paste0("#'   of each island/site.", "\n"))
  cat1(paste0("#' @param ", "start ", "NULL or custom parameter start values for the optimisation algorithm.", "\n"))
  cat1(paste0("#' @param ", "grid_start ", "NULL or the number of points sampled in the model parameter space", "\n"))
  cat1(paste0("#'   or FALSE to prevent any grid start after a fail in inital optimization", "\n"))
  cat1(paste0("#'   to run a grid search.", "\n"))
  cat1(paste0("#' @param ", "normaTest ", "The test used to test the normality of the residuals of the", "\n"))
  cat1(paste0("#'   model. Can be any of 'lillie' (Lilliefors Kolmogorov-Smirnov test; the", "\n"))
  cat1(paste0("#'   default), 'shapiro' (Shapiro-Wilk test of normality), 'kolmo'", "\n"))
  cat1(paste0("#'   (Kolmogorov-Smirnov test), or 'none' (no residuals normality test is undertaken).", "\n"))
  cat1(paste0("#' @param ", "homoTest ","The test used to check for homogeneity of the residuals of", "\n"))
  cat1(paste0("#'   the model. Can be any of 'cor.fitted' (a correlation of the residuals with", "\n"))
  cat1(paste0("#'   the model fitted values; the default), 'cor.area' (a correlation of the", "\n"))
  cat1(paste0("#'   residuals with the area values), or 'none' (no residuals homogeneity test is undertaken).","\n"))
  
  cat1(paste0("#' @details The model is fitted using non-linear regression. The model parameters are estimated", "\n"))
  cat1(paste0("#'   by minimizing the residual sum of squares with an unconstrained Nelder-Mead optimization algorithm", "\n")) 
  cat1(paste0("#'   and the \\code{\\link{optim}} function. To avoid numerical problems and speed up the convergence process,", "\n")) 
  cat1(paste0("#'   the starting values used to run the optimization algorithm are carefully chosen, or custom values can be provided", "\n")) 
  cat1(paste0("#'   using the argument \\code{start}. The fitting process also determines the observed shape of the model fit,", "\n")) 
  cat1(paste0("#'   and whether or not the observed fit is asymptotic (see Triantis et al. 2012 for further details).", "\n", "\n"))
  cat1(paste0("#'   Model validation is undertaken by assessing the normality (\\code{normaTest}) and homogeneity (\\code{homoTest})", "\n")) 
  cat1(paste0("#'   of the residuals and a warning is provided in \\code{\\link{summary.sars}} if either test is failed.", "\n", "\n")) 
  cat1(paste0("#'   A selection of information criteria (e.g. AIC, BIC) are returned and can be used to compare models", "\n")) 
  cat1(paste0("#'   (see also \\code{\\link{fit_collection}} and \\code{\\link{sar_multi}}).", "\n")) 
  
  cat1(paste0("#' @return ", "A list of class 'sars' with the following components: ", "\n"))
  cat1(paste0("#'   \\itemize{", "\n"))
  cat1(paste0("#'     \\item{par} { The model parameters}", "\n"))
  cat1(paste0("#'     \\item{value} { Residual sum of squares}", "\n"))
  cat1(paste0("#'     \\item{counts} {  The number of iterations for the convergence of the fitting algorithm}", "\n"))
  cat1(paste0("#'     \\item{convergence} { Numeric code indicating model convergence (0 = converged)}", "\n"))
  cat1(paste0("#'     \\item{message} { Any message from the model fit algorithm}", "\n"))
  cat1(paste0("#'     \\item{hessian} { A symmetric matrix giving an estimate of the Hessian at the solution found}", "\n"))
  cat1(paste0("#'     \\item{verge} { Logical code indicating model convergence}", "\n"))
  cat1(paste0("#'     \\item{startValues} { The start values for the model parameters used in the optimisation}", "\n"))
  cat1(paste0("#'     \\item{data} { Observed data}", "\n"))
  cat1(paste0("#'     \\item{model} { A list of model information (e.g. the model name and formula)}", "\n"))
  cat1(paste0("#'     \\item{calculated} {  The fitted values of the model}", "\n"))
  cat1(paste0("#'     \\item{residuals} { The model residuals}", "\n"))
  cat1(paste0("#'     \\item{AIC} { The AIC value of the model}", "\n"))
  cat1(paste0("#'     \\item{AICc} { The AICc value of the model}", "\n"))
  cat1(paste0("#'     \\item{BIC} { The BIC value of the model}", "\n"))
  cat1(paste0("#'     \\item{R2} { The R2 value of the model}", "\n"))
  cat1(paste0("#'     \\item{R2a} { The adjusted R2 value of the model}", "\n"))
  cat1(paste0("#'     \\item{sigConf} { The model coefficients table}", "\n"))
  cat1(paste0("#'     \\item{normaTest} { The results of the residuals normality test}", "\n"))
  cat1(paste0("#'     \\item{homoTest} { The results of the residuals homogeneity test}", "\n"))
  cat1(paste0("#'     \\item{observed_shape} { The observed shape of the model fit}", "\n"))
  cat1(paste0("#'     \\item{asymptote} { A logical value indicating whether the observed fit is asymptotic}}", "\n", "\n"))
  cat1(paste0("#'   The \\code{\\link{summary.sars}} function returns a more useful summary of", "\n"))
  cat1(paste0("#'   the model fit results, and the \\code{\\link{plot.sars}} plots the model fit.", "\n"))
  
  cat1(paste0("#' @references Triantis, K.A., Guilhaumon, F. & Whittaker, R.J. (2012) The island species-area", "\n")) 
  cat1(paste0("#'   relationship: biology and statistics. Journal of Biogeography, 39, 215-231.", "\n"))
  
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
  cat1(paste0(funName,' <- function(data = galap, start = NULL, grid_start = NULL, normaTest =  "lillie",
              homoTest = "cor.fitted"){',"\n"))
  
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
  
  cat1("fit <- get_fit(model = model, data = data, start = start, grid_start = grid_start, algo = 'Nelder-Mead', 
       normaTest =  normaTest, homoTest = homoTest, verb = TRUE)","\n")
  cat1("if(is.na(fit$value)){","\n")
  cat1("  return(list(value = NA))","\n")
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





