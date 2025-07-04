################################################################################
#
#              MODEL FACTORY
#
# This file contains the "model factory" machinery of the sars R package
# The function "model_factory" is applied over all the files present in
# the "inst/non_lin_models" directory. Each of these file define a particular
# non linear SAR model (see files in the "inst/non_lin_models" for examples).
# The function "model_factory" creates the individual non linear models .R files 
# present in the "R" directory, providing a mechanism for expanding the number
# of non linear SAR models included in the R package:
#
# 1. clone the "https://github.com/txm676/sars" repository
# 2. create a "model file" containing the (list) description of the model and
#    save it in the "inst/non_lin_models" directory
# 3. source the "cat_roxygen" and "model_factory" functions
# 4. run the "model_factory" fucntion (see commands at the end of this file)
# 5. build the sars package
# 
################################################################################


########################
#function to cat rOxygen comments
########################

cat_roxygen <- function(model, funName, fileName){
  
  cat1 <- function(...){cat(..., file = fileName, append = TRUE)}
  #cat1 <- function(...){cat(...)}
  
  cat1(paste0("#' ","Fit the ", model$name," model", "\n"))
  cat1("\n")
  if (funName == "sar_mmf"){
    cat1(paste0("#' @description ","Fit the ", model$name,
                " model to SAR data. This function has been deprecated.", "\n"))
  } else {
    cat1(paste0("#' @description ","Fit the ", model$name,
              " model to SAR data.", "\n"))
  }
  cat1(paste0("#' @usage ", funName, "(data, start = NULL,",
              " grid_start = 'partial',", "\n"))
  cat1(paste0("#'   grid_n = NULL, normaTest = 'none',", "\n"))
  cat1(paste0("#'   homoTest = 'none', homoCor = 'spearman', verb = TRUE)\n"))
  cat1(paste0("#' @param ", "data ", "A dataset in the form of a dataframe ",
                "with two columns: \n"))
  cat1(paste0("#'   the first with island/site areas, and the second with ",
                  "the species richness\n")) 
  cat1(paste0("#'   of each island/site.\n"))
  cat1(paste0("#' @param ", "start ", "NULL or custom parameter start", 
              " values for the optimisation algorithm.\n"))
  cat1(paste0("#' @param ", "grid_start ", "Should a grid search procedure be", 
              " implemented to test multiple starting parameter values.",
              " Can be one of 'none', 'partial' or 'exhaustive' The default",
              " is set to 'partial'.\n"))
  cat1(paste0("#' @param ", "grid_n ", "If \\code{grid_start = exhaustive}, the",
              " number of points sampled in the starting parameter space.\n"))
  cat1(paste0("#' @param ", "normaTest ", "The test used to test the", 
              " normality of the residuals of the\n"))
  cat1(paste0("#'   model. Can be any of 'lillie' (Lilliefors ", 
              "test\n"))
  cat1(paste0("#', 'shapiro' (Shapiro-Wilk test of normality),", 
              " 'kolmo'", "\n"))
  cat1(paste0("#'   (Kolmogorov-Smirnov test), or 'none' (no residuals ", 
              "normality test is undertaken; the default).\n"))
  cat1(paste0("#' @param ", "homoTest ","The test used to check for", 
              " homogeneity of the residuals of\n"))
  cat1(paste0("#'   the model. Can be any of 'cor.fitted' (a correlation ", 
              "of the residuals with\n"))
  cat1(paste0("#'   the model fitted values), 'cor.area'", 
              " (a correlation of the\n"))
  cat1(paste0("#'   residuals with the area values), or 'none' (no residuals", 
              " homogeneity test is undertaken; the default).\n"))
  cat1(paste0("#' @param ", "homoCor ","The correlation test to be used", 
              " when \\code{homoTest !='none'}. Can be any of 'spearman'\n"))
  cat1(paste0("#'   (the default), 'pearson', or 'kendall'.\n"))
  cat1(paste0("#' @param ", "verb ","Whether or not to print certain warnings ", 
              "(default = TRUE)\n"))
  
  cat1(paste0("#' @details The model is fitted using non-linear regression.", 
              " The model parameters are estimated", "\n"))
  cat1(paste0("#'   by minimizing the residual sum of squares with an", 
              " unconstrained Nelder-Mead optimization algorithm\n")) 
  cat1(paste0("#'   and the \\code{\\link{optim}} function. To avoid", 
              " numerical problems and speed up the convergence process,\n")) 
  cat1(paste0("#'   the starting values used to run the optimization", 
              " algorithm are carefully chosen. However, if this does\n")) 
  cat1(paste0("#' not work, custom values can be provided (using the", 
              " \\code{start} argument), or a more comprehensive search\n"))
  cat1(paste0("#'   can be undertaken using the \\code{grid_start} argument.",
              " See the vignette for more information.\n")) 
  cat1(paste0("#'   The fitting process", 
              " also determines the observed shape of the model fit,\n")) 
  cat1(paste0("#'   and whether or not the observed fit is asymptotic (see", 
              " Triantis et al. 2012 for further details).\n\n"))
  cat1(paste0("#'   Model validation can be undertaken by assessing the", 
              " normality (\\code{normaTest}) and homogeneity", 
              " (\\code{homoTest})\n")) 
  cat1(paste0("#'   of the residuals and a warning is provided in", 
              " \\code{\\link{summary.sars}} if either test is chosen and fails.",
              "\n\n")) 
  cat1(paste0("#'   A selection of information criteria (e.g. AIC, BIC) are", 
              " returned and can be used to compare models\n")) 
  cat1(paste0("#'   (see also \\code{\\link{sar_average}}).\n")) 
  cat1(paste0("#'   \n")) 
  cat1(paste0("#'   As grid_start has a random component, when",
              " \\code{grid_start != 'none'} in your model fitting, you can\n")) 
  cat1(paste0("#'    get slightly different results each time you fit a model\n")) 
  cat1(paste0("#'   \n")) 
  
  cat1(paste0("#'    The parameter confidence intervals returned in sigConf are just",
              " simple confidence intervals, calculated as 2 * standard error.\n"))
  
  if (funName == "sar_power"){
    cat1(paste0("#'   For the power model (and only this model) the returned object ",
                "(sigConf) and model summary also includes the parameter estimates\n",
                "#' generated from fitting the model using \\code{\\link{nls}} and",
                " using as starting parameter estimates the parameter values\n",
                "#' from our model fitting. This also returns the confidence",
                " intervals generated with \\code{\\link{confint}} (which\n",
                "#' calls MASS:::confint.nls), which should be more accurate",
                " than the default \\code{sars} CIs.\n"))
  }
  
  cat1(paste0("#' @importFrom ", "stats lm quantile\n")) 
  
  cat1(paste0("#' @return ", "A list of class 'sars' with the following", 
              " components: \n"))
  cat1(paste0("#'   \\itemize{\n"))
  cat1(paste0("#'     \\item \\strong{par}  The model parameters\n"))
  cat1(paste0("#'     \\item \\strong{value}  Residual sum of squares\n"))
  cat1(paste0("#'     \\item \\strong{counts}   The number of iterations for the", 
              " convergence of the fitting algorithm\n"))
  cat1(paste0("#'     \\item \\strong{convergence}  Numeric code returned from optim", 
              " indicating model convergence (0 = converged)\n"))
  cat1(paste0("#'     \\item \\strong{message}  Any message from the model fit", 
              " algorithm\n"))
  cat1(paste0("#'     \\item \\strong{hessian}  A symmetric matrix giving an", 
              " estimate of the Hessian at the solution found\n"))
  cat1(paste0("#'     \\item \\strong{verge}  Logical code indicating that optim model", 
              " convergence value is zero\n"))
  cat1(paste0("#'     \\item \\strong{startValues}  The start values for the model", 
              " parameters used in the optimisation\n"))
  cat1(paste0("#'     \\item \\strong{data}  Observed data\n"))
  cat1(paste0("#'     \\item \\strong{model}  A list of model information (e.g. the", 
              " model name and formula)\n"))
  cat1(paste0("#'     \\item \\strong{calculated}   The fitted values of the model",
              "\n"))
  cat1(paste0("#'     \\item \\strong{residuals}  The model residuals\n"))
  cat1(paste0("#'     \\item \\strong{AIC}  The AIC value of the model\n"))
  cat1(paste0("#'     \\item \\strong{AICc}  The AICc value of the model\n"))
  cat1(paste0("#'     \\item \\strong{BIC}  The BIC value of the model\n"))
  cat1(paste0("#'     \\item \\strong{R2}  The R2 value of the model\n"))
  cat1(paste0("#'     \\item \\strong{R2a}  The adjusted R2 value of the model\n"))
  cat1(paste0("#'     \\item \\strong{sigConf}  The model coefficients table\n"))
  cat1(paste0("#'     \\item \\strong{normaTest}  The results of the residuals", 
              " normality test", "\n"))
  cat1(paste0("#'     \\item \\strong{homoTest}  The results of the residuals", 
              " homogeneity test\n"))
  cat1(paste0("#'     \\item \\strong{observed_shape}  The observed shape of the", 
              " model fit\n"))
  cat1(paste0("#'     \\item \\strong{asymptote}  A logical value indicating whether", 
              " the observed fit is asymptotic\n"))
  cat1(paste0("#'     \\item \\strong{neg_check}  A logical value indicating whether",
              " negative fitted values have been returned}\n\n"))
  cat1(paste0("#'   The \\code{\\link{summary.sars}} function returns a more", 
              " useful summary of\n"))
  cat1(paste0("#'   the model fit results, and the \\code{\\link{plot.sars}}", 
              " plots the model fit.\n"))
  cat1(paste0("#' @references Triantis, K.A., Guilhaumon, F. & Whittaker,", 
              " R.J. (2012) The island species-area", "\n")) 
  cat1(paste0("#'   relationship: biology and statistics. Journal of", 
              " Biogeography, 39, 215-231.\n"))
  
  #grid_start makes example too long for betap, so just for this model turn
  #grid_start off. For mmf, we need to suppress the warning that it is
  #deprecated
  if (funName == "sar_betap"){
    cat1(paste0("#' @examples", "\n", 
                "#' #Grid_start turned off for speed (not recommended)", "\n",
                "#' data(galap)", "\n",
                "#' fit <- ", funName,
                "(galap, grid_start = 'none')", "\n", 
                "#' summary(fit)","\n", "#' plot(fit)\n"))
  } else if (funName== "sar_mmf"){
    cat1(paste0("#' @examples", "\n", "#' data(galap)", "\n",
                "#' fit <- ", "suppressWarnings(", funName,
                "(galap))", "\n", "#' summary(fit)","\n", "#' plot(fit)\n"))
  } else{
    cat1(paste0("#' @examples", "\n", "#' data(galap)", "\n",
                "#' fit <- ", funName,
                "(galap)", "\n", "#' summary(fit)","\n", "#' plot(fit)\n"))
  }

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
  # if(fExists & !overwrite) stop("a model function '", funName, "' 
  # already exists. If you know what you do, please set overwrite=TRUE")
  #if(!(fExists) | overwrite)
  
  file.create(file.path("R",fileName), overwrite = overwrite, 
              showWarnings = FALSE)
  #cat the roxygen2 comments
  cat_roxygen(model, funName, file.path("R",fileName))
  cat1("\n")
  cat1("\n")
  
  #function definition
  cat1(paste0(funName,' <- function(data, start = NULL,', 
                        ' \ngrid_start = "partial", grid_n = NULL,', 
                        ' \nnormaTest =  "none", homoTest =',
                        ' \n"none", homoCor = "spearman", verb = TRUE){',"\n"))
  
  #add deprecation warning to mmf function
  if (funName == "sar_mmf"){
    cat1(".Deprecated()\n")
  }
  
  #checks
  cat1("if (!(is.matrix(data) | is.data.frame(data)))", 
       " \nstop('data must be a matrix or dataframe')\n")
  cat1("data <- as.data.frame(data)\n")
  cat1("if (anyNA(data)) stop('NAs present in data')\n")
  cat1("normaTest <- match.arg(normaTest, c('none', 'shapiro', 'kolmo',\n")
  cat1("'lillie'))","\n")
  cat1("homoTest <- match.arg(homoTest, c('none', 'cor.area',\n")
  cat1("'cor.fitted'))","\n")
  cat1("if (homoTest != 'none'){\n")
  cat1("homoCor <- match.arg(homoCor, c('spearman', 'pearson',\n")
  cat1("'kendall'))","\n")
  cat1("}\n")
  cat1("if (!(grid_start %in% c('none', 'partial', 'exhaustive'))){\n")
  cat1("stop('grid_start should be one of none, partial or exhaustive')\n")
  cat1("}\n")
  cat1("if (grid_start == 'exhaustive'){\n")
  cat1("  if (!is.numeric(grid_n))\n")
  cat1("  stop('grid_n should be numeric if grid_start == exhaustive')\n")
  cat1("  }\n")
  cat1("if (!is.logical(verb)){\n")
  cat1("stop('verb should be logical')\n")
  cat1("}\n")
  #data ordering and column naming (assuming Area then Species Richness)
  cat1("data <- data[order(data[,1]),]\n")
  cat1("colnames(data) <- c('A','S')\n")
  
  cat1("#check for all equal richness values (particuarly zeros)\n")
  cat1("xr <- range(data$S)/mean(data$S)\n")
  cat1("if (isTRUE(all.equal(xr[1], xr[2]))) {\n")
  cat1("  if (data$S[1] == 0){\n")
  cat1("   warning('All richness values are zero: parameter estimates of',\n")
  cat1("           ' non-linear models should be interpreted with caution')\n")
  cat1("     } else{\n")
  cat1("       warning('All richness values identical')\n")
  cat1("     }}\n")

  #model definition (Appending it)
  #dump("model", file = filePath, append = TRUE) #THE DUMP IS NOT 
  #PASTING THE MODEL LIST
  file.append(file.path("R",fileName), 
              file.path(getwd(),"inst","non_lin_models",f))
                
  
  cat1("\n")
  
  cat1("model <- compmod(model)\n")
  
  cat1("fit <- get_fit(model = model, data = data, start = start,", 
      " \ngrid_start = grid_start, grid_n = grid_n, algo = 'Nelder-Mead', 
       normaTest =  normaTest, homoTest = homoTest, 
       homoCor = homoCor, verb = verb)\n")
  cat1("if(is.na(fit$value)){\n")
  cat1("  return(list(value = NA))\n")
  cat1("}else{","\n")
  cat1("  obs <- obs_shape(fit, verb = verb)\n")
  cat1("  fit$observed_shape <- obs$fitShape\n")
  cat1("  fit$asymptote <- obs$asymp\n")
  cat1("  fit$neg_check <- any(fit$calculated < 0)\n")
  cat1("  class(fit) <- 'sars'\n")
  cat1("  attr(fit, 'type') <- 'fit'\n")
  cat1("  return(fit)\n")
  cat1("}\n")
  
  #function end
  cat1(paste0("}#end of ",funName,"\n"))
  
}#eo model_factory


########################
#using the model factory
########################
#setwd("E:/Working directory/sars")
#setwd("C:/Users/Tom/Desktop/sars")
modFiles <- list.files(file.path("C:/Users/Tom/Desktop/sars","inst","non_lin_models"))
lapply(modFiles, model_factory, overwrite = TRUE)
