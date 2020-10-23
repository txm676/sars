
#' Fit the General Dynamic Model of Island Biogeography
#'
#' @description Fit the general dynamic model (GDM) of island biogeography using
#'   a variety of non-linear and linear SAR models. Functions are provided to
#'   compare the GDM fitted using different SAR models, and also, for a given
#'   SAR model, to compare the GDM with alternative nested candidate models
#'   (e.g. S ~ Area + Time).
#' @usage gdm(data, model = "linear", mod_sel = FALSE, AST = c(1, 2, 3),
#'   start_vals = NULL)
#' @param data A dataframe or matrix with at least three columns, where one
#'   column should include island area values, one island richness values and
#'   one island age values.
#' @param model Name of the SAR model to be used to fit the GDM. Can be any of
#'   'loga', 'linear', 'power', 'all', or 'lin_pow'.
#' @param mod_sel Logical argument specifying whether, for a given SAR model, a
#'   model comparison of the GDM with other nested candidate models should be
#'   undertaken.
#' @param AST The column locations in \code{data} for the area, richness and
#'   time values (in that order).
#' @param start_vals An optional dataframe with starting parameter values for
#'   the non-linear regression models (same format as in \code{\link{nls}}).
#'   Default is set to NULL.
#' @details The GDM models island species richness as a function of island area
#'   and island age, and takes the general form: S ~ A + T + T^2, where S =
#'   richness, A =area, and T = island age. The T^2 term is included as the GDM
#'   predicts a hump-shaped relationship between island richness and island age.
#'   However, a variety of different SAR models have been used to fit the GDM
#'   and four options are available here: three using non-linear regression and
#'   one using linear regression.
#'
#'   Non-linear models
#'
#'   Three SAR models can be used here to fit the GDM: the logarithmic
#'   (\code{model = "loga"}), linear (\code{model = "linear"}) and power
#'   (\code{model = "power"}) SAR models. Model fitting follows the procedure in
#'   Cardoso et al. (2015). For example, when the linear SAR model is used, the
#'   GDM can be fitted using the expression: S ~ Int + A*Area + Ti*T + Ti2*T^2,
#'   where Int, A, Ti and Ti2 are free parameters to be estimated. For all three
#'   models, the GDM is fitted using non-linear regression and the
#'   \code{\link{nls}} function. For ease of fitting, the logarithmic and power
#'   SAR models are included in their logarithmic form, e.g. the power model is
#'   fitted using: S ~ exp(Int + A*log(A)), where Int and A are parameters to be
#'   estimated.
#'
#'   For each model fit, the residual standard error (RSE), R2 and AIC and AICc
#'   values are reported. However, as the model fit object is returned, it is
#'   possible to calculate or extract various other measures of goodness of fit
#'   (see \code{\link{nls}}).
#'
#'   If \code{mod_sel = TRUE}, the GDM (using a particular SAR model) is fitted
#'   and compared with three other (nested) candidate models: area and time
#'   (i.e. no time^2 term), just area, and an intercept only model. The
#'   intercept only model is fitted using \code{lm} rather than \code{nls}. If
#'   \code{model = "all"}, the GDM is fitted three times (using the power, loga
#'   and linear SAR models), and the fits compared using \code{AIC} and
#'   \code{AICc}.
#'
#'   Non-linear regression models are sensitive to the starting parameter values
#'   selected. The defaults used here have been chosen as they provide a
#'   sensible general choice, but they will not work in all circumstances. As
#'   such, alternative starting values can be provided using the
#'   \code{start_vals} argument - this is done in the same way as for
#'   \code{\link{nls}}. The four parameter names are: Int (intercept), A (area),
#'   Ti (Time), Ti2 (Time^2) (see the example below). This only works for the
#'   full GDM non-linear models, and not for the nested models that are fitted
#'   when \code{mod_sel = TRUE} or for the linear models (where they are not
#'   needed). If used with \code{model = "all"}, the same starting parameter
#'   values will be provided to each of the three GDM models (power SAR,
#'   logarithmic SAR and linear SAR).
#'   
#'   Linear Power Model
#'
#'   As an alternative to fitting the GDM with the power SAR using non-linear
#'   regression, the model is often fitted using linear regression (in its
#'   log-log form). This can also be useful if you are having problems with the
#'   non-linear regression algorithms not converging for this model. If
#'   \code{model = "lin_pow"} is used, the GDM is fitted using the log-log
#'   version of the power model (with log(richness) and log(area)) using linear
#'   regression. A constant of 1 is added to all richness values if any zeros
#'   are detected. Using this option is basically a wrapper for the \code{lm}
#'   function. If \code{mod_sel == TRUE}, the GDM is fitted and compared with
#'   three other (nested) candidate models: area and time (i.e. no time^2 term),
#'   just area, and an intercept only model.
#'   
#' @return Different objects are returned depending on whether the non-linear or
#'   linear models are fitted.
#'
#'   Non-linear models
#'
#'   An object of class 'gdm'. If \code{model} is one of "loga", "linear" or
#'   "power" the returned object is a \code{\link{nls}} model fit object. If
#'   \code{model == "all"}, the returned object is a list with three elements;
#'   each element being a \code{nls} fit object. If \code{mod_sel == TRUE} and
#'   \code{model != "all"}, a list with four elements is returned; each element
#'   being a \code{lm} or \code{nls} fit object. When \code{model == "all"}, a
#'   list with three elements is returned; each element being a list of the four
#'   model fits for a particular SAR model.
#'
#'   Linear Power Model
#'
#'   If \code{model = "lin_pow"} is used, the returned object is
#'   of class 'gdm' and 'lm' and all of the method functions associated with
#'   standard 'lm' objects (e.g. plot and summary) can be used. If \code{mod_sel
#'   = TRUE} a list with four elements is returned; each element being a
#'   \code{lm} object.
#'
#' @note The intercept (Int) parameter that is returned in the power model
#'   (\code{model = "power"}) fits is on the log scale.
#'
#'   R2 is calculated using the same function as used in the main sars model
#'   functions.
#'
#' @importFrom stats nls lm
#' @references Whittaker, R. J., Triantis, K. A., & Ladle, R. J. (2008). A
#'   general dynamic theory of oceanic island biogeography. Journal of
#'   Biogeography, 35, 977-994.
#'
#'   Borregaard, M. K. et al. (2017). Oceanic island biogeography through the
#'   lens of the general dynamic model: assessment and prospect. Biological
#'   Reviews, 92, 830-853.
#'
#'   Cardoso, P., Rigal, F., & Carvalho, J. C. (2015). BAT–Biodiversity
#'   Assessment Tools, an R package for the measurement and estimation of alpha
#'   and beta taxon, phylogenetic and functional diversity. Methods in Ecology
#'   and Evolution, 6, 232-236.
#' @examples
#' #create an example dataset and fit the GDM using the logarithmic SAR model
#' data(galap)
#' galap$t <- rgamma(16, 5, scale = 2)
#' g <- gdm(galap, model = "loga", mod_sel = FALSE)
#'
#' #Compare the GDM (using the logarithmic model) with other nested candidate
#' #models
#' g2 <- gdm(galap, model = "loga", mod_sel = TRUE)
#'
#' #compare the GDM fitted using the linear, logarithmic and power SAR models
#' g3 <- gdm(galap, model = "all", mod_sel = FALSE)
#'
#' #fit the GDM using the log-log power model and linear regression, and
#' #compare it with other nested models
#' g4 <- gdm(galap, model = "lin_pow", mod_sel = TRUE)
#'
#' #provide different starting parameter values when fitting the non-linear
#' #power model GDM
#' g5 <- gdm(galap, model = "power",
#' start_vals = data.frame("Int" = 0, "A" = 1, Ti = 1, Ti2 = 0))
#' @export


gdm <- function(data, model = "linear", mod_sel = FALSE, AST = c(1, 2, 3),
                start_vals = NULL){
  if (anyNA(data)) stop("NAs present in data")
  if (!(is.matrix(data) | is.data.frame(data))) 
    stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)
  if (ncol(data) < 3) stop("Not enough columns/variables to fit GDM")
  if (ncol(data) > 3) {
    warning("More than three columns in dataframe: using the first three")
    data <- data[, 1:3]
  }
  if (!(model %in% c("power", "loga", "linear", "all", "lin_pow"))) {
    stop("provided model name not available")
  }
  if (!is.logical(mod_sel)) stop("mod_sel argument should be TRUE or FALSE")
  
  if (all(AST == c(1, 2, 3))){
    colnames(data) <-c("Area", "SR", "Time")
  } else{
    data <- data[, AST]
    colnames(data) <- c("Area", "SR", "Time")
  }
  
  if (!is.null(start_vals)){
    if (!is.data.frame(start_vals)){
      stop("start_vals should be a dataframe with 1 row and four columns")
    }
    if (!all(dim(start_vals) == c(1,4))){
      stop("start_vals should be a dataframe with 1 row and four columns")
    }
    if (!identical(colnames(start_vals), c("Int", "A", "Ti", "Ti2"))){
      stop("colnames of start_vals shoud be c(Int, A, Ti, Ti2)")
    }
  } else{
    start_vals <- data.frame(Int = 0, A = 1, Ti = 1, Ti2 = 0)
  }
  
  if (model == "all") allMods <- vector("list", length = 3)

  if (model == "loga" | model == "all"){
      
    #cat("\n","Fitting the GDM using the logarithmic model", "\n")
     fit <- nls(formula = 
                  SR ~ Int + A * log(Area) + Ti * Time + Ti2 * Time ^ 2, 
                         data = data, start = start_vals)
     if (mod_sel == TRUE){
       fitL <- vector("list", length = 4)
       fitL[[1]] <- fit
       fitL[[2]] <- nls(formula = SR ~ Int + A * log(Area) + Ti * Time, 
                        data = data, start = data.frame(Int = 0,  A = 1, 
                                                        Ti = 1))
       fitL[[3]] <- nls(formula = SR ~ Int + A * log(Area), 
                        data = data, start = data.frame(Int = 0, A = 1))
       fitL[[4]] <- lm(SR ~ 1, data = data) #intercept only is mean of Y; 
       #so can use lm as no functional form implied
       fit <- fitL
     }
     if (mod_sel){
       class(fit) <- c("gdm")
     } else {
       class(fit) <- c("gdm", "nls")
     }
     attr(fit, "Type") <- "loga"
     attr(fit, "mod_sel") <- mod_sel
     if (model == "all") allMods[[1]] <- fit
     
  } 
  if (model == "linear" | model == "all"){
   # cat("\n","Fitting the GDM using the linear model", "\n")
    fit <- nls(SR ~ Int + A * Area + Ti * Time + Ti2 * Time ^ 2, 
               data = data, 
               start = start_vals)
    
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 4)
      fitL[[1]] <- fit
      fitL[[2]] <- nls(formula = SR ~ Int + A * Area + Ti * Time, 
                       data = data, 
                       start = data.frame(Int = 0, A = 1, Ti = 1))
      fitL[[3]] <- nls(formula = SR ~ Int + A * Area, 
                       data = data, start = data.frame(Int = 0, A = 1))
      fitL[[4]] <- lm(SR ~ 1, data = data)
      fit <- fitL
    }
    if (mod_sel){
      class(fit) <- c("gdm")
    } else {
      class(fit) <- c("gdm", "nls")
    }
    attr(fit, "Type") <- "linear"
    attr(fit, "mod_sel") <- mod_sel
    if (model == "all") allMods[[2]] <- fit
  }
  
  if (model == "power" | model == "all"){
   
    fit <- nls(SR ~ exp(Int + A * log(Area) + Ti * Time + Ti2 * Time ^ 2), 
                         data = data, 
                         start = start_vals)
    
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 4)
      fitL[[1]] <- fit
      fitL[[2]] <- nls(formula = SR ~ exp(Int + A * log(Area) + Ti * Time), 
                       data = data, 
                       start = data.frame(Int = 0, A = 1, Ti = 1))
      
      fitL[[3]] <- nls(formula = SR ~ exp(Int + A * log(Area)), 
                       data = data, start = data.frame(Int = 0, A = 1))
      fitL[[4]] <- lm(SR ~ 1, data = data)
      fit <- fitL
    }
    if (mod_sel){
      class(fit) <- c("gdm")
    } else {
      class(fit) <- c("gdm", "nls")
    }
    attr(fit, "Type") <- "power"
    attr(fit, "mod_sel") <- mod_sel
    if (model == "all") allMods[[3]] <- fit
  }
  
  if (model == "lin_pow"){
    data$Time2 <- data$Time ^ 2
    data$log_Area <- log(data$Area)
    if (any(data$SR == 0)){
      data$SR <- data$SR + 1
    }
    data$log_SR <- log(data$SR)
    fit <- lm(log_SR ~ log_Area + Time + Time2, data = data)
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 4)
      fitL[[1]] <- fit
      fitL[[2]] <- lm(log_SR ~ log_Area + Time, data = data)
      fitL[[3]] <- lm(log_SR ~ log_Area, data = data)
      fitL[[4]] <- lm(log_SR ~ 1, data = data)
      fit <- fitL
    }
    class(fit) <- c("gdm", "lm")
    attr(fit, "Type") <- "lin_pow"
    attr(fit, "mod_sel") <- mod_sel
  }

  if (model == "all") {
    class(allMods) <- "gdm"
    attr(allMods, "Type") <- "allMods"
    attr(allMods, "mod_sel") <- mod_sel
    return(allMods)
  } else {
    return(fit)
  }
}
  




























