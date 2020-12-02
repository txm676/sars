
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
#'   'loga', 'linear', 'power_area', 'power_area_time', 'all', or 'ATT2'.
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
#'   and five options are available here: four using non-linear regression and
#'   one using linear regression.
#'
#'   Non-linear models
#'
#'   Four SAR models can be used here to fit the GDM: the logarithmic
#'   (\code{model = "loga"}), linear (\code{model = "linear"}) and power
#'   (\code{model = "power_area"}) SAR models. Another variant of the GDM
#'   includes power functions of both area and time (\code{model =
#'   "power_area_time"}). Model fitting follows the procedure in Cardoso et al.
#'   (2015). For example, when the linear SAR model is used, the GDM can be
#'   fitted using the expression: S ~ Int + A*Area + Ti*T + Ti2*T^2, where Int,
#'   A, Ti and Ti2 are free parameters to be estimated. When the power model is
#'   used just for area, the equivalent expression is: S ~ exp(Int + A*log(Area)
#'   + Ti*T + Ti2*T^2). For all four models, the GDM is fitted using non-linear
#'   regression and the \code{\link{nls}} function. It should be noted that the
#'   two power models are fitted using S ~ exp(...) to ensure the same response
#'   variable (i.e. S and not log(S)) is used in all GDM models and thus AIC etc
#'   can be used to compare them.
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
#'   \code{model = "all"}, the GDM is fitted four times (using the power_area,
#'   power_area_time, loga and linear SAR models), and the fits compared using
#'   \code{AIC} and \code{AICc}.
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
#'   values will be provided to each of the four GDM models (power_area,
#'   power_area_time, logarithmic and linear).
#'   
#'   Linear ATT2 Model
#'
#'   As an alternative to fitting the GDM using non-linear regression, the model
#'   can be fitted in various ways using linear regression. This can also be
#'   useful if you are having problems with the non-linear regression algorithms
#'   not converging. If \code{model = "ATT2"} is used, the GDM is fitted using
#'   the semi-log logarithmic SAR model using linear regression (with
#'   untransformed richness and time, and log(area)); this is the original GDM
#'   model fitted by Whittaker et al. (2008) and we have used their chosen name
#'   (ATT2) to represent it. Steinbauer et al. (2013) fitted variants of this
#'   model using linear regression by log-transforming richness and / or time.
#'   While we do not provide functionality for fitting these variants, this is
#'   easily done by simply providing the log-transformed variable values to the
#'   function rather than the untransformed values. Using \code{model = "ATT2"}
#'   is basically a wrapper for the \code{lm} function. If \code{mod_sel ==
#'   TRUE}, the GDM is fitted and compared with three other (nested) candidate
#'   models: log(area) and time (i.e. no time^2 term), just log(area), and an
#'   intercept only model.
#'   
#' @return Different objects are returned depending on whether the non-linear or
#'   linear regression models are fitted.
#'
#'   Non-linear models
#'
#'   An object of class 'gdm'. If \code{model} is one of "loga", "linear",
#'   "power_area" or "power_area_time" the returned object is a
#'   \code{\link{nls}} model fit object. If \code{model == "all"}, the returned
#'   object is a list with four elements; each element being a \code{nls} fit
#'   object. If \code{mod_sel == TRUE} and \code{model != "all"}, a list with
#'   four elements is returned; each element being a \code{lm} or \code{nls} fit
#'   object. When \code{model == "all"}, a list with four elements is returned;
#'   each element being a list of the four model fits for a particular SAR
#'   model.
#'
#'   Linear ATT2 Model
#'
#'   If \code{model = "ATT2"} is used, the returned object is
#'   of class 'gdm' and 'lm' and all of the method functions associated with
#'   standard 'lm' objects (e.g. plot and summary) can be used. If \code{mod_sel
#'   = TRUE} a list with four elements is returned; each element being a
#'   \code{lm} object.
#'
#' @note The intercept (Int) parameter that is returned in the power models fits
#'   (\code{model = "power_area" | "power_area_time"}) is on the log scale.
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
#'   
#'   Steinbauer, M.J., Dolos, K., Field, R., Reineking, B. & Beierkuhnlein, C.
#'   (2013) Re-evaluating the general dynamic theory of oceanic island
#'   biogeography. Frontiers of Biogeography, 5.
#'   
#'   Carey, M., Boland, J., Weigelt, P. & Keppel, G. (2020) Towards an extended
#'   framework for the general dynamic theory of biogeography. Journal of
#'   Biogeography, 47, 2554-2566.
#' @examples
#' #create an example dataset and fit the GDM using the logarithmic SAR model
#' data(galap)
#' galap$t <- c(4, 1, 13, 16, 15, 2, 6, 4, 5, 11, 3, 9, 8, 10, 12, 7)
#' g <- gdm(galap, model = "loga", mod_sel = FALSE)
#'
#' #Compare the GDM (using the logarithmic model) with other nested candidate
#' #models
#' g2 <- gdm(galap, model = "loga", mod_sel = TRUE)
#'
#' #compare the GDM fitted using the linear, logarithmic and both power models
#' g3 <- gdm(galap, model = "all", mod_sel = FALSE)
#'
#' #fit the GDM using the original ATT2 model of Whittaker et al. 2008 using lm,
#' #and compare it with other nested models
#' g4 <- gdm(galap, model = "ATT2", mod_sel = TRUE)
#'
#' #provide different starting parameter values when fitting the non-linear
#' #power model GDM
#' g5 <- gdm(galap, model = "power_area",
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
  if (!(model %in% c("power_area", "power_area_time",
                     "loga", "linear", "all", "ATT2"))) {
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
    start_vals <- data.frame(Int = 1, A = 1, Ti = 1, Ti2 = 0)
  }
  
  if (model == "all") allMods <- vector("list", length = 4)

  if (model == "loga" | model == "all"){
      
    #cat("\n","Fitting the GDM using the logarithmic model", "\n")
     fit <- nls(formula = 
                  SR ~ Int + A * log(Area) + Ti * Time + Ti2 * Time ^ 2, 
                         data = data, start = start_vals)
     if (mod_sel == TRUE){
       fitL <- vector("list", length = 4)
       fitL[[1]] <- fit
       fitL[[2]] <- nls(formula = SR ~ Int + A * log(Area) + Ti * Time, 
                        data = data, start = data.frame(Int = 1,  A = 1, 
                                                        Ti = 1))
       fitL[[3]] <- nls(formula = SR ~ Int + A * log(Area), 
                        data = data, start = data.frame(Int = 1, A = 1))
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
                       start = data.frame(Int = 1, A = 1, Ti = 1))
      fitL[[3]] <- nls(formula = SR ~ Int + A * Area, 
                       data = data, start = data.frame(Int = 1, A = 1))
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
  
  if (model == "power_area" | model == "all"){
   
    fit <- nls(SR ~ exp(Int + A * log(Area) + Ti * Time + Ti2 * Time ^ 2), 
                         data = data, 
                         start = start_vals)
    
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 4)
      fitL[[1]] <- fit
      fitL[[2]] <- nls(formula = SR ~ exp(Int + A * log(Area) + Ti * Time), 
                       data = data, 
                       start = data.frame(Int = 1, A = 1, Ti = 1))
      
      fitL[[3]] <- nls(formula = SR ~ exp(Int + A * log(Area)), 
                       data = data, start = data.frame(Int = 1, A = 1))
      fitL[[4]] <- lm(SR ~ 1, data = data)
      fit <- fitL
    }
    if (mod_sel){
      class(fit) <- c("gdm")
    } else {
      class(fit) <- c("gdm", "nls")
    }
    attr(fit, "Type") <- "power_area"
    attr(fit, "mod_sel") <- mod_sel
    if (model == "all") allMods[[3]] <- fit
  }
  
  if (model == "power_area_time" | model == "all"){
    
    fit <- nls(SR ~ exp(Int + A * log(Area) + Ti * log(Time) + 
                          Ti2 * (log(Time) ^ 2)), 
               data = data, 
               start = start_vals)
    
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 4)
      fitL[[1]] <- fit
      fitL[[2]] <- nls(formula = SR ~ exp(Int + A * log(Area) + Ti * log(Time)), 
                       data = data, 
                       start = data.frame(Int = 1, A = 1, Ti = 1))
      
      fitL[[3]] <- nls(formula = SR ~ exp(Int + A * log(Area)), 
                       data = data, start = data.frame(Int = 1, A = 1))
      
      fitL[[4]] <- lm(SR ~ 1, data = data)
      fit <- fitL
    }
    if (mod_sel){
      class(fit) <- c("gdm")
    } else {
      class(fit) <- c("gdm", "nls")
    }
    attr(fit, "Type") <- "power_area_time"
    attr(fit, "mod_sel") <- mod_sel
    if (model == "all") allMods[[4]] <- fit
  }
  
  if (model == "ATT2"){
    data$Time2 <- data$Time ^ 2
    data$log_Area <- log(data$Area)
   # if (any(data$SR == 0)){
   #   data$SR <- data$SR + 1
    #}
    #data$log_SR <- log(data$SR)
    fit <- lm(SR ~ log_Area + Time + Time2, data = data)
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 4)
      fitL[[1]] <- fit
      fitL[[2]] <- lm(SR ~ log_Area + Time, data = data)
      fitL[[3]] <- lm(SR ~ log_Area, data = data)
      fitL[[4]] <- lm(SR ~ 1, data = data)
      fit <- fitL
    }
    class(fit) <- c("gdm", "lm")
    attr(fit, "Type") <- "ATT2"
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
  




























