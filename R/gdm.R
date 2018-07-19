###########Fit the GDM#############################

#data <- data.frame("A" = c(10,40,80,160,160), "S" = c(1,3,5,8,10), Ti = c(1,2,3,4,5))

#' Fit the General Dynamic Model of Island Biogeography
#'
#' @description Fit the general dynamic model (GDM) of island biogeography using
#'   a variety of SAR models. Functions are provided to compare the GDM fitted
#'   using different SAR models, and also, for a given SAR model, to compare the
#'   GDM with alternative nested candidate models (e.g. S ~ A + T).
#' @usage gdm(data, model = "lin_pow", mod_sel = FALSE, A = 1, S = 2, Ti = 3)
#' @param data A dataframe or matrix with at least three columns, where one
#'   column should include island area values, one island richness values and
#'   one island age values.
#' @param model Name of the SAR model to be used to fit the GDM. Can be any of
#'   'expo', 'linear', 'power', 'all', or 'lin_pow'.
#' @param mod_sel Logical argument specifying whether, for a given SAR model, a
#'   model comparison of the GDM with other nested candidate models should be
#'   undertaken.
#' @param AST The column locations in \code{data} for the area, richness and
#'   time values (in that order).
#' @details The GDM models island species richness as a function of island area
#'   and island age, and takes the general form: S ~ A + T + T^2, where S =
#'   richness, A =area, and T = island age. The T^2 term is included as the GDM
#'   predicts a hump-shaped relationship between island richness and island age.
#'   However, a variety of different SAR models have been used to fit th GDM and
#'   four options are available here: the exponential, linear and power SAR
#'   models, and the log-log version of the power model. For example, when the
#'   linear SAR model is used, the GDM can be fitted using the expression: S ~ c
#'   + z*Area + k*T + j*T^2, where c,z,k,j are parameters to be estimated.
#'
#'   For the exponential, linear and power SAR models, the GDM is fitted using
#'   non-linear regression and the \code{\link{nls}} function. In the case of
#'   the log-log power model, the GDM is fitted using the \code{\link{lm}}
#'   function.
#'
#'   Residual standard error (RSE) is used instead of R2, and for each model the
#'   AIC and RSE is provided.
#'
#'
#' @return An object of class 'gdm'. If \code{model == "lin_pow"}, the returned
#'   object is a \code{\link{lm}} model fit. If \code{model %in% c("expo",
#'   "linear", "power")} the returned object is a \code{\link{nls}} model fit
#'   object. If \code{model == "all"}, the returned object is a list with three
#'   elements; each element being a \code{nls} fit object.
#'
#'   If \code{mod_sel == TRUE} and \code{model != "all"}, a list with four
#'   elements is returned; each element being a \code{lm} or \code{nls} fit
#'   object. When \code{model == "all"}, a list with three elements is returned;
#'   each element being a list of the four model fits for a particular SAR
#'   model.
#' @note AIC is calculated using the \code{\link{AIC}} function, which is based
#'   on the log-likelihood and not the residual sum of squares (the latter is
#'   used in the main functions of the sars package).
#'   
#'   A plot generic function enabling 3-d plotting of the GDM fit will be
#'   provided in a future version of the package.
#'
#' @import stats
#' @references Whittaker, R. J., Triantis, K. A., & Ladle, R. J. (2008). A
#'   general dynamic theory of oceanic island biogeography. Journal of
#'   Biogeography, 35, 977-994.
#'
#'   Borregaard, M. K. et al. (2017). Oceanic island biogeography through the
#'   lens of the general dynamic model: assessment and prospect. Biological
#'   Reviews, 92, 830-853.
#' @examples
#' #create an example dataset and fit the GDM using the exponential SAR model
#' data(galap)
#' galap$t <- rgamma(16, 5, scale = 2)
#' g <- gdm(galap, model = "expo", mod_sel = FALSE)
#'
#' #Compare the GDM (using the exponential model) with other nested candidate models
#' g2 <- gdm(galap, model = "expo", mod_sel = TRUE)
#'
#' #compare the GDM fitted using the linear, exponential and power SAR models
#' g3 <- gdm(galap, model = "all", mod_sel = FALSE)
#' @export

#no R2 provided for non-linear models as only for linear models, residual standard
#error provided

#null model calculated using lm not nls

#AIC calculated using AIC function (log likelihood) not the 
#RSS like rest of sars: check with FG this is OK

#all model comparison does not include log-log power as can't use AIC

gdm <- function(data, model = "lin_pow", mod_sel = FALSE, AST = c(1, 2, 3)){
  if (anyNA(data)) stop("NAs present in data")
  if (!(is.matrix(data) || is.data.frame(data))) stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)
  if (ncol(data) < 3) stop("Not enough columns/variables to fit GDM")
  if (ncol(data) > 3) {
    warning("More than three columns in dataframe: using the first three")
    data <- data[, 1:3]
  }
  if (!(model %in% c("lin_pow", "power", "expo", "linear", "all"))) {
    stop("provided model name not available")
  }
  if (!is.logical(mod_sel)) stop("mod_sel argument should be TRUE or FALSE")
  
  if (all(AST == c(1, 2, 3))){
    colnames(data) <-c("Area", "SR", "Time")
  } else{
    data <- data[, AST]
    colnames(data) <- c("Area", "SR", "Time")
  }
  
  if (model == "all") allMods <- vector("list", length = 2)

  if (model == "lin_pow"){
    #cat("\n","Fitting the GDM using the linear (log-log) power model", "\n")
    if (any(data$S == 0)) data$S <- data$S + 0.1
    data$Area <- log(data$Area)
    data$SR <- log(data$SR)
    data$Time2 <- data$Time ^ 2
    fit <- lm(SR ~ Area + Time + Time2, data = data)
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 4)
      fitL[[1]] <- fit
      fitL[[2]] <- lm(SR ~ Area + Time, data = data)
      fitL[[3]] <- lm(SR ~ Area, data = data)
      fitL[[4]] <- lm(SR ~ 1, data = data)
      fit <- fitL
    }
    class(fit) <- c("gdm", "lm")
    attr(fit, "Type") <- "lin_pow"
    attr(fit, "mod_sel") <- mod_sel
    
  } 
  if (model == "expo" || model == "all"){
      
    #cat("\n","Fitting the GDM using the exponential model", "\n")

     fit <- nls(formula = SR ~ Int + A * log(Area) + Ti * Time + Ti2 * Time ^ 2, 
                         data = data, start = data.frame(Int = 0, A = 1, Ti = 1, Ti2 = 0))
     if (mod_sel == TRUE){
       fitL <- vector("list", length = 4)
       fitL[[1]] <- fit
       fitL[[2]] <- nls(formula = SR ~ Int + A * log(Area) + Ti * Time, 
                        data = data, start = data.frame(Int = 0, A = 1, Ti = 1))
       fitL[[3]] <- nls(formula = SR ~ Int + A * log(Area), 
                        data = data, start = data.frame(Int = 0, A = 1))
       fitL[[4]] <- lm(SR ~ 1, data = data) #intercept only is mean of Y; so can use lm as no functional form implied
       fit <- fitL
     }
     if (mod_sel){
       class(fit) <- c("gdm")
     } else {
       class(fit) <- c("gdm", "nls")
     }
     attr(fit, "Type") <- "expo"
     attr(fit, "mod_sel") <- mod_sel
     if (model == "all") allMods[[1]] <- fit
     
  } 
  if (model == "linear" || model == "all"){
   # cat("\n","Fitting the GDM using the linear model", "\n")
    fit <- nls(SR ~ Int + A * Area + Ti * Time + Ti2 * Time ^ 2, 
               data = data, start = data.frame(Int = 0, A = 1, Ti = 1, Ti2 = 0))
    
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 4)
      fitL[[1]] <- fit
      fitL[[2]] <- nls(formula = SR ~ Int + A * Area + Ti * Time, 
                       data = data, start = data.frame(Int = 0, A = 1, Ti = 1))
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
  
  if (model == "all") {
    class(allMods) <- "gdm"
    attr(allMods, "Type") <- "allMods"
    attr(allMods, "mod_sel") <- mod_sel
    return(allMods)
  } else {
    return(fit)
  }
}
  

