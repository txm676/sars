###########Fit the GDM#############################

#data <- data.frame("A" = c(10,40,80,160,160), "S" = c(1,3,5,8,10), Ti = c(1,2,3,4,5))

#' Fit the General Dynamic Model of Island Biogeography 
#'
#' @description Fit Coleman's (1981) random placement model to a species-site
#'   abundance matrix: rows are species and columns are sites. Note that the
#'   data must be abundance data and not presence-absence data. According to
#'   this model, the number of species occurring on an island depends on the
#'   relative area of the island and the regional relative species abundances.
#'   The fit of the random placement model can be determined through use of a
#'   diagnostic plot (see \code{\link{plot.coleman}}) of island area (log
#'   transformed) against species richness, alongside the model’s predicted
#'   values (see Wang et al., 2010). Following Wang et al. (2010), the model is
#'   rejected if more than a third of the observed data points fall beyond one
#'   standard deviation from the expected curve.
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
#' @return A list of class "coleman" with four elements. The first element
#'   contains the fitted values of the model. The second element contains the
#'   standard deviations of the fitted values, and the third and fourth contain
#'   the relative island areas and observed richness values, respectively.
#'   \code{\link{plot.coleman}} plots the model.
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
#' #Compare the GDM (using the exponential model) with other **
#' g2 <- gdm(galap, model = "expo", mod_sel = TRUE)
#' 
#' #compare the GDM fitted using the linear, exponential and power SAR models
#' g3 <- gdm(galap, model = "all", mod_sel = FALSE)
#' @export

#no R2 provided for non-linear models as only for linear models, residual standard
#error provided

#AIC calculated using AIC function (log likelihood) not the 
#RSS like rest of sars: check with FG this is OK

#all model comparison does not include log-log power as can't use AIC

#3d plotting will be provided in a future version of the package

#import rgamma

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
       fitL <- vector("list", length = 3)
       fitL[[1]] <- fit
       fitL[[2]] <- nls(formula = SR ~ Int + A * log(Area) + Ti * Time, 
                        data = data, start = data.frame(Int = 0, A = 1, Ti = 1))
       fitL[[3]] <- nls(formula = SR ~ Int + A * log(Area), 
                        data = data, start = data.frame(Int = 0, A = 1))
      # fitL[[4]] <- nls(formula = SR, data = data)
       fit <- fitL
     }
     class(fit) <- c("gdm", "nls")
     attr(fit, "Type") <- "expo"
     attr(fit, "mod_sel") <- mod_sel
     if (model == "all") allMods[[1]] <- fit
     
  } 
  if (model == "linear" || model == "all"){
   # cat("\n","Fitting the GDM using the linear model", "\n")
    fit <- nls(SR ~ Int + A * Area + Ti * Time + Ti2 * Time ^ 2, 
               data = data, start = data.frame(Int = 0, A = 1, Ti = 1, Ti2 = 0))
    
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 3)
      fitL[[1]] <- fit
      fitL[[2]] <- nls(formula = SR ~ Int + A * Area + Ti * Time, 
                       data = data, start = data.frame(Int = 0, A = 1, Ti = 1))
      fitL[[3]] <- nls(formula = SR ~ Int + A * Area, 
                       data = data, start = data.frame(Int = 0, A = 1))
      # fitL[[4]] <- nls(formula = SR, data = data)
      fit <- fitL
    }
    class(fit) <- c("gdm", "nls")
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
  

