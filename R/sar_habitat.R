########################################################################
##internal function for fitting habitat nls models using ######
## grid search
########################################################################

#' function for using grid search of nls parameter space
#' @importFrom stats AIC
#' @importFrom minpack.lm nlsLM nls.lm.control 
#' @noRd
habitat_optim <- function(mod_nam, data){
  
  start.list <- list(
    seq(0,50,5),
    c(0.1,0.25,0.75,1),
    c(-1, 0.00001, 0.0001, 0.001, 
      0.01, 0.1, 1, 10, 50))
  
  names(start.list) <- c("c1", "z", "d")
  
  grid.start <- expand.grid(start.list)
  
  mod_nam2 <- switch(mod_nam,
                     "Kallimanis" = formula(S ~ c1 * A^(z + d * H)),
                     "jigsaw" = formula(S ~ (c1 * H^d) * ((A / H)^z)))
  
  fit.list <- suppressWarnings(apply(grid.start, 1, function(x){
    tryCatch(minpack.lm::nlsLM(mod_nam2,
                 start = x,
                 control = minpack.lm::nls.lm.control(maxiter = 1000, 
                                          maxfev = 100000),
                 data = data),
             error = function(e) NA)
  }))
  
  len.fit.list <- sapply(fit.list, length)
  if (any(len.fit.list > 1)){
    good.fit.list <- which(len.fit.list > 1)
    new.fit.list <- fit.list[good.fit.list]
    AIC.fit.list <- vapply(new.fit.list, AIC, 
                           FUN.VALUE = numeric(1))
    #if multiple min, it just picks the first
    best.fit <- new.fit.list[[which.min(AIC.fit.list)]]
  } else {
    best.fit <- NA 
  }
  return(best.fit)
}


#' Fit habitat SAR models
#'
#' @description Fit three SAR regression models that include habitat diversity:
#' the choros model, the Kallimanis model, and the jigsaw model.
#' @usage sar_habitat(data, modType = "power_log", con = NULL, logT = log)
#' @param data A dataset in the form of a dataframe with at least three columns:
#'   the first with island/site areas (A), the second with island / site habitat
#'   diversity (H), and the third with the species richness of each island/site
#'   (S).
#' @param modType What underlying SAR model form should be used. Should be one
#'   of "power" (non-linear power), "logarithmic" (logarithmic SAR), or
#'   "power_log" (log-log power; default).
#' @param con The constant to add to the species richness values in cases where
#'   at least one of the islands has zero species.
#' @param logT The log-transformation to apply to the area and richness values.
#'   Can be any of \code{log}(default), \code{log2} or \code{log10}.
#' @details These functions are described in more detail in the accompanying paper
#'   (Furness et al., 2023). The code to fit the models was also taken from this
#'   paper.
#'   
#'   Three habitat SAR models are available:
#'   \itemize{ \item{choros model:} { Proposes that species richness is better
#'   predicted by the product of habitat heterogeneity and area (S = c.(A.H)^z) }
#'   \item{Kallimanis model:} {
#'   Proposes that increasing habitat heterogeneity increases species richness
#'   by increasing the slope (on a log-log plot) of the Arrhenius power model
#'   (S = c1.A^(z + d.H)) }
#'   \item{jigsaw model:} {
#'   Models species richness in an area as the sum of the species richness
#'   values of several smaller component subareas, which can be visualised as
#'   pieces of a jigsaw puzzle, i.e., it partitions the species–area and
#'   species–heterogeneity scaling relationships (S = (c1.H^d).((A / H)^z)) }}
#'   
#'   In addition to these three models, a simple 'non-habitat' SAR model is also
#'   fit, which varies depending on \code{modType}: the non-linear power, the
#'   logarithmic or the log-log power model.
#'   
#'   The untransformed (\code{modType = "power"}) and logarithmic (\code{modType
#'   = "logarithmic"}) models are fitted using non-linear regression and the
#'   \code{\link{nlsLM}} function. For the jigsaw and Kallimanis
#'   models in untransformed space, a grid search process is used to test
#'   multiple starting parameter values for the \code{\link{nlsLM}} function - see
#'   details in the documentation for \code{\link{sar_average}} - if multiple
#'   model fits are returned, the fit with the lowest \code{AIC} is returned.
#'   Providing starting parameter estimates for multiple datasets is tricky, and
#'   thus you may find the jigsaw and Kallimanis models cannot be fitted in
#'   untransformed space or with the logarithmic models. If this is the case,
#'   please let the package maintainer know and we can edit the starting
#'   parameter values. The log-log models (\code{modType = "power_log"}) are all
#'   fitted using linear regression ( \code{\link{lm}} function).
#'   
#'   \code{sar_habitat()} uses the \code{\link{nlsLM}} from the
#'   \code{minpack.lm} package rather than \code{\link{nls}} as elsewhere in the
#'   package as we found that this resulted in better searches of the parameter
#'   space for the habitat models (and less convergence errors), particularly
#'   for the logarithmic models. \code{\link{nlsLM}} is a modified version of
#'   \code{\link{nls}} that uses the Levenberg-Marquardt fitting algorithm, but
#'   returns a standard \code{\link{nls}} object and thus all the normal
#'   subsequent \code{\link{nls}} functions can be used. Note also that
#'   occasionally a warning is returned of NaNs being present, normally relating
#'   to the jigsaw model (logarithmic version). We believe this mostly relates
#'   to models fitted during the optimisation process rather than the final
#'   returned model. Nonetheless, users are still recommended to check the
#'   convergence information of the returned model fits.
#'   
#' @return A list of class "habitat" and "sars" with up to four elements, each
#'   holding one of the individual model fit objects (either \code{\link{nls}}
#'   or \code{\link{lm}} class objects). \code{\link{summary.sars}} provides a
#'   more user-friendly ouput (including a model summary table ranked by AICc
#'   and presenting the model coefficients, and R2 and information criteria
#'   values etc.) and \code{\link{plot.habitat}} provides a simple bar of
#'   information criteria weights. For the models fitted using non-linear
#'   regression, the R2 and adjusted R2 are 'pseudo R2' values and are
#'   calculated using the same approach as in the rest of the package (e.g.,
#'   \code{\link{sar_power}}.
#'   
#'   Note that if any of the models cannot be fitted - this is particularly the
#'   case when fitting the untransformed or logarithmic models which use
#'   non-linear regression (see above) - they are removed from the returned
#'   object.
#' @note The jigsaw model is equivalent to the trivariate power-law model of 
#' Tjørve (2009), see Furness et al. (2023).
#' 
#' The jigsaw model (power-law form) cannot have a poorer fit than the choros or
#' power model based on RSS and thus R2. Comparing models using information
#' criteria is thus advised.
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom stats lm
#' @references Furness, E.N., Saupe, E.E., Garwood, R.J., Mannion, P.D. &
#'   Sutton, M.D. (2023) The jigsaw model: a biogeographic model that partitions
#'   habitat heterogeneity from area. Frontiers of Biogeography, 15, e58477.
#'   
#'   Kallimanis, A.S., Mazaris, A.D., Tzanopoulos, J., Halley, J.M., Pantis,
#'   J.D., & Sgardelis, S.P. (2008) How does habitat diversity affect the
#'   species–area relationship? Global Ecology and Biogeography, 17, 532-538
#'   
#'   Tjørve, E. (2009) Shapes and functions of species– area curves (II): a
#'   review of new models and parameterizations. Journal of Biogeography, 36,
#'   1435-1445.
#'   
#'   Triantis, K.A., Mylonas, M., Lika, K. & Vardinoyannis, K. (2003) A model
#'   for the species-area-habitat relationship. Journal of Biogeography, 30,
#'   19–27.
#' @author Euan N. Furness and Thomas J. Matthews
#' @examples
#' data(habitat)
#' #Fit the models in log-log space
#' s <- sar_habitat(data = habitat, modType = "power_log", 
#' con = NULL, logT = log)
#' #Look at the model comparison summary
#' s2 <- summary(s)
#' s2
#' #Make a simple plot of AICc weights
#' plot(s, IC = "AICc", col = "darkred")
#' 
#' #Fit the logarithmic version of the models
# s3 <- sar_habitat(data = habitat, modType = "logarithmic",
# con = NULL, logT = log)
# summary(s3)
# plot(s, IC = "BIC", col = "darkblue")
#' @export

sar_habitat <- function(data, modType = "power_log", 
                          con = NULL, logT = log){
  
  if (!(is.matrix(data) | is.data.frame(data)))
    stop('data must be a matrix or dataframe')
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop('NAs present in data')
  if (!any(c("power", "logarithmic", "power_log") %in% 
           modType)){
    stop("modType should be one of 'power', 'logarithmic' or 'power_log'")
  }
  if (!is.primitive(logT)) stop("logT should be a (primitive) function,
                                specifically: log, log2 or log10")

  data <- data[,1:3]
  data <- data[order(data[,1]),]
  colnames(data) <- c('A','H', 'S')
  
  ##Creating additional variables
  data$choros <- data$A * data$H #both untransformed
  data$choros_log <- logT(data$choros)
  data$HlogA <- data$H * logT(data$A) #untransformed H
  
  ##log conversion (if needed)
  if (modType == "power_log"){
    data$A <- logT(data$A)
    data$H <- logT(data$H)
    if (any(data$S == 0)){
      if (any(length(con) > 1 | !(is.numeric(con)))){
        stop("The dataset has richness values of zero, ",
        "con should be a numeric vector of length 1")
      }
      message("\nThe dataset has zero richness values, ", con,
              " has been added to all richness values.\n\n")
      data$S <- logT(data$S + con)
    } else{
      data$S <- logT(data$S)
    }
  }
  
  if (modType == "power"){
  #Fit nls for each of the four tested models in untransformed space
    choros <- tryCatch(minpack.lm::nlsLM(S ~ c1 * choros^z,
                           start = list("c1" = 5, 
                                        "z" = 0.25),
                           control = minpack.lm::nls.lm.control(maxiter = 1000, 
                                                    maxfev = 100000),
                           data = data),
                       error = function(e) NA)
    
    jigsaw <- habitat_optim("jigsaw", data)

    Kallimanis <- habitat_optim("Kallimanis", data)

    classical <- tryCatch(minpack.lm::nlsLM(S ~ c1 * A^z,
                              start = list("c1" = 5,
                                           "z" = 0.25),
                              control = minpack.lm::nls.lm.control(maxiter = 1000, 
                                                     maxfev = 100000),
                              data = data),
                          error = function(e) NA)

  } else if (modType == "logarithmic"){
    # Fit nls for each of the four tested models in semi-log space
    choros <- tryCatch(minpack.lm::nlsLM(S ~ c1 + z*choros_log,
                              start = list("c1" = 5, 
                                           "z" = 0.25),
                              control = minpack.lm::nls.lm.control(maxiter = 1000, 
                                                       maxfev = 100000),
                              data = data),
                          error = function(e) NA)
    
    jigsaw <- tryCatch(minpack.lm::nlsLM(S ~ (H^d) * logT(c1 * (A / H)^z),
                           start = list("c1" = 5,
                                        "z" = 1,
                                        "d" = 0.6),
                           control = minpack.lm::nls.lm.control(maxiter = 1000, 
                                                    maxfev = 100000),
                           data = data),
                       error = function(e) NA)

    Kallimanis <- tryCatch(minpack.lm::nlsLM(S ~ c1 + (z + d * H) * logT(A),
                           start = list("c1" = 5,
                                        "z" = 1,
                                        "d" = 0.6),
                           control = minpack.lm::nls.lm.control(maxiter = 1000, 
                                                    maxfev = 100000),
                           data = data),
                       error = function(e) NA)

    classical <- tryCatch(minpack.lm::nlsLM(S ~ c1 + z * logT(A),
                             start = list("c1" = 5,
                                          "z" = 0.25),
                             control = minpack.lm::nls.lm.control(maxiter = 1000, 
                                                      maxfev = 100000),
                             data = data),
                         error = function(e) NA)

  } else if (modType == "power_log"){
    # Fit four models in log-log space
    choros <- lm(S ~ choros_log, data = data)
    jigsaw <- lm(S ~ A + H, data = data)
    #I've checked and this log-log form of Kallimanis
    #matches if you take the log of both sides of the
    #untransformed model
    Kallimanis <- lm(S ~ A + HlogA, data = data)
    classical <- lm(S ~ A, data = data)
  }

  res <- list("choros" = choros, "jigsaw" = jigsaw,
              "Kallimanis" = Kallimanis, 
              "power" = classical)
  if (modType == "logarithmic"){
    names(res)[which(names(res) == "power")] <-
      "logarithmic"
  }
  
  attr(res, "failedMods") <- "none"
  
  #check for any NAs in the nls models
  res_len <- vapply(res, length, FUN.VALUE = numeric(1))
  if (all(res_len == 1)){
    stop("No model could be fitted given the starting parameters")
  }
  if (any(res_len == 1)){
    wNA <- which(res_len == 1)
    message("The following models could not be fitted",
            " given the starting parameters and have been excluded: ",
            paste(names(res_len)[wNA], collapse = ", "))
    res <- res[-wNA]
    attr(res, "failedMods") <- names(res_len)[wNA]
  }

  class(res) <- c("habitat", "sars", "list")
  attr(res, "type") <- "habitat"
  attr(res, "modType") <- modType
  return(res)
}
