#' Fit habitat SAR models
#'
#' @description Fit three SAR regression models that include habitat diversity:
#' the choros model, the Kallimanis model, and the jigsaw model.
#' @usage sar_habitat(data, modType = "power_log", con = 1, logT = log)
#' @param data A dataset in the form of a dataframe with at least three columns:
#'   the first with island/site areas, the second with island / site habitat
#'   diversity, and the third with the species richness of each island/site.
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
#'   predicted by the product of habitat heterogeneity and area (S = c(AH)^z) } 
#'   \item{Kallimanis model:} { 
#'   Proposes that increasing habitat heterogeneity increases species richness 
#'   by increasing the slope (on a log-log plot) of the Arrhenius model } 
#'   \item{jigsaw model:} {
#'   Models species richness in an area as the sum of the species richness 
#'   values of several smaller component subareas, which can be visualised as 
#'   pieces of a jigsaw puzzle, i.e., it partitions the species–area and 
#'   species–heterogeneity scaling relationships }}
#'   
#'   In addition to these three models, a simple 'non-habitat' SAR model is also
#'   fit, which varies depending on \code{modType}: the non-linear power, the
#'   logarithmic or the log-log power model.
#'
#' @return A list of class "habitat" and "sars" with four elements, each holding
#'   one of the individual model fit objects (either ** or lm objects).
#'   \code{\link{summary.sars}} provides a more user-friendly ouput (including a
#'   model summary table ranked by AICc and presenting the model coefficients,
#'   and R2 and information criteria values etc.) and \code{\link{plot.habitat}}
#'   provides a simple bar of information criteria weights. For the models
#'   fitted using non-linear regression, the R2 and adjusted R2 are 'pseudo R2'
#'   values and are calculated using the same approach as in the
#'   \code{\link{sar_average}} function.
#' @note The jigsaw model is equivalent to the trivariate power-law model of 
#' Tjørve (2009), see Furness et al. (2023).
#' 
#' The jigsaw model (power-law form) cannot have a poorer fit than the choros or
#' power model based on RSS and thus R2. Comparing models using information
#' criteria is thus advised.
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
#' con = 1, logT = log)
#' #Look at the model comparison summary
#' s2 <- summary(s)
#' s2
#' #Make a simple plot of AICc weights
#' plot(s, IC = "AICc", col = "darkred")
#' @export

sar_habitat <- function(data, modType = "power_log", 
                          con = 1, logT = log){
  
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
  if (any(length(con) > 1 | !(is.numeric(con)))) 
    stop("con should be a numeric vector of length 1")
  
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
      data$S <- logT(data$S + con)
    } else{
      data$S <- logT(data$S)
    }
  }

  ##NEED TO ADD the modType to the summary print text
  if (modType == "logarithmic"){
    # Fit nls for each of the four tested models in semi-log space
    choros <- tryCatch(nls(S ~ c1 + z*choros_log,
                              start = list("c1" = 5, 
                                           "z" = 0.25),
                              data = data),
                          error = function(e) NA)
    
    jigsaw <- tryCatch(nls(S ~ (H^d) * logT(c1 * (A / H)^z),
                           start = list("d" = 0.6,
                                        "c1" = 5,
                                        "z" = 1),
                           data = data),
                       error = function(e) NA)

    Kallimanis <- tryCatch(nls(S ~ c1 + (z + d * H) * logT(A),
                           start = list("c1" = 5,
                                        "z" = 1,
                                        "d" = 0.6),
                           data = data),
                       error = function(e) NA)

    classical <- tryCatch(nls(S ~ c1 + z * logT(A),
                             start = list("c1" = 5,
                                          "z" = 0.25),
                             data = data),
                         error = function(e) NA)
    
    res <- list("choros" = choros, "jigsaw" = jigsaw,
                "Kallimanis" = Kallimanis, 
                "logarithmic" = classical)

  } else if (modType == "power_log"){
    # Fit four models in log-log space
    choros <- lm(S ~ choros_log, data = data)
    jigsaw <- lm(S ~ A + H, data = data)
    Kallimanis <- lm(S ~ A + HlogA, data = data)
    classical <- lm(S ~ A, data = data)
    
    res <- list("choros" = choros, "jigsaw" = jigsaw,
                "Kallimanis" = Kallimanis, 
                "power" = classical)
  }
  
  attr(res, "failedMods") <- "none"
  
  #check for any NAs in the nls models
  res_len <- vapply(res, length, FUN.VALUE = numeric(1))
  if (all(res_len == 1)){
    stop("No nls model could be fitted given the starting parameters")
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
  
  #other habitat models: countryside BG model, matrix and edge-calibrated
  #models, two-habitat SAR, lost-habitat SAR, any in Carey et al.
}
