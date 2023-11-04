#' Fit habitat SAR models
#'
#' @description Fit three SAR regression models including habitat divesity data.
#' @usage sar_habitat(data, logAxes = "both", con = 1, logT = log)
#' @param data A dataset in the form of a dataframe with at least ** columns:
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @param logAxes What log-transformation should be applied to the area and
#'   richness values. Should be one of  "area" (only area is log-transformed;
#'   default) or "both" (both area and richness log-transformed).
#' @param con The constant to add to the species richness values in cases where
#'   one of the islands has zero species.
#' @param logT The log-transformation to apply to the area and richness values.
#'   Can be any of \code{log}(default), \code{log2} or \code{log10}.
#' @details This function is described in more detail in the accompanying paper
#'   (Matthews & Rigal, 2020).
#'
#'   Fitting the continuous and left-horizontal piecewise models (particularly
#'   the two-threshold models) can be time consuming if the range in area is
#'   large and/or the \code{interval} argument is small. For the two-threshold
#'   continuous slope and left-horizontal models, the use of parallel processing
#'   (using the \code{parallel} argument) is recommended. The number of cores
#'   (\code{cores}) must be provided.
#'
#' @return A list of class "threshold" and "sars" with five elements. The first
#'   element contains the different model fits (lm objects). The second element
#'   contains the names of the fitted models, the third  contains the threshold
#'   values, the fourth element the dataset (i.e. a dataframe with area and
#'   richness values), and the fifth contains details of any axes
#'   log-transformations undertaken. \code{\link{summary.sars}} provides a more
#'   user-friendly ouput (including a model summary table) and
#'   \code{\link{plot.threshold}} plots the model fits.
#' @note Due to the
#' @references Furness, E.N., Saupe, E.E., Garwood, R.J., Mannion, P.D. &
#'   Sutton, M.D. (2023) The jigsaw model: a biogeographic model that partitions
#'   habitat heterogeneity from area. Frontiers of Biogeography, 15, e58477.
#'
#'   Triantis, K.A., Mylonas, M., Lika, K. & Vardinoyannis, K. (2003) A model
#'   for the species-area-habitat relationship. Journal of Biogeography, 30,
#'   19–27.
#' @author Euan N. Furness and Thomas J. Matthews
#' @examples
#' data(aegean2)
#' a2 <- aegean2[1:168,]
#' fitT <- sar_threshold(data = a2, mod = c("ContOne", "DiscOne"), 
#' interval = 0.1, non_th_models = TRUE, logAxes = "area", logT = log10)
#' summary(fitT)
#' plot(fitT)
#' #diagnostic plots for the ContOne model
#' par(mfrow=c(2, 2))
#' plot(fitT[[1]][[1]])
#' @export

sar_habitat <- function(data, logAxes = "both", 
                          con = 1, logT = log){
  
  if (!(is.matrix(data) | is.data.frame(data)))
    stop('data must be a matrix or dataframe')
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop('NAs present in data')
  if (!any(c("area", "both") %in% logAxes)){
    stop("logAxes should be one of 'area', or 'both'")
  }
  if (!is.primitive(logT)) stop("logT should be a (primitive) function,
                                specifically: log, log2 or log10")
  if (any(length(con) > 1 | !(is.numeric(con)))) 
    stop("con should be a numeric vector of length 1")
  
  data <- data[order(data[,1]),]
  colnames(data) <- c('A','S')##############
  
  #log conversion (if needed)
  if (logAxes == "area"){
    data$A <- logT(data$A)
  } else if (logAxes == "both"){
    data$A <- logT(data$A)
    if (any(data$S == 0)){
      data$S <- logT(data$S + con)
    } else{
      data$S <- logT(data$S)
    }
  }

  if (logAxes == "area"){
    # Fit nlsLM for each of the four tested models in semi-log space
    choros<-nlsLM(Species~c+z*log(Choros),start=list(c=1,z=0.3),
                  control=nls.lm.control(maxiter=1000,maxfev=100000),data=all)
    jigsaw<-nlsLM(Species~(Het^d)*log(c*(Area/Het)^z),
                  start=list(d=1,c=1,z=0.3),control=nls.lm.control(maxiter=1000, maxfev=100000),data=all)
    Kallimanis<-nlsLM(Species~c+(z+d*Het)*log(Area),
                      start=list(c=1,z=0.3,d=0.1),control=nls.lm.control(maxiter=1000, maxfev=100000),data=all)
    classical<-nlsLM(Species~c+z*log(Area),
                     start=list(c=1,z=0.3),control=nls.lm.control(maxiter=1000, maxfev=100000),data=all)

  } else if (logAxes == "both"){
    # Fit four models in log-log space
    choros <- lm(ln.S. ~ ln.Choros., data = all)
    jigsaw <- lm(ln.S. ~ ln.Area.+ ln.Het., data = all)
    Kallimanis <- lm(ln.S. ~ ln.Area.+ H.ln.A., data = all)
    classical <- lm(ln.S. ~ ln.Area., data = all)
  }
  
  ##Input - could manually provide some parameters (e.g., d of 
  #jigsaw)
  
  ##Return - fitted model objects in a list
  
  ##Summary, then returns the AICc, BIC, AIC and R2 etc in model
  #summary table
  
  #plot: just a plot of AICc weights? 
  
  #other habitat models: countryside BG model, matrix and edge-calibrated
  #models, two-habitat SAR, lost-habitat SAR, any in Carey et al.
  
  ##Questions:
  #should we include an option to manually increase the number of
  #parameters for the the choros model
}
