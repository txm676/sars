#'  Plot Model Fits for an mmSAR2 Object
#'
#' @description S3 method for class 'mmSAR2'. plot.mmSAR2 creates plots for
#'   objects of class mmSAR2, using the ggplot2 R package. The exact plot(s)
#'   constructed depends on the 'Type' attribute (e.g. 'lin_pow') of the mmSAR2
#'   object.
#' @param object An object of class 'mmSAR2'.
#' @param ...	further arguments passed to or from other methods.**see ggplot2 book for examples
#' @references Chang.
#' @examples
#' data(galap)
#' fit <- lin_pow(galap, a = 1, s = 2, con = 1)
#' plot(fit)
#' @import ggplot2
#' @export

plot.mmSAR2 <- function(object, ...){

  if (attributes(object)$Type == "lin_pow"){
    logDat <- log(object$Area)
    if (any(object$Richness == 0)){
      lor = log(object$Richness + con)
    } else {
      lor = log(object$Richness)
    }
    gdf <- data.frame(Area = logDat, fit_ric = object$Fitted, obs_ric = lor)
    g <- ggplot(data = gdf) + geom_line(aes(Area, fit_ric), size = 1.3) +
      geom_point(aes(Area, obs_ric), size = 1.8, colour = "lightseagreen") +
      ggtitle(attributes(object)$Dataset) +
      theme_bw() + xlab("Area (log transformed)") + ylab("Richness (log transformed)") +
      theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
            plot.title = element_text(size = 15))
  }
  g
}






