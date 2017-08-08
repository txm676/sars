#'  Plot Model Fits for an mmSAR2 Object
#'
#' @description S3 method for class 'mmSAR2'. plot.mmSAR2 creates plots for
#'   objects of class mmSAR2, using the ggplot2 R package. The exact plot(s)
#'   constructed depends on the 'Type' attribute (e.g. 'lin_pow') of the mmSAR2
#'   object.
#'
#'   For an mmSAR2 object of Type 'lin_pow' (i.e. linear power model fit), the
#'   plot.mmSAR2 function returns a plot of the model fit (black line) and the
#'   observed richness values (coloured circles).
#' @param object An object of class 'mmSAR2'.
#' @param lsi Argument for line width (default = 3).
#' @param ps Argument for point size (default = 1).
#' @param pc Argument for point colour (default = darkgreen).
#' @param ...	further arguments passed to or from other methods.
#' @examples
#' data(galap)
#' fit <- lin_pow(galap, a = 1, s = 2, con = 1)
#' plot(fit)
#' @export




plot.mmsar2 <- function(object, lsi = 3, ps = 1, pc = "darkgreen", title = NULL,...){

    if (is.null(title)) title <- object$model$name




  int_plot <- function(x, title, s1 = 2, s2 = 1, s3 = 14, s4 = 13, s5 = 15,
                       c1 = "darkred", c2 = "black", xl = "Area", yl = "Species richness",
                       p1 = 0){
    df <- x$data
    xx <- df$A
    yy <- df$S
    ff <- x$calculated

    g <- ggplot2::ggplot(data = df) +
      ggplot2::geom_point(ggplot2::aes(x = xx, y = yy), size = s1, col = c1) +
      ggplot2::geom_line(ggplot2::aes(x = xx, y = ff), size = s2, col = c2) +
      ggplot2::xlab(xl) + ggplot2::ylab(yl) +
      ggplot2::ggtitle(title) + ggplot2::theme_bw() +
      ggplot2::theme(axis.title = ggplot2::element_text(size = s3), axis.text =
                       ggplot2::element_text(size = s4)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = p1, size = s5))

    return(g)
  }



  }









  if (length(object == 1))





  if (attributes(object)$Type == "lin_pow"){
    logDat <- log(object$Area)
    con <- attributes(object)$Constant
    if (any(object$Richness == 0)){
      lor = log(object$Richness + con)
    } else {
      lor = log(object$Richness)
    }
    gdf <- data.frame(Area = logDat, fit_ric = object$Fitted, obs_ric = lor)
    plot(gdf$Area, gdf$fit_ric, col = "white", xlab = "Area (log transformed)",
         ylab = "Richness (log transformed)", cex.lab = 1.3, main = attributes(object)$Dataset, ...)
    lines(gdf$Area, gdf$fit_ric, lwd = lsi)
    points(gdf$Area, gdf$obs_ric, pch = 16, col = pc, cex = ps)
  }
}















plot.mmSAR2 <- function(object, lsi = 3, ps = 1, pc = "darkgreen", ...){

  if (attributes(object)$Type == "lin_pow"){
    logDat <- log(object$Area)
    con <- attributes(object)$Constant
    if (any(object$Richness == 0)){
      lor = log(object$Richness + con)
    } else {
      lor = log(object$Richness)
    }
    gdf <- data.frame(Area = logDat, fit_ric = object$Fitted, obs_ric = lor)
    plot(gdf$Area, gdf$fit_ric, col = "white", xlab = "Area (log transformed)",
         ylab = "Richness (log transformed)", cex.lab = 1.3, main = attributes(object)$Dataset, ...)
    lines(gdf$Area, gdf$fit_ric, lwd = lsi)
    points(gdf$Area, gdf$obs_ric, pch = 16, col = pc, cex = ps)
  }
}





