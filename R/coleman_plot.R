#' Plot Model Fits for a 'coleman' Object
#'
#' @description S3 method for class 'coleman'. \code{plot.coleman} creates a
#'   plot for objects of class coleman, using the R base plotting framework.
#' @param x An object of class 'coleman'.
#' @param xlab Title for the x-axis.
#' @param ylab Title for the y-axis.
#' @param pch Plotting character (for points).
#' @param cex A numerical vector giving the amount by which plotting symbols
#'   (points) should be scaled relative to the default.
#' @param pcol Colour of the points.
#' @param cex.lab The amount by which the the axis titles should be scaled
#'   relative to the default.
#' @param cex.axis The amount by which the the axis labels should be scaled
#'   relative to the default.
#' @param lwd Line width.
#' @param lcol1 Line colour of the fitted model curve.
#' @param lcol2 Line colour of the model standard deviation curves.
#' @param ModTitle Plot title (default is null, which equates to no main
#'   title).
#' @param TiAdj Which way the plot title (if included) is justified.
#' @param TiLine Places the plot title (if included) this many lines outwards
#'   from the plot edge.
#' @param cex.main The amount by which the the plot title (if included)
#'   should be scaled relative to the default.
#' @param \dots Further graphical parameters (see
#'   \code{\link[graphics]{par}},
#'   \code{\link[graphics]{plot}},\code{\link[graphics]{title}},
#'   \code{\link[graphics]{lines}}) may be supplied as arguments.
#' @details The resultant plot contains the observed richness values with the
#'   model fit and confidence intervals. Following Wang et al. (2010), the
#'   model is rejected if more than a third of the observed data points fall
#'   beyond one standard deviation from the expected curve.
#' @usage coleman(data, area)
#' @importFrom dplyr arrange_
#' @import graphics
#' @examples
#' data(cole_sim)
#' fit <- coleman(cole_sim[[1]], cole_sim[[2]])
#' plot(fit, ModTitle = "Hetfield")
#' @export


plot.coleman <- function(x, xlab = "Relative area (log transformed)", 
                         ylab = "Species richness", pch = 16, cex = 1.2, 
           pcol = 'black', cex.lab = 1.3, cex.axis = 1,
           lwd = 2, lcol1 = 'black', lcol2 = "darkgrey", ModTitle = NULL, 
           TiAdj = 0, TiLine = 0.5, cex.main = 1.5, ...)
{
    df <- data.frame(x$Predicted_values, x$Standard_deviation,
                   x$Relative_areas, x$Species_richness)
    colnames(df) <- c("pv", "sd", "ra", "os")
    
    df <- arrange_(df, ~ra)#using standard evaluation
    
    plot(x = log(df$ra), y = df$os, ylim = c(min((df$pv - df$sd)),
                                             max((df$pv + df$sd))),
         xlab = xlab, ylab = ylab, pch = pch, cex = cex, col = pcol,
         cex.lab = cex.lab, cex.axis = cex.axis, ...)
    lines(x = log(df$ra), y = (df$pv + df$sd), col = lcol2, lwd = lwd, ...)
    lines(x = log(df$ra), y = (df$pv - df$sd), col = lcol2, lwd = lwd, ...)
    lines(x = log(df$ra), y = df$pv, col = lcol1, lwd = lwd, ...)
    if (!is.null(ModTitle)) title(main = ModTitle, adj = TiAdj,
                                  line = TiLine, 
                                  cex.main = cex.main, ...)
    }

