#'  Plot Model Fits for an mmsar2 Object
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
#' @export


plot.mmsar2 <- function(x, title = NULL, s1 = 2, s2 = 1, s3 = 14, s4 = 13, s5 = 15,
                        c1 = "darkred", c2 = "black", xl = "Area", yl = "Species richness",
                        p1 = 0)
{

    if (is.null(title)) title <- object$model$name

    g1 <- int_plot(x, title, s1, s2, s3, s4, s5,
                   c1, c2, xl, yl, p1)

    return(g1)
}




