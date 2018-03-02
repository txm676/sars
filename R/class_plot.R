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
#' @param th an additional ggplot2 theme can be added and instructions can be passed.
#'    As ggplot2 is not *loaded*, th needs to use ggplot2::, e.g. ggplot2::theme_bw
#' @export


plot.sars <- function(x, xlab = NULL, ylab = NULL, pch = 16, cex = 1.2, 
                      pcol = "black", ModTitle = NULL, TiAdj = 0, TiLine = 0.5, lwd = 2,
                      lcol = "black", di = NULL, ...)
{
  
  if (is.null(xlab)){
    if (attributes(x)$type == "fit" || attributes(x)$type == "fit_collection"){
        xlab = "Area"
    } else if (attributes(x)$type == "lin_pow"){
      xlab = "Log(Area)"
    } else {
      stop ("Type attribute not recognised")
    }
  }
  
  if (is.null(ylab)){
    if (attributes(x)$type == "fit" || attributes(x)$type == "fit_collection"){
      ylab = "Species richness"
    } else if (attributes(x)$type == "lin_pow"){
      ylab = "Log(Species richness)"
    }
  }
  
  if (attributes(x)$type == "fit"){
    if (is.null(ModTitle)) ModTitle <- x$model$name
    
    df <- x$data
    xx <- df$A
    yy <- df$S
    ff <- x$calculated
    
    plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, cex = cex, ...)
    title(main = ModTitle, adj = TiAdj, line = TiLine, ...)
    lines(x = xx, y = ff, lwd = lwd, col = lcol, ...)
    
  }
  
  if (attributes(x)$type == "fit_collection"){

    if (is.null(di)) {
      di <- ceiling(length(x) / 2)
      par(mfrow = c(di, di))
    } else {
      par(mfrow = di)
    }
      lapply(x, function(x){
      
      df <- x$data
      xx <- df$A
      yy <- df$S
      ff <- x$calculated
      ModTitle <- x$model$name 
      
      plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, cex = cex, ...)
      title(main = ModTitle, adj = TiAdj, line = TiLine, ...)
      lines(x = xx, y = ff, lwd = lwd, col = lcol, ...)
    })
  }
  
  if (attributes(x)$type == "lin_pow"){
    if (is.null(ModTitle)) ModTitle <- "Log-log power"
 
    df <- x$data
    xx <- df$A
    yy <- df$S
    ff <- x$calculated
    
    plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, cex = cex, ...)
    title(main = ModTitle, adj = TiAdj, line = TiLine, ...)
    lines(x = xx, y = ff, lwd = lwd, col = lcol, ...)
  }
}



