#'  Plot Model Fits for a 'sars' Object
#'
#' @description S3 method for class 'sars'. \code{plot.sars} creates plots for
#'   objects of class sars, using the R base plotting framework. The exact
#'   plot(s) constructed depends on the 'Type' attribute of the sars object.For
#'   example, for a sars object of Type 'fit', the \code{plot.sars} function
#'   returns a plot of the model fit (line) and the observed richness values
#'   (points).
#' @param x An object of class 'sars'.
#' @param mfplot Logical argument specifying whether the model fits in a
#'   fit_collection should be plotted on one single plot (\code{mfplot = TRUE}) or
#'   separate plots (\code{mfplot = FALSE}; the default).
#' @param xlab Title for the x-axis (default depends on the Type attribute).
#' @param ylab Title for the y-axis (default depends on the Type attribute).
#' @param pch Plotting character (for points).
#' @param cex A numerical vector giving the amount by which plotting symbols
#'   (points) should be scaled relative to the default.
#' @param pcol Colour of the points.
#' @param ModTitle Plot title (default is null, which reverts to the model
#'   name). For no title, use ModTitle = "". For a sars object of type
#'   fit_collection, a vector of names can be provided (e.g. \code{letters[1:3]}).
#' @param TiAdj Which way the plot title is justified.
#' @param TiLine Places the plot title this many lines outwards from the plot
#'   edge.
#' @param cex.main The amount by which the the plot title should be scaled
#'   relative to the default.
#' @param cex.lab The amount by which the the axis titles should be scaled
#'   relative to the default.
#' @param cex.axis The amount by which the the axis labels should be scaled
#'   relative to the default.
#' @param yRange The range of the y-axis.
#' @param lwd Line width.
#' @param lcol Line colour.
#' @param di Dimensions to be passed to \code{par(mfrow=())} to specify the size
#'   of the plotting window, when plotting multiple plots from a sars object of
#'   Type fit_collection. For example, \code{di = c(1, 3)} creates a plotting
#'   window with 1 row and 3 columns. The default (null) creates a square
#'   plotting window of the correct size.
#' @param pLeg Logical argument specifying whether or not the legend should be
#'   plotted for fit_collection plots (when \code{mfplot = TRUE}) or sar_multi
#'   plots. When a large number of model fits are plotted the legend takes up a
#'   lot of space, and thus the default is \code{pLeg = FALSE}.
#' @param \dots Further graphical parameters (see \code{\link[graphics]{par}},
#'   \code{\link[graphics]{plot}},\code{\link[graphics]{title}},
#'   \code{\link[graphics]{lines}}) may be supplied as arguments.
#' @importFrom graphics plot lines title
#' @examples
#' data(galap)
#' #fit and plot a sars object of Type fit.
#' fit <- sar_power(galap)
#' plot(fit, ModTitle = "A)", lcol = "blue")
#'
#' #fit and plot a sars object of Type fit_collection.
#' fit2 <- sar_expo(galap)
#' fit3 <- sar_epm1(galap)
#' fc <- fit_collection(fit, fit2, fit3)
#' plot(fc, ModTitle = letters[1:3], xlab = "Size of island")
#' @rdname plot.sars
#' @export



#xlab = NULL; ylab = NULL; pch = 16; cex = 1.2; 
#pcol = 'dodgerblue2'; ModTitle = NULL; TiAdj = 0; TiLine = 0.5; cex.main = 1.5;
#cex.lab = 1.3; cex.axis = 1; yRange = NULL;
#lwd = 2; lcol = 'dodgerblue2'; di = NULL




plot.sars <- function(x, mfplot = FALSE, xlab = NULL, ylab = NULL, pch = 16, cex = 1.2, 
                      pcol = 'dodgerblue2', ModTitle = NULL, TiAdj = 0, TiLine = 0.5, cex.main = 1.5,
                      cex.lab = 1.3, cex.axis = 1, yRange = NULL,
                      lwd = 2, lcol = 'dodgerblue2', di = NULL, pLeg = FALSE, ...)
{
  
  
  if (mfplot && attributes(x)$type != "fit_collection") stop("mfplot argument only for use with Type 'fit_collection'")
  
  
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
    if (is.null(yRange)){
    yMax <- max(c(yy,ff))#fitted line can be above the largest observed data point
    yMin <- 0
    yRange = c(yMin, yMax)
    }
    
    plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
         cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange, ...)
    title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main, ...)
    lines(x = xx, y = ff, lwd = lwd, col = lcol,  ...)
    
  }
  
  if (attributes(x)$type == "fit_collection"){
    
    if (!mfplot){

    if (!is.null(ModTitle)){
      if (length(ModTitle) == 1 && ModTitle == "") ModTitle <- rep("", length(x))
      if (length(ModTitle) != length(x)) stop("The length of ModTitle does not match the length of x")
      for (i in seq_along(x)){
        x[[i]]$model$name <- ModTitle[i]
      }
    }
    
    if (is.null(di)) {
      if (length(x) == 2){ #of length(x) = 2 the dividing by two does not work
        par(mfrow = c(1, 2))
      }else{
      di <- ceiling(length(x) / 2)
      par(mfrow = c(di, di))
      }#eo if x==2
    } else {
      par(mfrow = di)
    }#eo is.null if
      lapply(x, function(x){
      
      df <- x$data
      xx <- df$A
      yy <- df$S
      ff <- x$calculated
      ModTitle <- x$model$name 
      if (is.null(yRange)){
        yMax <- max(c(yy,ff))#fitted line can be above the largest observed data point
        yMin <- 0
        yRange = c(yMin, yMax)
      }
      
      plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
           cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,ylim = yRange, ...)
      title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main, ...)
      lines(x = xx, y = ff, lwd = lwd, col = lcol, ...)
    })
     par(mfrow = c(1,1))#change par back to default
    }# !mfplot
    
    if (mfplot){
      #observed data
      df <- x[[1]]$data 
      xx <- df$A
      yy <- df$S
      nams <-  vapply(x, function(x) x$model$name, FUN.VALUE = character(1))
      
      #fitted values for each model
      mf <- lapply(x, function(x) x$calculated)
      mf2 <- matrix(unlist(mf), ncol = length(x), byrow = FALSE)
      mf2 <- as.data.frame(mf2)
      colnames(mf2) <- nams
      ###
      if (is.null(yRange)){
        yMax <- max(c(yy,unlist(mf)))#fitted line can be above the largest observed data point
        yMin <- 0
        yRange = c(yMin, yMax)
      }
      
      #if legend to be included, work out size of plot
      if (pLeg == TRUE){
        xMax <- max(xx)*0.05
        lSiz <- legend(max(xx) +xMax, max(yy), legend = nams, horiz = F, lty = 1:ncol(mf2), col=1:ncol(mf2), plot = F)
        legWid <- lSiz$rect$left + lSiz$rect$w
        xMAX <- legWid + (legWid * 0.01)
        plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
             cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, xlim = c(min(xx), xMAX),
             ylim = yRange, bty = "L")
      } else {
      plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
           cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,ylim = yRange, bty = "L")
      }
      matlines(xx, mf2, lwd = lwd)
      title(main = "MultiModel Fits", adj = TiAdj, line = TiLine,cex.main = cex.main)
     if (pLeg == TRUE) legend(max(xx) + xMax, max(yy), legend = nams, horiz = F, lty = 1:ncol(mf2), 
                              col=1:ncol(mf2))
    }#eo mfplot
  
  }#eo if fit_collection
  
  if (attributes(x)$type == "lin_pow"){
    if (is.null(ModTitle)) ModTitle <- "Log-log power"
 
    df <- x$data
    xx <- df$A
    yy <- df$S
    ff <- x$calculated
    if (is.null(yRange)){
      yMax <- max(c(yy,ff))#fitted line can be above the largest observed data point
      yMin <- 0
      yRange = c(yMin, yMax)
    }
    
    plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
         cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,ylim = yRange, ...)
    title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main, ...)
    lines(x = xx, y = ff, lwd = lwd, col = lcol, ...)
  }
 }




#'  Plot Model Fits for a 'multi_sars' Object
#' @importFrom graphics plot lines title
#' @export

#NOT FINISHED

#need to add in ylim argument


plot.sar.multi <- function(x, type = "both", allCurves = TRUE,
                            xlab = NULL, ylab = NULL, pch = 16, cex = 1.2, 
                      pcol = 'dodgerblue2', ModTitle = NULL, TiAdj = 0, TiLine = 0.5, cex.main = 1.5,
                      cex.lab = 1.3, cex.axis = 1,
                      lwd = 2, lcol = 'dodgerblue2', di = c(1, 2), ...)
{
  ic <- x[[2]]$ic 
  dat <- x$details$fits
  
  #filter out bad models
  bad <- vapply(dat, function(x) any(is.na(x$sigConf)), FUN.VALUE = logical(1))
  dat2 <- dat[-which(bad)]
  
  #observed data
  df <- dat[[1]]$data 
  xx <- df$A
  yy <- df$S
  nams <-  vapply(dat2, function(x) x$model$name, FUN.VALUE = character(1))
  
  #fitted values for each model
  mf <- lapply(dat2, function(x) x$calculated)
  mf2 <- matrix(unlist(mf), ncol = length(dat2), byrow = FALSE)
  mf2 <- as.data.frame(mf2)
  colnames(mf2) <- nams
  
  #get correct IC info
  wh <- which(names(dat2[[1]]) == ic)
  icv <- vapply(dat2, function(x) unlist(x[wh]), FUN.VALUE = double(1))
  #delta
  delt <- icv - min(icv)
  #weight
  akaikesum <- sum(exp( -0.5*(delt)))
  aw <- exp(-0.5*delt) / akaikesum
  if (sum(aw) != 1) stop("IC weights do not sum to 1")
  
  #get weighted fitted values for each model
  mf3 <- matrix(nrow = nrow(mf2), ncol = ncol(mf2))
  for (i in seq_along(aw)) {mf3[ ,i] <- mf2[ ,i] * aw[i]}
  wfv <- rowSums(mf3)

  if (allCurves){
    mf2$MultiModel <- wfv
    nams2 <- c(nams, "MultiModel")
  }
  
  if (type == "both") par(mfrow = di)
  
  if (type == "both" || type == "multi"){
  #first plot with all curves
  if (allCurves){
      plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
      cex = cex, cex.lab = cex.lab, cex.axis = cex.axis)
      matlines(xx, mf2, lwd = lwd)
      legend("top", legend = nams2, inset=c(-0.2,0), lty = 1:ncol(mf2), col=1:ncol(mf2)) 
      title(main = "MultiModel Fits", adj = TiAdj, line = TiLine,cex.main = cex.main)
  } else if (!allCurves){
    #multimodel SAR curve
    plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
         cex = cex, cex.lab = cex.lab, cex.axis = cex.axis)
    title(main = "MultiModel Fits", adj = TiAdj, line = TiLine, cex.main = cex.main)
    lines(x = xx, y = wfv, lwd = lwd, col = lcol)
  }
  }
  
  if (type == "both" || type == "bar"){  
  ##barplot of IC weights
  barplot(aw, ylim=c(0, max(aw) + 0.05), cex.names=.68, ylab="IC weights", cex.lab = 1, 
          names.arg = nams)
  title(main = "Model weights", cex.main = 1.5, adj = 0, line = 0.5)
  }

}










