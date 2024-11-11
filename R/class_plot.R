#' Plot Model Fits for a 'sars' Object
#'
#' @description S3 method for class 'sars'. \code{plot.sars} creates plots
#'   for objects of class 'sars' (type = 'fit', "lin_pow' and
#'   'fit_collection'), using the R base plotting framework. The exact
#'   plot(s) constructed depends on the 'Type' attribute of the 'sars'
#'   object. For example, for a 'sars' object of Type 'fit', the
#'   \code{plot.sars} function returns a plot of the model fit (line) and the
#'   observed richness values (points). For a 'sars' object of Type
#'   'fit_collection' the \code{plot.sars} function returns either a grid
#'   with n individual plots (corresponding to the n model fits in the
#'   fit_collection), or a single plot with all n model fits included.
#'
#'   For plotting a 'sar_average' object, see \code{\link{plot.multi}}.
#' @param x An object of class 'sars'.
#' @param mfplot Logical argument specifying whether the model fits in a
#'   fit_collection should be plotted on one single plot (\code{mfplot =
#'   TRUE}) or separate plots (\code{mfplot = FALSE}; the default).
#' @param xlab Title for the x-axis (default depends on the Type attribute).
#' @param ylab Title for the y-axis (default depends on the Type attribute).
#' @param pch Plotting character (for points).
#' @param cex A numerical vector giving the amount by which plotting symbols
#'   (points) should be scaled relative to the default.
#' @param pcol Colour of the points.
#' @param ModTitle Plot title (default is \code{ModTitle = NULL}, which
#'   reverts to a default name depending on the type of plot). For no title,
#'   use \code{ModTitle = ""}. For a sars object of type fit_collection, a
#'   vector of names can be provided (e.g. \code{letters[1:3]}).
#' @param TiAdj Which way the plot title is justified.
#' @param TiLine Places the plot title this many lines outwards from the plot
#'   edge.
#' @param cex.main The amount by which the plot title should be scaled
#'   relative to the default.
#' @param cex.lab The amount by which the axis titles should be scaled
#'   relative to the default.
#' @param cex.axis The amount by which the axis labels should be scaled
#'   relative to the default.
#' @param yRange The range of the y-axis.
#' @param lwd Line width.
#' @param lcol Line colour.
#' @param di Dimensions to be passed to \code{par(mfrow=())} to specify the
#'   size of the plotting window, when plotting multiple plots from a sars
#'   object of Type fit_collection. For example, \code{di = c(1, 3)} creates
#'   a plotting window with 1 row and 3 columns. The default (null) creates a
#'   square plotting window of the correct size.
#' @param pLeg Logical argument specifying whether or not the legend should be
#'   plotted for fit_collection plots (when \code{mfplot = TRUE}) or. When a
#'   large number of model fits are plotted the legend takes up a lot of space,
#'   and thus the default is \code{pLeg = FALSE}.
#' @param \dots Further graphical parameters (see
#'   \code{\link[graphics]{par}},
#'   \code{\link[graphics]{plot.default}},\code{\link[graphics]{title}},
#'   \code{\link[graphics]{lines}}) may be supplied as arguments.
#' @importFrom graphics par plot legend barplot
#' @importFrom graphics points lines polygon title matlines matplot
#' @examples
#' data(galap)
#' #fit and plot a sars object of Type fit.
#' fit <- sar_power(galap)
#' plot(fit, ModTitle = "A)", lcol = "blue")
#'
#' #fit and plot a sars object of Type fit_collection.
#' fc <- sar_multi(data = galap, obj = c("power", "loga", "epm1"), 
#' grid_start = "none")
#' plot(fc, ModTitle = letters[1:3], xlab = "Size of island")
#' @rdname plot.sars
#' @export


plot.sars <- function(x, mfplot = FALSE, xlab = NULL, ylab = NULL,
                      pch = 16, cex = 1.2,
                      pcol = 'dodgerblue2', ModTitle = NULL, TiAdj = 0,
                      TiLine = 0.5, cex.main = 1.5,
                      cex.lab = 1.3, cex.axis = 1, yRange = NULL,
                      lwd = 2, lcol = 'dodgerblue2', di = NULL,
                      pLeg = FALSE, ...)
{

  if (attributes(x)$type == "pred"){
    return(cat("\nNo plot method for 'pred' object of class 'sars'\n", 
               sep = ""))
  }
  if (attributes(x)$type == "threshold_ci"){
    return(cat("\nNo plot method for 'threshold_ci' object of class 'sars'\n", 
               sep = ""))
  }
  if (attributes(x)$type == "threshold_coef"){
    return(cat("\nNo plot method for 'threshold_coef' object of class 'sars'\n", 
               sep = ""))
  }

  if (mfplot & attributes(x)$type != "fit_collection")
    stop("mfplot argument only for use with Type 'fit_collection'")


  if (is.null(xlab)){
    if (attributes(x)$type == "fit" |
        attributes(x)$type == "fit_collection"){
        xlab <- "Area"
    } else if (attributes(x)$type == "lin_pow"){
      xlab <- "Log(Area)"
    } else {
      stop("Type attribute not recognised")
    }
  }

  if (is.null(ylab)){
    if (attributes(x)$type == "fit" |
        attributes(x)$type == "fit_collection"){
      ylab <- "Species richness"
    } else if (attributes(x)$type == "lin_pow"){
      ylab <- "Log(Species richness)"
    }
  }

  if (attributes(x)$type == "fit"){
    if (is.null(ModTitle)) ModTitle <- x$model$name

    df <- x$data
    xx <- df$A
    yy <- df$S
    xx2 <- seq(min(xx), max(xx), length.out = 1000)
    if (x$model$name == "Linear model"){
      c1 <- x$par[1]
      m <- x$par[2]
      ff <- c1 + (m * xx2)#have to have 1000 points for linear so it matches
                               #with xx2 for plotting
    } else {
      #create a set of 1000 fitted points using fitted model parameters
      #to ensure a smooth curve
      if (!all(x$model$mod.fun(xx, x$par) == x$calculated)){
        stop("Error in plotting, contact package owner")
      }
      ff <- x$model$mod.fun(xx2, x$par)
      if (anyNA(ff)){
        stop("Error in plotting, contact package owner")
      }
    }
    if (is.null(yRange)){
      yMax <- max(c(yy,ff))#fitted line can be above the largest observed point
      yMin <- min(c(yy,ff))
      yRange <- c(yMin, yMax)
    }

    plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol,
         cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,
         ylim = yRange, ...)
    title(main = ModTitle, adj = TiAdj, line = TiLine,
          cex.main = cex.main, ...)
    lines(x = xx2, y = ff, lwd = lwd, col = lcol,  ...)

  }

  if (attributes(x)$type == "fit_collection"){

    if (!mfplot){

    if (!is.null(ModTitle)){
      if (length(ModTitle) == 1 && ModTitle == "")
        ModTitle <- rep("", length(x))
      if (length(ModTitle) != length(x))
        stop("The length of ModTitle does not match the length of x")
      for (i in seq_along(x)){
        x[[i]]$model$name <- ModTitle[i]
      }
    }

    if (is.null(di)) {
      if (length(x) == 2){ #of length(x) = 2 the dividing by two doesnt work
        par(mfrow = c(1, 2))
      } else{
      di <- ceiling(length(x) / 2)
      par(mfrow = c(di, di))
      }#eo if x==2
    } else {
      par(mfrow = di)
    }#eo is.null if
      
      df <- x[[1]]$data
      xx <- df$A
      yy <- df$S
      xx2 <- seq(min(xx), max(xx), length.out = 1000)
      
      lapply(x, function(y){
      if (y$model$name == "Linear model"){
        c1 <- y$par[1]
        m <- y$par[2]
        ff <- c1 + (m * xx2)#have to have 1000 points for linear so it matches
                            #with xx2 for plotting
      } else {
        #create a set of 1000 fitted points using fitted model parameters
        #to ensure a smooth curve
        if (!all(y$model$mod.fun(xx, y$par) == y$calculated)){
          stop("Error in plotting, contact package owner")
        }
        ff <- y$model$mod.fun(xx2, y$par)
        if (anyNA(ff)){
          stop("Error in plotting, contact package owner")
        }
      }
      
      ModTitle <- x$model$name
      if (is.null(yRange)){
        yMax <- max(c(yy,ff))#fitted line can be above the largest observed
        yMin <- min(c(yy,ff))
        yRange <- c(yMin, yMax)
      }

      plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol,
           cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,
           ylim = yRange, ...)
      title(main = ModTitle, adj = TiAdj, line = TiLine,
            cex.main = cex.main, ...)
      lines(x = xx2, y = ff, lwd = lwd, col = lcol, ...)
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
     # mf <- lapply(x, function(x) x$calculated)
      xx2 <- seq(min(xx), max(xx), length.out = 1000)
      mf <- lapply(x, function(y) {
        if (y$model$name == "Linear model"){
          c1 <- y$par[1]
          m <- y$par[2]
          ff <- c1 + (m * xx2)#have to have 1000 points for linear so it fits
                              #in mf2 matrix
        } else {
          #create a set of 1000 fitted points using fitted model parameters
          #to ensure a smooth curve
          ff <- y$model$mod.fun(xx2, y$par)
          if (anyNA(ff)){
            stop("Error in plotting, contact package owner")
          }
          #ff <- x$calculated
        }
        ff
        })#eo lapply
      mf2 <- matrix(unlist(mf), ncol = length(x), byrow = FALSE)
      mf2 <- as.data.frame(mf2)
      colnames(mf2) <- nams
      ###
      if (is.null(yRange)){
        yMax <- max(c(yy,unlist(mf)))
        yMin <- min(c(yy,unlist(mf)))
        yRange <- c(yMin, yMax)
      }

      #main title
      if (is.null(ModTitle)) ModTitle <- ""

      #if legend to be included, work out size of plot
      if (pLeg == TRUE){
        #xMax <- max(xx)*0.05
        #lSiz <- legend(max(xx) +xMax, max(yy), legend = nams,
        #horiz = F, lty = 1:ncol(mf2), col=1:ncol(mf2), plot = F)
        #legWid <- lSiz$rect$left + lSiz$rect$w
        xMAX <- max(xx) + max(xx) * 0.5
        matplot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch,
                col = pcol,
             cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,
             xlim = c(min(xx), xMAX),
             ylim = yRange, bty = "L")
      } else {
        plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol,
           cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,
           ylim = yRange, bty = "L")
      }
      matlines(xx2, mf2, lwd = lwd, lty = seq_along(mf2), col=seq_along(mf2))
      title(main = ModTitle, adj = TiAdj, line = TiLine,cex.main = cex.main)
     if (pLeg == TRUE) legend(max(xx) + (max(xx) * 0.05), yMax,
                              legend = nams, horiz = FALSE,
                              lty = seq_along(mf2),col=seq_along(mf2))
    }#eo mfplot

  }#eo if fit_collection

  if (attributes(x)$type == "lin_pow"){
    if (is.null(ModTitle)) ModTitle <- "Log-log power"

    df <- x$data
    xx <- df$A
    yy <- df$S
    ff <- x$calculated
    if (is.null(yRange)){
      yMax <- max(c(yy,ff))#
      yMin <- min(c(yy,ff))
      yRange <- c(yMin, yMax)
    }

    plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch,
         col = pcol,
         cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,
         ylim = yRange, ...)
    title(main = ModTitle, adj = TiAdj, line = TiLine,
          cex.main = cex.main, ...)
    lines(x = xx, y = ff, lwd = lwd, col = lcol, ...)
  }
 }



#' Plot Model Fits for a 'multi' Object
#'
#' @description S3 method for class 'multi'. \code{plot.multi} creates plots
#'   for objects of class multi, using the R base plotting framework. Plots
#'   of all model fits, the multimodel SAR curve (with confidence intervals)
#'   and a barplot of the information criterion weights of the different
#'   models can be constructed.
#' @param x An object of class 'multi'.
#' @param type The type of plot to be constructed: either \code{type = multi}
#'   for a plot of the multimodel SAR curve, or \code{type = bar} for a
#'   barplot of the information criterion weights of each model.
#' @param allCurves A logical argument for use with \code{type = multi} that
#'   specifies whether all the model fits should be plotted with the
#'   multimodel SAR curve (\code{allCurves = TRUE}; the default) or that only
#'   the multimodel SAR curve should be plotted (\code{allCurves = FALSE}).
#' @param xlab Title for the x-axis. Only for use with \code{type = multi}.
#' @param ylab Title for the y-axis.
#' @param pch Plotting character (for points). Only for use with \code{type =
#'   multi}.
#' @param cex A numerical vector giving the amount by which plotting symbols
#'   (points) should be scaled relative to the default.
#' @param pcol Colour of the points. Only for use with \code{type = multi}.
#' @param ModTitle Plot title (default is \code{ModTitle = NULL}, which
#'   reverts to "Multimodel SAR" for \code{type = multi} and to "Model
#'   weights" for \code{type = bar}). For no title, use \code{ModTitle = ""}.
#' @param TiAdj Which way the plot title is justified.
#' @param TiLine Places the plot title this many lines outwards from the plot
#'   edge.
#' @param cex.main The amount by which the plot title should be scaled
#'   relative to the default.
#' @param cex.lab The amount by which the axis titles should be scaled
#'   relative to the default.
#' @param cex.axis The amount by which the axis labels should be scaled
#' relative to the default.
#' @param yRange The range of the y-axis. Only for use with \code{type = multi}.
#' @param lwd Line width. Only for use with \code{type = multi}.
#' @param lcol Line colour. Only for use with \code{type = multi}.
#' @param mmSep Logical argument of whether the multimodel curve should be
#'   plotted as a separate line (default = FALSE) on top of the others, giving
#'   the user more control over line width and colour. Only for use with
#'   \code{type = multi} and \code{allCurves = TRUE}.
#' @param lwd.Sep If \code{mmSep = TRUE}, the line width of the multimodel
#'   curve.
#' @param col.Sep If \code{mmSep = TRUE}, the colour of the multimodel curve.
#' @param pLeg Logical argument specifying whether or not the legend should be
#'   plotted  (when \code{type = multi} and \code{allCurves = TRUE}).
#' @param modNames A vector of model names for the barplot of weights (when
#'   \code{type = bar}). The default (\code{modNames = NULL}) uses abbreviated
#'   versions (see below) of the names from the \code{sar_average} function.
#' @param cex.names The amount by which the axis labels (model names) should be
#'   scaled relative to the default. Only for use with \code{type = bar}.
#' @param subset_weights Only create a barplot of the model weights for models
#'   with a weight value above a given threshold (\code{subset_weights}). Only
#'   for use with \code{type = bar}.
#' @param confInt A logical argument specifying whether confidence intervals
#'   should be plotted around the multimodel curve. Can only be used if
#'   confidence intervals have been generated in the \code{sar_average}
#'   function.
#' @param \dots Further graphical parameters (see
#'   \code{\link[graphics]{par}},
#'   \code{\link[graphics]{plot.default}},\code{\link[graphics]{title}},
#'   \code{\link[graphics]{lines}}) may be supplied as arguments.
#' @note In some versions of R and R studio, when plotting all model fits on the
#'   same plot with a legend it is necessary to manually extend your plotting
#'   window (height and width; e.g. the 'Plots' window of R studio) before
#'   plotting to ensure the legend fits in the plot. Extending the plotting
#'   window after plotting sometimes just stretches the legend.
#'
#'   Occasionally a model fit will converge and pass the model fitting checks
#'   (e.g. residual normality) but the resulting fit is nonsensical (e.g. a
#'   horizontal line with intercept at zero). Thus, it can be useful to plot
#'   the resultant 'multi' object to check the individual model fits. To
#'   re-run the \code{sar_average} function without a particular model, simply
#'   remove it from the \code{obj} argument.
#'
#'   For visual interpretation of the model weights barplot it is necessary
#'   to abbreviate the model names when plotting the weights of several
#'   models. To plot fewer bars, use the \code{subset_weights} argument to
#'   filter out models with lower weights than a threshold value. To provide
#'   a different set of names use the \code{modNames} argument. The model
#'   abbreviations used as the default are: \itemize{ 
#'   \item \strong{Pow} =   Power
#'   \item \strong{PowR} =   PowerR 
#'   \item \strong{E1} =   Extended_Power_model_1 
#'   \item \strong{E2} =   Extended_Power_model_2 
#'   \item \strong{P1} =   Persistence_function_1
#'   \item \strong{P2} =   Persistence_function_2 
#'   \item \strong{Loga} =   Logarithmic
#'   \item \strong{Kob} =   Kobayashi 
#'   \item \strong{MMF} =   MMF 
#'   \item \strong{Mon} =   Monod
#'   \item \strong{NegE} =   Negative_exponential 
#'   \item \strong{CR} =   Chapman_Richards
#'   \item \strong{CW3} =   Cumulative_Weibull_3_par. 
#'   \item \strong{AR} =  Asymptotic_regression 
#'   \item \strong{RF} =   Rational_function 
#'   \item \strong{Gom} =  Gompertz 
#'   \item \strong{CW4} =   Cumulative_Weibull_4_par. 
#'   \item \strong{BP} =  Beta-P_cumulative 
#'   \item \strong{Logi} =   Logistic(Standard)
#'   \item \strong{Hel} =   Heleg(Logistic) 
#'   \item \strong{Lin} =   Linear_model}
#' @examples
#' data(galap)
#' #plot a multimodel SAR curve with all model fits included
#' fit <- sar_average(data = galap, grid_start = "none")
#' plot(fit)
#'
#' #remove the legend
#' plot(fit, pLeg = FALSE)
#'
#' #plot just the multimodel curve
#' plot(fit, allCurves = FALSE, ModTitle = "", lcol = "black")
#' 
#' #plot all model fits and the multimodel curve on top as a thicker line
#' plot(fit, allCurves = TRUE, mmSep = TRUE, lwd.Sep = 6, col.Sep = "orange")
#'
#' #Plot a barplot of the model weights
#' plot(fit, type = "bar")
#' #subset to plot only models with weight > 0.05
#' plot(fit, type = "bar", subset_weights = 0.05)
#' @rdname plot.multi
#' @export


plot.multi <- function(x, type = "multi", allCurves = TRUE,
                            xlab = NULL, ylab = NULL, pch = 16, cex = 1.2,
                      pcol = 'dodgerblue2', ModTitle = NULL,
                      TiAdj = 0, TiLine = 0.5, cex.main = 1.5,
                      cex.lab = 1.3, cex.axis = 1, yRange = NULL,
                      lwd = 2, lcol = 'dodgerblue2', mmSep = FALSE,
                      lwd.Sep = 6, col.Sep = "black", pLeg = TRUE,
                      modNames = NULL, cex.names=.88,
                      subset_weights = NULL, confInt = FALSE, ...)
{

  if (confInt){
    if (length(x$details$confInt) == 1)
      stop("No confidence interval information in the fit object")
    CI <- x$details$confInt
  }


  ic <- x[[2]]$ic
  dat <- x$details$fits

  #filter out bad models
  #bad <- vapply(dat, function(x) any(is.na(x$sigConf)),
  #FUN.VALUE = logical(1))
  #dat2 <- dat[-which(bad)]
  dat2 <- dat

  #observed data
  df <- dat[[1]]$data
  xx <- df$A
  yy <- df$S
  nams <-  vapply(dat2, function(x) x$model$name, FUN.VALUE = character(1))
  xx2 <- seq(min(xx), max(xx), length.out = 1000)

  #fitted values for each model
  mf <- lapply(dat2, function(y) {
    if (y$model$name == "Linear model"){
      c1 <- y$par[1]
      m <- y$par[2]
      ff <- c1 + (m * xx2)#have to have 1000 points for linear so it fits
      #in mf2 matrix
    } else {
      #create a set of 1000 fitted points using fitted model parameters
      #to ensure a smooth curve
      ff <- y$model$mod.fun(xx2, y$par)
      if (anyNA(ff)){
        stop("Error in plotting, contact package owner")
      }
      #ff <- x$calculated
    }
    ff
  })
  
  mf2 <- matrix(unlist(mf), ncol = length(dat2), byrow = FALSE)
  mf2 <- as.data.frame(mf2)
  colnames(mf2) <- nams

  #get correct IC info
  icv <- vapply(dat2, function(x) unlist(x[[ic]]), FUN.VALUE = double(1))
  #delta
  delt <- icv - min(icv)
  #weight
  akaikesum <- sum(exp( -0.5*(delt)))
  aw <- exp(-0.5*delt) / akaikesum
  if (round(sum(aw), 0) != 1) stop("IC weights do not sum to 1")#have to
  #round as sometimes fractionally different to 1

  #get weighted fitted values for each model
  mf3 <- matrix(nrow = nrow(mf2), ncol = ncol(mf2))
  for (i in seq_along(aw)) mf3[ ,i] <- mf2[ ,i] * aw[i]
  wfv <- rowSums(mf3)

  #this is a test error for development
 # if (!all(round(wfv) == round(x$mmi)))
   # stop("Multimodel fitted values do not match between functions")

  if (allCurves & (!mmSep)){
    mf2$MultiModel <- wfv
    nams2 <- c(nams, "Multimodel SAR")
  } else{
    nams2 <- nams
  }

  if (type == "multi"){

    #set axis names
    if (is.null(xlab)) xlab <- "Area"
    if (is.null(ylab)) ylab <- "Species richness"


    #set y axis range
    if (is.null(yRange)){
      if (allCurves){
        if (confInt){ #CIs larger so need to add to ymax and ymin
          yMax <- max(c(yy,unlist(mf2), CI$U))
          yMin <- min(c(yy,unlist(mf2), CI$L))
        } else {
        yMax <- max(c(yy,unlist(mf2)))
        yMin <- min(c(yy,unlist(mf2)))
      }
        }else{
        if (confInt){ #CIs larger so need to add to ymax and ymin
          yMax <- max(c(yy,wfv, CI$U))
          yMin <- min(c(yy,wfv, CI$L))
        } else {
        yMax <- max(c(yy,wfv))
        yMin <- min(c(yy,wfv))
        }
      }
      yRange <- c(yMin, yMax)
    }

  #main title
  if (is.null(ModTitle)) ModTitle <- "Multimodel SAR"

  #first plot with all curves
  if (allCurves){
    #if legend to be included, work out size of plot
    if (pLeg == TRUE){
     # xMax <- max(xx)*0.05
      #plot(x = xx, y = yy, xlim = c(min(xx), xMAX), ylim = yRange)
     # lSiz <- legend(max(xx) + xMax, max(c(yy,unlist(mf2))),
      #legend = nams2, horiz = F, lty = 1:ncol(mf2),
      #col=1:ncol(mf2), plot = F)
     # legWid <- lSiz$rect$left + lSiz$rect$w
      #xMAX <- legWid + (legWid * 0.01)

      xMAX <- max(xx) + max(xx) * 0.5

     # legHeight <- lSiz$rect$h
    #  yMAX <- legHeight + (legHeight * 0.3)
     # yRange <- c(yMin, yMAX)

      if (confInt){
        matplot(x = xx, y = yy, xlab = xlab, ylab = ylab,
                cex.lab = cex.lab, cex.axis = cex.axis,
                xlim = c(min(xx), xMAX), ylim = yRange, pch = pch,
                col = "white")
        polygon(c(xx,rev(xx)),c(CI$L,rev(CI$U)),col="grey87",border=NA)
        points(x = xx, y = yy, pch = pch, col = pcol,
               cex = cex)
      } else {
        matplot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch,
                col = pcol,
                cex = cex, cex.lab = cex.lab, cex.axis = cex.axis,
                xlim = c(min(xx), xMAX), ylim = yRange)
      }#eo confInt
    } else{ #no legend

      if (confInt){
        plot(x = xx, y = yy, xlab = xlab, ylab = ylab,
             cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange)
        polygon(c(xx,rev(xx)),c(CI$L,rev(CI$U)),col="grey87",border=NA)
        points(x = xx, y = yy, pch = pch, col = pcol,
               cex = cex)
      } else {
        plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol,
        cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange)
      }#eo confInt
    }#eo no legend
      matlines(xx2, mf2, lwd = lwd, lty = seq_along(mf2), col=seq_along(mf2))
      title(main = ModTitle, adj = TiAdj, line = TiLine,cex.main = cex.main)
      if (pLeg){
        if (!mmSep) {
          legend(max(xx) + (max(xx) * 0.05), yMax,
                 legend = nams2, horiz = FALSE,
                 lty = seq_along(mf2), col= seq_along(mf2))
        } else {
          legend(max(xx) + (max(xx) * 0.05), yMax,
                 legend = c(nams2, "Multimodel SAR"), horiz = FALSE,
                 lty = c(seq_along(mf2), 1), col=c(seq_along(mf2), col.Sep),
                 lwd = c(rep(lwd, length(nams2)), lwd.Sep))
        }
      }#eo if PLeg
      if (mmSep) lines(xx2, wfv, lwd = lwd.Sep, col = col.Sep)
  } else if (!allCurves){
    #just multimodel SAR curve
    if (confInt){
      plot(x = xx, y = yy, xlab = xlab, ylab = ylab,
            cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange)
      polygon(c(xx,rev(xx)),c(CI$L,rev(CI$U)),col="grey87",border=NA)
      points(x = xx, y = yy, pch = pch, col = pcol,
             cex = cex)
      title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main)
      lines(x = xx2, y = wfv, lwd = lwd, col = lcol)
    } else {
      plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol,
      cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange)
      title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main)
      lines(x = xx2, y = wfv, lwd = lwd, col = lcol)
    }#eo confint
  }
  }

  if (type == "bar"){
  ##barplot of IC weight

  if (!is.null(subset_weights)) aw <- aw[aw > subset_weights]

  if (is.null(ylab)) ylab <- "IC weights"
  if (is.null(ModTitle)) ModTitle <- "Model weights"
  if (is.null(modNames)){
    modNames <- names(aw)
    modNames <- mod_abbrev(modNames)
  }

  barplot(aw, ylim=c(0, max(aw) + 0.05), cex.names= cex.names,
          ylab = ylab, cex.lab = cex.lab,
          cex.axis = cex.axis, names.arg = modNames)
  title(main = ModTitle, cex.main = cex.main, adj = TiAdj, line = TiLine)
  }

}

#########################################################################
#############Threshold plots###############################################
####################################################################

##############internal plotting functions

#' Discontinuous one thr plot
#'@importFrom graphics par plot lines title
#'@noRd

discOnePlot <- function (xx, yy, multPlot, xypred.disc, data, th, xlab = xlab, 
                         ylab = ylab, pch = pch, 
                         pcol = pcol, cex = cex, cex.lab = cex.lab, 
                         cex.axis = cex.axis, 
                         yRange = yRange, ModTitle = ModTitle, TiAdj = TiAdj,
                         TiLine = TiLine, 
                         cex.main = cex.main, lwd = lwd, lcol = lcol, ...) 
{
  disc1 <- xypred.disc[xypred.disc$x <= th["DiscOne"][1], ]
  disc2 <- xypred.disc[xypred.disc$x > th["DiscOne"][1], ]
  if (multPlot){
    plot(xx, yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
       cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange, 
       ...)
    title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main, 
        ...)
  }
  lines(disc1$x, disc1$DiscOne, lwd = lwd, col = lcol, ...)
  lines(disc2$x, disc2$DiscOne, lwd = lwd, col = lcol, ...)
}


#' Discontinuous two thr plot
#'@importFrom graphics par plot lines title
#'@noRd

discTwoPlot <- function (xx, yy, multPlot, xypred.disc, data, th, xlab = xlab, 
                         ylab = ylab, pch = pch, 
                         pcol = pcol, cex = cex, cex.lab = cex.lab, 
                         cex.axis = cex.axis, 
                         yRange = yRange, ModTitle = ModTitle, TiAdj = TiAdj, 
                         TiLine = TiLine, 
                         cex.main = cex.main, lwd = lwd, lcol = lcol, ...) 
{
  disc1 <- xypred.disc[xypred.disc$x <= th["DiscTwo"][[1]][1], ]
  disc2 <- xypred.disc[xypred.disc$x > th["DiscTwo"][[1]][1] & 
                         xypred.disc$x <= th["DiscTwo"][[1]][2], ]
  disc3 <- xypred.disc[xypred.disc$x > th["DiscTwo"][[1]][2], ]
  if (multPlot){
    plot(xx, yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
       cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange, 
       ...)
    title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main, 
        ...)
  }
  lines(disc1$x, disc1$DiscTwo, lwd = lwd, col = lcol, ...)
  lines(disc2$x, disc2$DiscTwo, lwd = lwd, col = lcol, ...)
  lines(disc3$x, disc3$DiscTwo, lwd = lwd, col = lcol, ...)
}


#' continuous models plot
#'@importFrom graphics par plot lines title
#'@noRd

contsPlot <- function (xx, yy, multPlot, xypred.cont, data, column, xlab = xlab, 
                       ylab = ylab, pch = pch, 
                       pcol = pcol, cex = cex, cex.lab = cex.lab,
                       cex.axis = cex.axis, 
                       yRange = yRange, ModTitle = ModTitle, TiAdj = TiAdj, 
                       TiLine = TiLine, 
                       cex.main = cex.main, lwd = lwd, lcol = lcol, ...) 
{
  if (multPlot){
    plot(xx, yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
       cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange, 
       ...)
    title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main, 
        ...)
  }
  lines(xypred.cont$x, xypred.cont[column][, 1], lwd = lwd, col = lcol, ...)
}


#'Plot Model Fits for a 'threshold' Object
#'
#'@description S3 method for class 'threshold'. \code{plot.threshold} creates
#'  plots for objects of class threshold, using the R base plotting framework.
#'  Plots of single or multiple threshold models can be constructed.
#'@param x An object of class 'threshold'.
#'@param xlab Title for the x-axis. Defaults will depend on any axes
#'  log-transformations.
#'@param ylab Title for the y-axis.Defaults will depend on any axes
#'  log-transformations.
#'@param multPlot Whether separate plots should be built for each model fit
#'  (default = TRUE) or all model fits should be printed on the same plot
#'  (FALSE)
#'@param pch Plotting character (for points).
#'@param cex A numerical vector giving the amount by which plotting symbols
#'  (points) should be scaled relative to the default.
#'@param pcol Colour of the points.
#'@param ModTitle Plot title (default is \code{ModTitle = NULL}), which reverts
#'  to the model names. For no title, use \code{ModTitle = ""}.
#'@param TiAdj Which way the plot title is justified.
#'@param TiLine Places the plot title this many lines outwards from the plot
#'  edge.
#'@param cex.main The amount by which the plot title should be scaled relative
#'  to the default.
#'@param cex.lab The amount by which the axis titles should be scaled relative
#'  to the default.
#'@param cex.axis The amount by which the axis labels should be scaled relative
#'  to the default.
#'@param yRange The range of the y-axis. Default taken as the largest value
#'  bacross the observed and fitted values.
#'@param lwd Line width.
#'@param lcol Line colour. If \code{multPlot = TRUE}, just a single colour
#'  should be given, If  \code{multPlot = FALSE}, either a single colour, or a
#'  vector of colours the same length as the number of model fits in \code{x}.
#'@param di Dimensions to be passed to \code{par(mfrow=())} to specify the size
#'  of the plotting window, when plotting multiple plots. For example, \code{di
#'  = c(1, 3)} creates a plotting window with 1 row and 3 columns. The default
#'  (\code{NULL}) creates a plotting window large enough to fit all plots in.
#'@param \dots Further graphical parameters (see \code{\link[graphics]{par}},
#'  \code{\link[graphics]{plot.default}},\code{\link[graphics]{title}},
#'  \code{\link[graphics]{lines}}) may be supplied as arguments.
#'@note The raw \code{\link{lm}} model fit objects are returned with the
#'  \code{\link{sar_threshold}} function if the user wishes to construct their
#'  own plots.
#'  
#'  Use \code{par(mai = c())} prior to calling plot, to set the graph margins,
#'  which can be useful when plotting multiple models in a single plot to ensure
#'  space within the plot taken up by the individual model fit plots is
#'  maximised.
#' @examples
#' data(aegean)
#'
#'#fit two threshold models (in logA-S space) and the linear and 
#'#intercept only models
#'fct <- sar_threshold(aegean, mod = c("ContOne", "DiscOne"),
#'                     non_th_models = TRUE,  interval = 5,
#'                     parallel = FALSE, logAxes = "area")
#'
#'#plot using default settings
#'plot(fct)
#'
#'#change various plotting settings, and set the graph margins prior to
#'#plotting
#'par(mai = c(0.7,0.7, 0.4, 0.3))
#'plot(fct, pcol = "blue", pch = 18, lcol = "green",
#'     ModTitle = c("A", "B", "C", "D"), TiAdj = 0.5, xlab = "Yorke")
#'     
#'#Plot multiple model fits in the same plot, with different colour for each 
#'#model fit
#'plot(fct, multPlot = FALSE, lcol = c("yellow", "red", "green", "purple"))
#' #Making a legend. First extract the model order:
#' fct[[2]]
#' #Then use the legend function -  note you may need to play around with the 
#' #legend location depending on your data.
#' legend("topleft", legend=c("ContOne", "DiscOne","Linear", "Intercept"), fill =
#' c("yellow", "red", "green", "purple"))
#'@rdname plot.threshold
#'@importFrom stats predict
#'@importFrom graphics par
#'@export


plot.threshold <- function(x, xlab = NULL, ylab = NULL, multPlot = TRUE,
                           pch = 16, cex = 1.2,
                           pcol = 'black', ModTitle = NULL, TiAdj = 0,
                           TiLine = 0.5, cex.main = 1.5,
                           cex.lab = 1.3, cex.axis = 1, yRange = NULL,
                           lwd = 2, lcol = 'red', di = NULL, ...) {
  
  if (length(lcol) > 1){
    if (length(lcol) != length(x[[1]])){
      stop("lcol should be a single colour or a vector equal to no. of models")
    }
  }
  
  if (!is.logical(multPlot)) stop("multPlot should be logical")
  
  data <- x[[4]]#nb. this will already be log-transformed if user has selected
  colnames(data) <- c("A", "S")
  data <- data[order(data$A),]
  mods <- x[[1]]
  names <- x[[2]]
  th <- x[[3]]
  names(mods) <- names
  names(th) <- names
  yy <- data$S
  xx <- data$A
  
  cont <- c("ContOne","ZslopeOne","ContTwo","ZslopeTwo","Linear","Intercept")
  xypred.cont <- data.frame("x" = seq(min(data$A), max(data$A), 
                                      (max(data$A) - min(data$A))/1000))
  for (i in which(names %in% cont)) {
    xypred.cont[, names[i]] <- stats::predict(mods[[i]], xypred.cont)
  }
  
  disc <- c("DiscOne","DiscTwo")
  xypred.disc <- data.frame("x" = xx)
  for (i in which(names %in% disc)) {
    xypred.disc[, names[i]] <- stats::predict(mods[[i]])
  }
  
  ##get max y-value across model fits and observed to 
  #set y-axis range (unless provided)
  if (is.null(yRange)){
    all_values <- data$S #start with observed richness values
    #then add in any predicted values from cont or disc models (if they are included)
    #The >1 argument is because the first column is just the area values
    if (ncol(xypred.cont) > 1) all_values <- c(all_values, 
                                               unlist(xypred.cont[,2:ncol(xypred.cont)]))
    if (ncol(xypred.disc) > 1) all_values <- c(all_values, 
                                               unlist(xypred.disc[,2:ncol(xypred.disc)]))
    ff <- range(all_values)
    yRange <- c(ff[1], ff[2])
  }
  if (is.null(xlab)) {
    xlab <- switch(x[[5]][[1]], 
                   none = "Area", 
                   area = "Log(Area)", 
                   both = "Log(Area)")
  }
  if (is.null(ylab)) {
    ylab <- switch(x[[5]][[1]], 
                   none = "Species richness", 
                   area = "Species richness", 
                   both = "Log(Species richness)")
  }
  if (is.null(ModTitle)) {
    if (multPlot){
    ModTitle <- sapply(names, function(x) {
      switch(x, ContOne = "Continuous one-threshold", 
             ZslopeOne = "Left-horizontal one-threshold", 
             DiscOne = "Discontinuous one-threshold", 
             ContTwo = "Continuous two-threshold", 
             ZslopeTwo = "Left-horizontal two-threshold", 
             DiscTwo = "Discontinuous two-threshold", 
             Linear = "Linear", 
             Intercept = "Intercept only")
    })
    }
  }
  
  ###############################################
  ##if only one model in object##################
  #################################################
  if (length(x[[1]]) == 1) {
    if (!multPlot) multPlot <- TRUE #no sense being FALSE if only one model
    if (any(c("DiscOne", "DiscTwo") %in% names)) {
      if (names == "DiscOne") {
        discOnePlot(xx, yy, multPlot = multPlot,
                    xypred.disc, data, th, xlab = xlab, ylab = ylab, 
                    pch = pch, pcol = pcol, cex = cex, cex.lab = cex.lab, 
                    cex.axis = cex.axis, yRange = yRange, ModTitle = ModTitle, 
                    TiAdj = TiAdj, TiLine = TiLine, cex.main = cex.main, 
                    lwd = lwd, lcol = lcol, ...)
      }
      else {
        discTwoPlot(xx, yy, multPlot = multPlot,
                    xypred.disc, data, th, xlab = xlab, ylab = ylab, 
                    pch = pch, pcol = pcol, cex = cex, cex.lab = cex.lab, 
                    cex.axis = cex.axis, yRange = yRange, ModTitle = ModTitle, 
                    TiAdj = TiAdj, TiLine = TiLine, cex.main = cex.main, 
                    lwd = lwd, lcol = lcol, ...)
      }
    }
    else {
      contsPlot(xx, yy, multPlot = multPlot,
                xypred.cont, data, column = names, xlab = xlab, 
                ylab = ylab, pch = pch, pcol = pcol, cex = cex, 
                cex.lab = cex.lab, cex.axis = cex.axis, yRange = yRange, 
                ModTitle = ModTitle, TiAdj = TiAdj, TiLine = TiLine, 
                cex.main = cex.main, lwd = lwd, lcol = lcol, 
                ...)
    }
  } else {
    if (multPlot){
    ##############################################################
    ##if more than one model###################################
    #######################################################
    if (is.null(di)) {
      if (length(mods) == 2) {
        par(mfrow = c(1, 2))
      }
      else {
        lmn <- as.character(length(mods))
        di <- switch(lmn, 
                     `3` = c(1, 3), 
                     `4` = c(2, 2), 
                     `5` = c(2, 3), 
                     `6` = c(2, 3), 
                     `7` = c(3, 3), 
                     `8` = c(3, 3))
        par(mfrow = c(di[1], di[2]))
      }
    }
    else {
      par(mfrow = di)
    }
    for (i in 1:length(mods)) {
      if (any(c("DiscOne", "DiscTwo") %in% names[[i]])) {
        if (names[[i]] == "DiscOne") {
          discOnePlot(xx, yy, multPlot = multPlot,
                      xypred.disc, data, th[i], xlab = xlab, 
                      ylab = ylab, pch = pch, pcol = pcol, cex = cex, 
                      cex.lab = cex.lab, cex.axis = cex.axis, yRange = yRange, 
                      ModTitle = ModTitle[i], TiAdj = TiAdj, TiLine = TiLine, 
                      cex.main = cex.main, lwd = lwd, lcol = lcol, 
                      ...)
        }
        else {
          discTwoPlot(xx, yy, multPlot = multPlot,
                      xypred.disc, data, th[i], xlab = xlab, 
                      ylab = ylab, pch = pch, pcol = pcol, cex = cex, 
                      cex.lab = cex.lab, cex.axis = cex.axis, yRange = yRange, 
                      ModTitle = ModTitle[i], TiAdj = TiAdj, TiLine = TiLine, 
                      cex.main = cex.main, lwd = lwd, lcol = lcol, 
                      ...)
        }
      }
      else {
        contsPlot(xx, yy, multPlot = multPlot, 
                  xypred.cont, data, column = names[[i]], 
                  xlab = xlab, ylab = ylab, pch = pch, pcol = pcol, 
                  cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, 
                  yRange = yRange, ModTitle = ModTitle[i], TiAdj = TiAdj, 
                  TiLine = TiLine, cex.main = cex.main, lwd = lwd, 
                  lcol = lcol)
      }
    }
    par(mfrow = c(1, 1))#change par back to default
    } else { #plot all model fits on the same plot
      
      plot(xx, yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol, 
           cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange,
           ...)
      if (!is.null(ModTitle) & length(ModTitle) == 1){
        title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main, 
            ...)
      }

      lcol2 <- lcol

      for (i in 1:length(mods)) {
        if (length(lcol2) > 1) lcol <- lcol2[i] #if lcol is a vector of colours
        if (any(c("DiscOne", "DiscTwo") %in% names[[i]])) {
          if (names[[i]] == "DiscOne") {
            discOnePlot(xx, yy, multPlot = multPlot, xypred.disc, data, 
                        th[i], xlab = xlab, 
                        ylab = ylab, pch = pch, pcol = pcol, cex = cex, 
                        cex.lab = cex.lab, cex.axis = cex.axis, yRange = yRange, 
                        ModTitle = ModTitle[i], TiAdj = TiAdj, TiLine = TiLine, 
                        cex.main = cex.main, lwd = lwd, lcol = lcol, 
                        ...)
          }
          else {
            discTwoPlot(xx, yy, multPlot = multPlot, xypred.disc, 
                        data, th[i], xlab = xlab, 
                        ylab = ylab, pch = pch, pcol = pcol, cex = cex, 
                        cex.lab = cex.lab, cex.axis = cex.axis, yRange = yRange, 
                        ModTitle = ModTitle[i], TiAdj = TiAdj, TiLine = TiLine, 
                        cex.main = cex.main, lwd = lwd, lcol = lcol, 
                        ...)
          }
        }
        else {
          contsPlot(xx, yy, multPlot = multPlot, 
                    xypred.cont, data, column = names[[i]], 
                    xlab = xlab, ylab = ylab, pch = pch, pcol = pcol, 
                    cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, 
                    yRange = yRange, ModTitle = ModTitle[i], TiAdj = TiAdj, 
                    TiLine = TiLine, cex.main = cex.main, lwd = lwd, 
                    lcol = lcol, 
                    ...)
        }#eo if DiscOne / DiscTwon
      }# eo for i
    }# eo if multPlot
  }#eo if one plot
}#eo function




#'Barplot of information criteria weights for a 'habitat' Object
#'
#'@description S3 method for class 'habitat'. \code{plot.habitat} creates
#'  a simple barplot of information criteria weights for the different model fits
#'@param x An object of class 'habitat'.
#'@param IC The information criterion weights to present (must be one of 'AIC',
#''BIC' or 'AICc')
#'@param \dots Further graphical parameters may be supplied as arguments.
#' @examples
#' #data(habitat)
#' ###TO FILL IN
#'@importFrom graphics barplot
#'@export

plot.habitat <- function(x,  
                        IC = "AICc",
                        type = 1,
                        powFit = TRUE,
                        lcol = NULL,
                        pLeg = TRUE,
                        legPos = "bottomright",
                        ...){
  
  if (attributes(x)$type == "habitat"){
    if (!IC %in% c("AIC", "BIC", "AICc")){
      stop("IC must be one of 'AIC', 'BIC', 'AICc'")
    }
    x2 <- summary(x)$Model_table
    #get delta ICs
    delta_ICs <- x2[,IC] - min(x2[,IC])
    #get akaike weights
    akaikesum <- sum(exp(-0.5*(delta_ICs)))
    x2$Weights <- exp(-0.5*delta_ICs) / akaikesum
    IC_nam <- paste0(IC, " weight")
    barplot(x2$Weights, names.arg = x2$Model,
            xlab = "Model", ylab = IC_nam, ...)
  } else if (attributes(x)$type == "countryside"){#eo if habitat
    
    if (length(x[[4]]) == 1){
      return("Plot not generated as some models could not be fitted")
    }
    
    dd <- x[[6]]
    
    dd_Area <- length(which(grepl("Area", colnames(dd))))
    dd_SR <- ncol(dd) - dd_Area
    
    dd2_Area <- dd[,1:dd_Area]
    dd2_SR <- dd[,(dd_Area + 1):ncol(dd)]
    
    dd_Ran <- range(dd[,1:dd_Area])
    
    if (type == 1){
    
      ##total predicted richness for the actual
      #add total area and totR columns in
      dd3_Area <- as.data.frame(dd2_Area)
      dd3_Area$totA <- rowSums(dd3_Area)
      dd3_Area$totR <-  f[[4]]
      #same for main dataframe
      ddTot <- dd
      ddTot$totA <- rowSums(dd2_Area)
      ddTot$totR <- rowSums(dd2_SR)
      
      if(!identical(ddTot$totA, dd3_Area$totA)){
        stop("rownames mismatch in plot.countryside,",
             " contact the package author")
      }
      
      plot(ddTot$totR, dd3_Area$totR,
           xlab = "Observed total richness",
           ylab = "Predicted total richness",
           ...)
      abline(0,1)
      
      #extract power model values
      if (powFit){
      if (length(x[[7]]) > 1){
        points(x[[7]]$data$S, x[[7]]$calculated,
               col = "red")
      } else {
        cat("\n\nPower model could not be fitted\n\n")
      }#eo if f7
      }#ep if powFit

    } else if (type == 2) {
    #check if ubiquitous sp included
    UB <- ifelse(ncol(x[[6]])%%2==0, FALSE, TRUE)
      
    ##predicted curves for each land-use
    Ar_seq <- seq(dd_Ran[1], dd_Ran[2],
                  length.out = 1000)
    #convert in N tables, where in each you can N columns,
    #where N = number of land-use types. In each all columns,
    #except the focal habitat are zeros
    ar_ls <- vector("list", length = dd_Area)
    totR_i <- matrix(ncol = dd_Area, nrow = length(Ar_seq))
    UB_i <- matrix(ncol = dd_Area, nrow = length(Ar_seq))
    colnames(totR_i) <- names(f[[2]][[1]])
    for (i in 1:dd_Area){
      m_ls <- matrix(0, ncol = dd_Area,
                     nrow = length(Ar_seq))
      m_ls[,i] <- Ar_seq
      totR_R <- apply(m_ls[,1:dd_Area],1,function(y){
        v <- as.vector(y)
        vc <- countryside_extrap(x, area = v)
        vc1 <- vc$Indiv_mods[i]
        if (UB){
          vc2 <- vc$Indiv_mods[length(vc$Indiv_mods)]
        } else{
          vc2 <- 0
        }
        c(vc1, vc2)
      })
      totR_i[,i] <- totR_R[1,]
      UB_i[,i] <- totR_R[2,]#ubiquitous species vals
    }#eo for i
    
    totR_i <- as.data.frame(totR_i)
    #take row means for ubiquitous species
    UB_i2 <- data.frame("Ub" = rowMeans(UB_i))
    
    COLs <- c("black", "red", "darkgreen", "blueviolet",
              "brown", "cornflowerblue","darkorange4",
              "deeppink4", "gold4", "gray16", "lightgreen")
    
    if (length(lcol) > 1){
      spC <- ifelse(UB == TRUE, ncol(totR_i) + 1, ncol(totR_i))
      if (length(lcol) != spC){
        warning("Length of lcol does not match number of habitat",
                " types - using randomly selected colours")
        lcol <- COLs[1:spC]
      }
    } else {
      lcol <- COLs[1:spC]
    }
    
    plot(Ar_seq, totR_i[,1], type = "l", 
         col = lcol[1],
         xlab = "Area", ylab = "Total richness",
         ylim = c(min(totR_i), max(totR_i)),
         ...)
    k <- 1
    apply(totR_i[,2:ncol(totR_i), drop = FALSE],
          2, function(y){
            k <<- k + 1
            lines(Ar_seq, y, col = lcol[k], ...)
          })
    k <- k + 1
    lines(Ar_seq, UB_i2[,1], col = lcol[k], ...)
    if (pLeg){
      CNz <- colnames(totR_i)
      if (UB) {CNz <- c(CNz, "UB")}
      legend(legPos, 
             legend = CNz,
             col = lcol, lty = 1)
    }
  } else {
      stop("Type should be either 1 or 2")
  }#eo if type 1
}#eo if countryside
}


#function to convert vector of model
#names into abbreviated versions depending on which models are provided
mod_abbrev <- function(nams){

x1 <-  c("power", "powerR","epm1","epm2","p1","p2","loga","koba","mmf",
         "monod","negexpo","chapman",
         "weibull3","asymp","ratio","gompertz","weibull4","betap","logistic",
         "heleg","linear")

x2 <- c("Pow", "PowR", "E1", "E2", "P1", "P2", "Loga", "Kob", "MMF",
        "Mon", "NegE",
        "CR", "CW3", "AR", "RF", "Gom", "CW4", "BP", "Logi", "Hel", "Lin")

df <- data.frame("Full_name" = x1, "Abbreviated_name" = x2)
dfb <- vapply(nams, function(x) which(df$Full_name == x), FUN.VALUE = numeric(1))
df2 <- df[dfb,]
if (nrow(df2) != length(nams)) stop("Not enough matched model names")
as.vector(df2$Abbreviated_name)
}
