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
#'   \code{\link[graphics]{plot}},\code{\link[graphics]{title}},
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
#' fc <- sar_multi(data = galap, obj = c("power", "loga", "epm1"))
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
    return(cat("\nNo plot method for 'pred' object of class 'sars'\n", sep = ""))
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
    ff <- x$calculated
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
    lines(x = xx, y = ff, lwd = lwd, col = lcol,  ...)

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
      lapply(x, function(x){

      df <- x$data
      xx <- df$A
      yy <- df$S
      ff <- x$calculated
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
      matlines(xx, mf2, lwd = lwd, lty = seq_along(mf2), col=seq_along(mf2))
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
#'   relative to the default.
#' @param yRange The range of the y-axis. Only for use with \code{type =
#'   multi}.
#' @param lwd Line width. Only for use with \code{type = multi}.
#' @param lcol Line colour. Only for use with \code{type = multi}.
#' @param pLeg Logical argument specifying whether or not the legend should
#'   be plotted  (when \code{type = multi} and \code{allCurves = TRUE}).
#' @param modNames A vector of model names for the barplot of weights (when
#'   \code{type = bar}). The default (\code{modNames = NULL}) uses
#'   abbreviated versions (see below) of the names from the \code{sar_average}
#'   function.
#' @param cex.names The amount by which the axis labels (model names) should
#'   be scaled relative to the default. Only for use with \code{type = bar}.
#' @param subset_weights Only create a barplot of the model weights for
#'   models with a weight value above a given threshold
#'   (\code{subset_weights}). Only for use with \code{type = bar}.
#' @param confInt A logical argument specifying whether confidence intervals
#'   should be plotted around the multimodel curve. Can only be used if
#'   confidence intervals have been generated in the \code{sar_average}
#'   function.
#' @param \dots Further graphical parameters (see
#'   \code{\link[graphics]{par}},
#'   \code{\link[graphics]{plot}},\code{\link[graphics]{title}},
#'   \code{\link[graphics]{lines}}) may be supplied as arguments.
#' @note When plotting all model fits on the same plot with a legend it is
#'   necessary to manually extend your plotting window (height and width;
#'   e.g. the 'Plots' window of R studio) before plotting to ensure the
#'   legend fits in the plot. Extending the plotting window after plotting
#'   simply stretches the legend.
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
#'   abbreviations used as the default are: \itemize{ \item{Pow = } { Power}
#'   \item{PowR = } { PowerR} \item{E1 = } { Extended_Power_model_1} \item{E2
#'   = } { Extended_Power_model_2} \item{P1 = } { Persistence_function_1}
#'   \item{P2 = } { Persistence_function_2} \item{Loga = } { Logarithmic}
#'   \item{Kob = } { Kobayashi} \item{MMF = } { MMF} \item{Mon = } { Monod}
#'   \item{NegE = } { Negative_exponential} \item{CR = } { Chapman_Richards}
#'   \item{CW3 = } { Cumulative_Weibull_3_par.} \item{AR = } {
#'   Asymptotic_regression} \item{RF = } { Rational_function} \item{Gom = } {
#'   Gompertz} \item{CW4 = } { Cumulative_Weibull_4_par.} \item{BP = } {
#'   Beta-P_cumulative} \item{Hel = } { Heleg(Logistic)} \item{Lin = } {
#'   Linear_model}}
#' @examples
#' data(galap)
#' #plot a multimodel SAR curve with all model fits included
#' fit <- sar_average(data = galap)
#' plot(fit)
#'
#' #remove the legend
#' plot(fit, pLeg = FALSE)
#'
#' #plot just the multimodel curve
#' plot(fit, allCurves = FALSE, ModTitle = "", lcol = "black")
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
                      lwd = 2, lcol = 'dodgerblue2', pLeg = TRUE,
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

  #fitted values for each model
  mf <- lapply(dat2, function(x) x$calculated)
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
  if (!all(round(wfv) == round(x$mmi)))
    stop("Multimodel fitted values do not match between functions")

  if (allCurves){
    mf2$MultiModel <- wfv
    nams2 <- c(nams, "Multimodel SAR")
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
                xlim = c(min(xx), xMAX), ylim = yRange)
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
      matlines(xx, mf2, lwd = lwd, lty = seq_along(mf2), col=seq_along(mf2))
      if (pLeg) legend(max(xx) + (max(xx) * 0.05), yMax,
                               legend = nams2, horiz = FALSE,
                               lty = seq_along(mf2), col=seq_along(mf2))
      title(main = ModTitle, adj = TiAdj, line = TiLine,cex.main = cex.main)
  } else if (!allCurves){
    #just multimodel SAR curve
    if (confInt){
      plot(x = xx, y = yy, xlab = xlab, ylab = ylab,
            cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange)
      polygon(c(xx,rev(xx)),c(CI$L,rev(CI$U)),col="grey87",border=NA)
      points(x = xx, y = yy, pch = pch, col = pcol,
             cex = cex)
      title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main)
      lines(x = xx, y = wfv, lwd = lwd, col = lcol)
    } else {
      plot(x = xx, y = yy, xlab = xlab, ylab = ylab, pch = pch, col = pcol,
      cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, ylim = yRange)
      title(main = ModTitle, adj = TiAdj, line = TiLine, cex.main = cex.main)
      lines(x = xx, y = wfv, lwd = lwd, col = lcol)
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


#function to convert vector of model
#names into abbreviated versions depending on which models are provided
mod_abbrev <- function(nams){

x1 <-  c("power", "powerR","epm1","epm2","p1","p2","loga","koba","mmf",
         "monod","negexpo","chapman",
         "weibull3","asymp","ratio","gompertz","weibull4","betap","heleg",
         "linear")

x2 <- c("Pow", "PowR", "E1", "E2", "P1", "P2", "Loga", "Kob", "MMF",
        "Mon", "NegE",
        "CR", "CW3", "AR", "RF", "Gom", "CW4", "BP", "Hel", "Lin")

df <- data.frame("Full_name" = x1, "Abbreviated_name" = x2)
df2 <- df[(which(df$Full_name %in% nams)),]
if (nrow(df2) != length(nams)) stop("Not enough matched model names")
as.vector(df2$Abbreviated_name)
}
