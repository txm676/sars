

########################################################################
##     internal functions to fit threshold models, using lm      #######
########################################################################

#' function to find the threshold for one threshold cont and zero slope
#' @noRd
find_one_threshold_cont <- function(x, y, fct, interval, nisl = NULL){
  if (is.null(nisl)){
    sequence <- seq(min(x), max(x), interval)
    } else {
    n <- nisl-1
    if (n == 0) {n <- 1}
    sequence <- seq(min(sort(x)[-(1:n)]), max(rev(sort(x))[-(1:n)]), interval)
  }
  s1 <- lapply(sequence, fct, x = x, y = y)
  w <- which.min(s1)
  if (length(w) == 1){
    threshold <- sequence[w]
  } else{
    w2 <- w[sample(1:length(w), 1)]
    threshold <- sequence[w2]
  }
  return(threshold)
}

#' function to find the threshold for one threshold discon
#' @noRd
find_one_threshold_disc <- function(x, y, fct, nisl = NULL){
  if (is.null(nisl)){
    x1 <- x
    } else {
    n <- nisl-1
    if (n == 0) {n <- 1}
    x1 <- x[x >= min(sort(x)[-(1:n)]) & x <= max(rev(sort(x))[-(1:n)])]
  }
  rss <- vector(length = length(x1))
  for (i in 1:length(x1)){
    rss[i] <- fct(x[i], x, y)
  }
  w <- which.min(rss)
  if (length(w) == 1){
    threshold <- x1[w]
  } else{
    w2 <- w[sample(1:length(w), 1)]
    threshold <- x1[w2]
  }
  return(threshold)
}

#' function to find the threshold for two threshold cont and zero
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach "%dopar%"
#' @noRd
find_two_thresholds_cont <- function(x, y, fct, interval, nisl = NULL, parallel, 
                                     cores){
  if (is.null(nisl)){
    sequence <- seq(min(x), max(x), interval)
    } else {
    n <- nisl-1
    if (n == 0) {n <- 1}
    sequence <- seq(min(sort(x)[-(1:n)]), max(rev(sort(x))[-(1:n)]), interval)
  }
  N <- length(sequence)
  if (parallel){
    cores <- cores
    cl <- makeCluster(cores); on.exit(stopCluster(cl))
    registerDoParallel(cl)
    i <- 1 #Dummy line 
    ssr_t1 <- foreach(i=seq(from= 1, to=N-1, by=1))  %dopar% { 
      ssr_t2 <- vector("list", length = length((i+1):N))
      k <- 1
      for (j in (i + 1):N){
        ssr_t2[[k]] <- c(fct(sequence[i], sequence[j], x, y), 
                         sequence[i], sequence[j])
      }
      k <- k + 1
      ssr_t2
    }#eo dopar
  } else{
    ssr_t1 <- vector("list", length = N - 1)
    for (i in 1:(N - 1)){
      ssr_t2 <- vector("list", length = length((i + 1):N))
      j <- 1
      for (j in (i + 1):N){
        ssr_t2[[j]] <- c(fct(sequence[i], sequence[j], x, y), sequence[i], 
                         sequence[j])
      }
      j <- j + 1
      ssr_t1[[i]] <- ssr_t2
    }
  }
  l2 <- do.call(rbind, lapply(ssr_t1, function(x) do.call(rbind, x)))
  thb <- l2[which(l2[,1] == min(l2[,1])), , drop = FALSE]
  
  th <- l2[which(l2[,1] == min(l2[,1])),][2:3]
  
  if (nrow(thb) == 1){
    th <- as.vector(thb)[2:3]
  } else {
    rr <- sample(1:nrow(th), 1)
    th <- as.vector(thb[rr,])[2:3]
  }
  return(th)
}


#' function to find the threshold for two threshold disc
#' @noRd
find_two_thresholds_disc <- function(x, y, fct, nisl = NULL){
  if(is.null(nisl)){n <- 0}
  if(is.null(nisl)){
    x1 <- x
    } else {
    n <- nisl - 1
    if (n == 0) {n <- 1}
    x1 <- x[x >= min(sort(x)[-(1:n)]) & x <= max(rev(sort(x))[-(1:n)])]
  }
  N <- length(x1) - 1
  ssr_t1 <-  vector("list", length = N)
  for(i in 1:(N)){
    ssr_t2 <- vector("list", length = length((i + 1):(N)))
    j <- 1
    for (j in (i + 1):(N)){
      ssr_t2[[j]] <- c(fct(x1[i], x1[j], x, y), x1[i], x1[j])
    }
    j <- j + 1
    ssr_t1[[i]] <- ssr_t2
  }
  l2 <- do.call(rbind, lapply(ssr_t1, function(x) do.call(rbind, x)))
  thb <- l2[which(l2[,1] == min(l2[,1])), , drop = FALSE]
  if (nrow(thb) == 1){
    th <- as.vector(thb)[2:3]
  } else {
    rr <- sample(1:nrow(thb), 1)
    th <- as.vector(thb[rr,])[2:3]
  }
  return(th)
}


#rss functions for each model

#' one thr continuous rss
#' @importFrom stats predict
#' @noRd
fct_cont_one <- function(th, x, y) {
  sum((y-stats::predict(lm(y ~ x + I((x - th)*as.numeric(x > th)))))^2)
}

#' zero slope rss
#' @importFrom stats predict
#' @noRd
fct_zslope_one <- function(th, x, y) {
  sum((y-stats::predict(lm(y~1+I((x-th)*as.numeric(x > th)))))^2)
}

#' one thr discontinuous rss
#' @importFrom stats predict
#' @noRd
fct_disc_one <- function(th, x, y) {
  sum((y-predict(lm(y ~ x*(x<=th) + x*(x>th))))^2)
} # modified

#' two thr continuous rss
#' @importFrom stats predict
#' @noRd
fct_cont_two <- function(th, th2, x, y) {
  sum((y-stats::predict(lm(y ~ x + I((x - th)*as.numeric(x > th)) + 
                      I((x - th2)*as.numeric(x > th2)))))^2)
}

#' zero slope 2 thr rss
#' @importFrom stats predict
#' @noRd
fct_zslope_two <- function(th, th2, x, y) {
  sum((y-stats::predict(lm(y ~ 1 + I((x - th)*as.numeric(x > th)) + 
                      I((x - th2)*as.numeric(x > th2)))))^2)
}

#' two thr discontinuous rss
#' @importFrom stats predict
#' @noRd
fct_disc_two <- function(th, th2, x, y) {
  sum((y-predict(lm(y ~ x*(x<=th) + x*(x>th & x<=th2) + x*(x>th2))))^2)
} # modified


#' Fit threshold SAR models
#'
#' @description Fit up to six piecewise (threshold) regression models to SAR
#'   data.
#' @usage sar_threshold(data, mod = "All", interval = NULL, nisl = NULL,
#'   non_th_models = TRUE, logAxes = "none", con = 1, logT = log, parallel =
#'   FALSE, cores = NULL)
#' @param data A dataset in the form of a dataframe with two columns: the first
#'   with island/site areas, and the second with the species richness of each
#'   island/site.
#' @param mod A vector of model names: an individual model, a set of models, or
#'   all models. Can be any of 'All' (fit all models), 'ContOne' (continuous
#'   one-threshold), 'ZslopeOne' (left-horizontal one-threshold), 'DiscOne'
#'   (discontinuous one-threshold), 'ContTwo' (continuous two-threshold),
#'   'ZslopeTwo' (left-horizontal two-threshold), or 'DiscTwo' (discontinuous
#'   two-threshold).
#' @param interval The amount to increment the threshold value by in the
#'   iterative model fitting process (not applicable for the discontinuous
#'   models). The default for non-transformed area reverts to 1, while for
#'   log-transformed area it is 0.01. However, these values may not be suitable
#'   depending on the range of area values in a dataset, and thus users are
#'   advised to manually set this argument.
#' @param nisl ****(n.b. it needs to be > xxx)
#' @param non_th_models Logical argument (default = TRUE) of whether two
#'   non-threshold models (i.e. a simple linear regression: y ~ x; and an
#'   intercept only model: y ~ 1) should also be fitted.
#' @param logAxes What log-transformation (if any) should be applied to the area
#'   and richness values. Should be one of "none" (no transformation), "area"
#'   (only area is log-transformed) or "both" (both area and richness
#'   log-transformed).
#' @param con The constant to add to the species richness values in cases where
#'   one of the islands has zero species.
#' @param logT The log-transformation to apply to the area and richness values.
#'   Can be any of \code{log}(default), \code{log2} or \code{log10}.
#' @param parallel Logical argument for whether parallel processing should be
#'   used. Only applicable when the continuous two-threshold and left-horizontal
#'   two-threshold models are being fitted.
#' @param cores Number of cores to use. Only applicable when \code{parallel =
#'   TRUE}.
#' @details Fitting the continuous and left-horizontal piecewise models
#'   (particularly the two-threshold models) can be time consuming if the
#'   range in area is large and/or the \code{interval} argument is small. For
#'   the two-threshold continuous slope and left-horizontal models, the use of
#'   parallel processing (using the \code{parallel} argument) is recommended.
#'   The number of cores (\code{cores}) must be provided.
#'
#'   Note that the interval argument is not used to fit discontinuous models,
#'   as, in these cases, the breakpoint must be at a datapoint.
#'
#'   There has been considerable debate regarding the number of parameters that
#'   are included in different piecewise models. Here (and thus in our
#'   calculation of AIC, AICc, BIC etc) we consider ContOne to have five
#'   parameters, ZslopeOne - 4, DiscOne - 6, ContTwo - 7, ZslopeTwo - 6, DiscTwo
#'   - 8. The standard linear model and the intercept model are considered to
#'   have 3 and 2 parameters, respectively. The raw \code{\link{lm}} model fits
#'   are provided in the output, however, if users want to calculate information
#'   criteria using different numbers of parameters.
#'   
#'   The raw \code{\link{lm}} model fits can also be used to explore classic
#'   diagnostic plots for linear regression analysis in R using the function
#'   \code{\link{plot}} or other diagnostic tests such \code{outlierTest},
#'   \code{leveragePlots} or \code{influencePlot}, available in the \code{car}
#'   package. This is advised as currently there are no model validation checks
#'   undertaken automatically, unlike elsewhere in the sars package. 
#'   
#'   Confidence intervals around the breakpoints in the one-threshold continuous
#'   and left- horizontal models can be calculated using the
#'   \code{\link{threshold_ci}} function. The intercepts and slopes of the
#'   different segments in the fitted breakpoint models can be calculated using
#'   the \code{\link{get_coef}} function.
#'   
#'   Rarely, multiple breakpoint values can return the same minimum rss (for a
#'   given model fit). In these cases, we just randomly choose and return one
#'   and also produce a warning. If this occurs it is worth checking the data
#'   and model fits carefully.
#'   
#' @return A list of class "threshold" and "sars" with five elements. The first
#'   element contains the different model fits (lm objects). The second element
#'   contains the names of the fitted models, the third  contains the threshold
#'   values, the fourth element the dataset (i.e. a dataframe with area and
#'   richness values), and the fifth contains details of any axes
#'   log-transformations undertaken. \code{\link{summary.sars}} provides a more
#'   user-friendly ouput (including a model summary table) and
#'   \code{\link{plot.threshold}} plots the model fits.
#' @note Due to the increased number of parameters, fitting piecewise regression
#'   models to datasets with few islands is not recommended. In particular, we
#'   would advise against fitting the two-threshold models to small SAR datasets
#'   (e.g. fewer than 10 islands for the one threshold models, and 20 islands
#'   for the two threshold models).
#' @references Lomolino, M.V. & Weiser, M.D. (2001) Towards a more general
#'   species-area relationship: diversity on all islands, great and small.
#'   Journal of Biogeography, 28, 431-445.
#'
#'   Gao, D., Cao, Z., Xu, P. & Perry, G. (2019) On piecewise models and
#'   species-area patterns. Ecology and Evolution, 9, 8351-8361.
#'
#'   Matthews, T.J. et al. (In review) Unravelling the small-island effect
#'   through phylogenetic community ecology
#' @author Francois Rigal and Thomas J. Matthews
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

sar_threshold <- function(data, mod = "All", interval = NULL, nisl = NULL,
                          non_th_models = TRUE, logAxes = "none", 
                          con = 1, logT = log,
                          parallel = FALSE, cores = NULL){
  
  if (!(is.matrix(data) | is.data.frame(data)))
    stop('data must be a matrix or dataframe')
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop('NAs present in data')
  if (!(is.character(mod) | is.vector(mod))) 
    stop("mod should be character vector")
  if(!is.null(interval)){
    if (!is.numeric(interval) | length(interval) != 1)
      stop("interval should be a numeric vector of length 1")
    if (interval > max(data[ ,1])) stop("interval must be smaller than max area")
  }
  if ("All" %in% mod & length(mod) > 1) stop("length(mod) > 1 and contains All")
  if (!all(mod %in% c("All", "ContOne", "ZslopeOne", 
                      "DiscOne", "ContTwo", "ZslopeTwo", "DiscTwo")))
    stop ("Incorrect model names provided; see help for 'mod' argument")
  if (parallel) {
    if (is.null(cores)) 
      stop("cores argument must be provided for parallel processing")
    if (!is.numeric(cores) | length(cores) > 1) 
      stop("cores should a numeric vector of length 1")
  }
  if (!any(c("none", "area", "both") %in% logAxes)){
    stop("logAxes should be one of 'none', 'area', or 'both'")
  }
  if (!is.primitive(logT)) stop("logT should be a (primitive) function,
                                specifically: log, log2 or log10")
  if (any(length(con) > 1 | !(is.numeric(con)))) 
    stop("con should be a numeric vector of length 1")
  if (nrow(data) < 10) warning("Sample size is quite small for fitting",
                               " threshold models, see documentation")
  if ((!is.null(nisl)) && nisl > length(data[,1]) / 2){
    stop('nisl is larger than half of the total number of islands')
  }
  
  #create list to hold final model fits;
  #and create a names vector to provide, as naming the list elements 
  #seems to delete the lm content
  #and one for optimum threshold values
  if ("All" %in% mod & non_th_models){
    res <- vector("list", length = 8)
    names <- vector(length = 8)
    th <- vector("list", length = 8)
  } else if ("All" %in% mod & (!non_th_models)){
    res <- vector("list", length = 6)
    names <- vector(length = 6)
    th <-vector("list", length = 6)
  } else if (non_th_models){
    res <- vector("list",length = length(mod) + 2)
    names <- vector(length = length(mod) + 2)
    th <-vector("list", length = length(mod) + 2)
  } else{
    res <- vector("list",length = length(mod))
    names <- vector(length = length(mod))
    th <-vector("list", length = length(mod))
  }
  r <- 1 #counter
  
  data <- data[order(data[,1]),]
  colnames(data) <- c('A','S')
  
  #log conversion (if needed)
  if (logAxes == "area"){
    data$A <- logT(data$A)
    if (is.null(interval)) interval <- 0.01
  } else if (logAxes == "both"){
    data$A <- logT(data$A)
    if (is.null(interval)) interval <- 0.01
    if (any(data$S == 0)){
      data$S <- logT(data$S + con)
    } else{
      data$S <- logT(data$S)
    }
  }
  
  if (is.null(interval)) interval <- 1

  y <- data$S
  x <- data$A
  
  ###get the thresholds############
  
  if (any(c("All", "ContOne") %in% mod)) {
    #one breakpoint continuous
    threshold_cont <- find_one_threshold_cont(x, y, fct_cont_one, 
                                              interval, nisl)
    res[[r]] <- lm(y ~ x + I((x - threshold_cont) * 
                               as.numeric(x > threshold_cont)))
   # names(res[[r]]$coefficients) <- c("Intercept", "z1", 
    #                                  "z2")
    names[r] <- "ContOne"
    th[[r]] <- threshold_cont
    r <- r + 1
  }
  #one breakpoint zero slope
  if (any(c("All", "ZslopeOne") %in% mod)) {
    threshold_zslope <- find_one_threshold_cont(x, y, fct_zslope_one, 
                                                interval, nisl)
    res[[r]] <- lm(y ~ 1 + I((x - threshold_zslope) * 
                               as.numeric(x > threshold_zslope)))
   # names(res[[r]]$coefficients) <- c("Intercept", "z2")
    names[r] <- "ZslopeOne"
    th[[r]] <- threshold_zslope
    r <- r + 1
  }
  
  #one breakpoint discontinuous
  if (any(c("All", "DiscOne") %in% mod)) {
    threshold_disc <- find_one_threshold_disc(x, y, fct_disc_one, nisl)
    
    res[[r]] <- lm(y ~ x*(x <= threshold_disc[1]) + x*(x > threshold_disc[1]))
    
    names[r] <- "DiscOne"
    th[[r]] <- threshold_disc
    r <- r + 1
  }
  #two threshold continuous
  if (any(c("All", "ContTwo") %in% mod)) {
    threshold_two_cont <- find_two_thresholds_cont(x, y, 
                                                   fct_cont_two, interval,
                                                   nisl,
                                                   parallel = parallel, 
                                                   cores = cores)
    res[[r]] <- lm(y ~ x + I((x - threshold_two_cont[1]) * 
                               as.numeric(x > threshold_two_cont[1])) + 
                     I((x - threshold_two_cont[2]) * 
                         as.numeric(x > threshold_two_cont[2])))
    names[r] <- "ContTwo"
    th[[r]] <- threshold_two_cont
    r <- r + 1
  }
  #two threshold zero slope
  if (any(c("All", "ZslopeTwo") %in% mod)) {
    threshold_two_zslope <- find_two_thresholds_cont(x, y, 
                                                     fct_zslope_two, 
                                                     interval, nisl,
                                                     parallel = parallel, 
                                                     cores = cores)
    res[[r]] <- lm(y ~ 1 + I((x - threshold_two_zslope[1]) * 
                               as.numeric(x > threshold_two_zslope[1])) + 
                     I((x - threshold_two_zslope[2]) * 
                         as.numeric(x > threshold_two_zslope[2])))
    names[r] <- "ZslopeTwo"
    th[[r]] <- threshold_two_zslope
    r <- r + 1
  }
  
  #two breakpoint discontinuous
  if (any(c("All", "DiscTwo") %in% mod)) {
    threshold_two_disc <- find_two_thresholds_disc(x, y, 
                                                   fct_disc_two, nisl)
    res[[r]] <- lm(y ~ x * (x <= threshold_two_disc[1]) + 
                     x*(x > threshold_two_disc[1] & x<=threshold_two_disc[2]) + 
                     x * (x > threshold_two_disc[2]))
    
    names[r] <- "DiscTwo"
    th[[r]] <- threshold_two_disc
    r <- r + 1
  }

  ##fit the final non-threshold models for comparison
  if (non_th_models){
    res[[r]] <- model_lm <- lm(y ~ x)
    th[[r]] <- NA
    r <- r + 1
    res[[r]]<- model_null <- lm(y ~ 1)
    th[[r]] <- NA
    names[(r-1):r] <- c("Linear", "Intercept")
  }
  if (logAxes == "none") logT = "none"
  res2 <- list(res, names, th, data, c(logAxes, logT))
  class(res2) <- c("threshold", "sars", "list")
  attr(res2, "type") <- "threshold"
  return(res2)
}


######################################################################
#' Calculate confidence intervals around breakpoints
#'
#' @description Generate confidence intervals around the breakpoints of the
#' one-threshold continuous and left-horizontal models, using bootstrapping.
#' @usage threshold_ci(object, method = "boot", interval = NULL, Nboot = 100,
#'   verb = TRUE)
#' @param object An object of class 'thresholds', generated using the
#'   \code{\link{sar_threshold}} function. The object must contain fits of
#'   either (or both) of the one-threshold continuous or the one-threshold
#'   left-horizontal model.
#' @param method *** should be one of \code{boot} or \code{F}.
#' @param interval The amount to increment the threshold value by in the
#'   iterative model fitting process. The default for non-transformed area
#'   reverts to 1, while for log-transformed area it is 0.01. It is advised that
#'   the same interval value used when running \code{\link{sar_threshold}} is
#'   used here.
#' @param Nboot Number of bootstrap samples (for use with \code{method =
#'   "boot"}).
#' @param verb Should progress be reported. If \code{TRUE}, every 50th bootstrap
#'   sample is reported (for use with \code{method = "boot"}).
#' @details *F* method
#'
#'   Calculates confidence intervals (95 percent) around the breakpoint
#'   estimates using the bootstrap approach of Toms & Lesperance (2003). If the
#'   number of bootstrap samples (\code{Nboot}) is large, the function can take
#'   a while to run.
#'
#'   Currently only available for the one-threshold continuous and left-
#'   horizontal threshold models.
#'
#' @return A list of class "sars" with two elements. The first element contains
#'   the full set of bootstrapped breakpoint estimates for each model. The
#'   second contains the 95 percent confidence interval values.
#' @references Toms, J.D. & Lesperance, M.L. (2003) Piecewise regression: a tool
#'   for identifying ecological thresholds. Ecology, 84, 2034-2041.
#' @examples
#' data(aegean2)
#' a2 <- aegean2[1:168,]
#' fitT <- sar_threshold(data = a2, mod = "ContOne", 
#' interval = 0.1, non_th_models = TRUE, logAxes = "area", logT = log10)
#' #calculate confidence intervals (using very low Nboot just as an example)
#' CI <- threshold_ci(fitT, method = "boot", interval = NULL, Nboot = 3)
#' CI
#' ##F EXAMPLE
#' @importFrom stats fitted resid
#' @export

threshold_ci <- function(object, method = "boot", interval = NULL, 
                         Nboot = 100, verb = TRUE) 
{
  if (!"threshold" %in% class(object)) 
    stop("Object should be of class 'threshold'")
  names <- object[[2]]
  if (!any(names %in% c("ContOne", "ZslopeOne"))) 
    stop("Confidence interval is only available for models ContOne and ZslopeOne")
  if (!is.null(interval)) {
    if (!is.numeric(interval) | length(interval) != 1) 
      stop("interval should be a numeric vector of length 1")
    if (interval > max(object[[4]]$A)) 
      stop("interval must be smaller than max", "area")
  }
  if (is.null(interval)) {
    logAxes <- object[[5]][1]
    if (logAxes == "none") {
      interval <- 1
    }
    else {
      interval <- 0.01
    }
  }
  if (method == "boot") {
    w1 <- which(names %in% c("ContOne", "ZslopeOne"))
    bootR <- vector("list", length = 3)
    bootR[[1]] <- vector("list", length = length(w1))
    names(bootR[[1]]) <- names[w1]
    bootR[[2]] <- vector("list", length = length(w1))
    names(bootR[[2]]) <- names[w1]
    k <- 1
    names(bootR) <- c("Bootstrap values", "CIs")
    for (j in w1) {
      n1 <- names[[j]]
      if (n1 == "ContOne") {
        fct <- fct_cont_one
      }
      else {
        fct <- fct_zslope_one
      }
      mods <- object[[1]][[j]]
      x <- object[[4]]$A
      y <- object[[4]]$S
      res <- resid(mods)
      fit <- fitted(mods)
      new.df <- data.frame(res, x)
      boot <- vector(length = Nboot)
      for (i in 1:Nboot) {
        new.res <- sample(new.df$res, replace = T)
        xbt <- new.df$x
        ybt <- fit + new.res
        sequence <- seq(min(xbt), max(xbt), interval)
        s1 <- lapply(sequence, fct, x = xbt, y = ybt)
        w <- which.min(s1)
        if (length(w) == 1) {
          boot[i] <- sequence[w]
        }
        else {
          w2 <- w[sample(1:length(w), 1)]
          boot[i] <- sequence[w2]
        }
        if (verb) {
          if (i%%50 == 0) 
            cat("Bootstrap sample", i, "out of", Nboot, 
                "for model", k, "of", length(w1), "\n")
        }
      }
      CI <- quantile(boot, c(0.025, 0.975))
      bootR[[1]][[k]] <- boot
      bootR[[2]][[k]] <- CI
      k <- k + 1
    }
    bootR[[3]] <- "boot"
    names(bootR)[3] <- "Method"
  } else if (method == "F") {
    w1 <- which(names %in% c("ContOne", "ZslopeOne"))
    Fconf <- vector("list", length = length(w1) + 1)
    names(Fconf) <- c(names[w1], "Method")
    k <- 1
    
    for (j in w1) {
      n1 <- names[[j]]
      if (n1 == "ContOne") {
        fct = fct_cont_one
      }
      else {
        fct = fct_zslope_one
      }
      
      mods <- object[[1]][[j]]
      x <- object[[4]]$A
      y <- object[[4]]$S
      x1 <- seq(min(x), max(x), interval)
      mod <- object[[1]][[j]]
      res.lm <- summary(mod)
      left.f <- deparse(formula(mod)[[3]])
      S <- sapply(x1, fct, x = x, y = y)
      s2.opt <- res.lm$sigma^2
      S.opt <- min(S)
      Fstat <- (S - S.opt) / s2.opt
      alpha.ci <- 0.05
      z <- qf(1-alpha.ci, 1, res.lm$df[3])
      CI <- which(Fstat <= z)
      CI.lower <- x1[min(CI)]
      CI.upper <- x1[max(CI)]
      Fconf[[k]] <- c(CI.lower, CI.upper)
      k <- k + 1
    }
    bootR <- Fconf
    bootR[[(length(w1) + 1)]] <- "F"
  } else{
    stop("method should be one of 'boot' or 'F'")
  }
  class(bootR) <- "sars"
  attr(bootR, 'type') <- 'threshold_ci'
  return(bootR)
}

#################################################################################
#' Calculate the intercepts and slopes of the different segments
#'
#' @description Calculate the intercepts and slopes of the different segments in
#'   any of the fitted breakpoint regression models available in the package.
#' @usage get_coef(fit)
#' @param fit An object of class 'thresholds', generated using the
#'   \code{\link{sar_threshold}} function.
#' @details The coefficients in the fitted breakpoint regression models do not
#'   all represent the intercepts and slopes of the different segments; to get
#'   these it is necessary to add different coefficients together.
#' @return A dataframe with the intercepts (ci) and slopes (zi) of all segments
#'   in each fitted model. The numbers attached to c and z relate to the
#'   segment, e.g. c1 and z1 are the intercept and slope of the first segment.
#'   For the left-horizontal models, the slope of the first segment (i.e. the
#'   horizontal segment) is not returned. NA values represent cases where a
#'   given parameter is not present in a particular model.
#' @examples
#' data(aegean2)
#' a2 <- aegean2[1:168,]
#' fitT <- sar_threshold(data = a2, mod = c("ContOne", "DiscOne", "ZslopeOne"),
#' interval = 0.1, non_th_models = TRUE, logAxes = "area", logT = log10)
#' #get the slopes and intercepts for these three models
#' coefs <- get_coef(fitT)
#' coefs
#' @importFrom stats coef
#' @export

get_coef <- function(fit){
  if (!"threshold" %in% class(fit))
    stop("fit object should be of class 'threshold'")
  names <- fit[[2]]
  wn <- which(names %in% c("ContOne", "ZslopeOne", 
                           "DiscOne", "ContTwo", "ZslopeTwo", "DiscTwo"))
  resM <- matrix(nrow = length(wn), ncol = 6)
  colnames(resM) <-  c("c1", "z1", "c2", "z2", "c3", "z3")
  rownames(resM) <- names[wn]
  k <- 1
  
  if ("ContOne" %in% names) {
    w <- which(names == "ContOne")
    ContOne <- fit[[1]][[w]]
    c1 <- coef(ContOne)[1]
    z1 <- coef(ContOne)[2]
    z2 <- coef(ContOne)[2] + coef(ContOne)[3]
    resM[k,] <- c(c1, z1, NA, z2, NA, NA)
    k <- k + 1
  }
  
  #one breakpoint zero slope
  if ("ZslopeOne" %in% names) {
    w <- which(names == "ZslopeOne")
    ZslopeOne <- fit[[1]][[w]]
    c1 <- coef(ZslopeOne)[1]
    z2 <- coef(ZslopeOne)[2]
    resM[k,] <- c(c1, NA, NA, z2, NA, NA)
    k <- k + 1
  }
  
  #one breakpoint discontinuous
  if ("DiscOne" %in% names) {
    w <- which(names == "DiscOne")
    DiscOne <- fit[[1]][[w]]
    c1 <- coef(DiscOne)[1] + coef(DiscOne)[3]
    z1 <- coef(DiscOne)[2] + coef(DiscOne)[5]
    c2 <- coef(DiscOne)[1]
    z2 <- coef(DiscOne)[2]
    resM[k,] <- c(c1, z1, c2, z2, NA, NA)
    k <- k + 1
  }
  #two threshold continuous
  if ("ContTwo" %in% names) {
    w <- which(names == "ContTwo")
    ContTwo <- fit[[1]][[w]]
    c1 <- coef(ContTwo)[1]
    z1 <- coef(ContTwo)[2]
    z2 <- coef(ContTwo)[2] + coef(ContTwo)[3]
    z3 <- coef(ContTwo)[2] + coef(ContTwo)[3] + coef(ContTwo)[4]
    resM[k,] <- c(c1, z1, NA, z2, NA, z3)
    k <- k + 1
  }
  #two threshold zero slope
  if ("ZslopeTwo" %in% names) {
    w <- which(names == "ZslopeTwo")
    ZslopeTwo <- fit[[1]][[w]]
    c1 <- coef(ZslopeTwo)[1]
    z2 <- coef(ZslopeTwo)[2]
    z3 <- coef(ZslopeTwo)[2] + coef(ZslopeTwo)[3]
    resM[k,] <- c(c1, NA, NA, z2, NA, z3)
    k <- k + 1
  }
  
  #two breakpoint discontinuous
  if ("DiscTwo" %in% names) {
    w <- which(names == "DiscTwo")
    DiscTwo <- fit[[1]][[w]]
    c1 <- coef(DiscTwo)[1] + coef(DiscTwo)[3]
    z1 <- coef(DiscTwo)[2] + coef(DiscTwo)[6]
    c2 <- coef(DiscTwo)[1] + coef(DiscTwo)[4]
    z2 <- coef(DiscTwo)[2] + coef(DiscTwo)[7]
    c3 <- coef(DiscTwo)[1]
    z3 <- coef(DiscTwo)[2]
    resM[k,] <- c(c1, z1, c2, z2, c3, z3)
    k <- k + 1
  }
  
  resM2 <- as.data.frame(round(resM, 2))
  class(resM2) <- c("sars", "data.frame")
  attr(resM2, 'type') <- "threshold_coef"
  return(resM2)
}

