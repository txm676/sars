#' @export
#' #LINK TO SAR MULTI HELP FILE
sars_models <- function() {
  c("power","powerR","epm1","epm2","p1","p2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg", "linear")
}


#' Fit a multimodel SAR curve
#'
#' @description Construct a multimodel species-area relationship curve using
#'   information criterion weights and up to twenty SAR models.
#' @usage sar_multi(data, obj = c("power",
#'   "powerR","epm1","epm2","p1","p2","expo","koba",
#'   "mmf","monod","negexpo","chapman","weibull3","asymp",
#'   "ratio","gompertz","weibull4","betap","heleg", "linear"), crit = "Info",
#'   normaTest = "lillie", homoTest = "cor.fitted", neg_check = TRUE,
#'   alpha_normtest = 0.05, alpha_homotest = 0.05, confInt = FALSE, ciN = 100,
#'   verb = TRUE)
#' @param data A dataset in the form of a dataframe with two columns: the first
#'   with island/site areas, and the second with the species richness of each
#'   island/site.
#' @param obj Either a vector of model names or a fit_collection object.
#' @param crit The criterion used to compare models and compute the model
#'   weights. The default \code{crit = "Info"} switches to AIC or AICc depending
#'   on the number of data points in the dataset. For BIC, use \code{crit =
#'   "Bayes"}.
#' @param normaTest The test used to test the normality of the residuals of each
#'   model. Can be any of "lillie" (Lilliefors Kolmogorov-Smirnov test; the
#'   default), "shapiro" (Shapiro-Wilk test of normality), "kolmo"
#'   (Kolmogorov-Smirnov test), or "none" (no residuals normality test is
#'   undertaken).
#' @param homoTest The test used to check for homogeneity of the residuals of
#'   each model. Can be any of "cor.fitted" (a correlation of the residuals with
#'   the model fitted values; the default), "cor.area" (a correlation of the
#'   residuals with the area values), or "none" (no residuals homogeneity test
#'   is undertaken).
#' @param neg_check Whether or not a check should be undertaken to flag any
#'   models that predict negative richness values.
#' @param alpha_normtest The alpha value used in the residual normality test
#'   (default = 0.05, i.e. any test with a P value < 0.05 is flagged as failing
#'   the test).
#' @param alpha_homotest The alpha value used in the residual homogeneity test
#'   (default = 0.05, i.e. any test with a P value < 0.05 is flagged as failing
#'   the test).
#' @param confInt A logical argument specifying whether confidence intervals
#'   should be calculated for the multimodel curve using bootstrapping.
#' @param ciN The number of bootstrap samples to be drawn to calculate the
#'   confidence intervals (if \code{confInt == TRUE}).
#' @param verb verbose (default: \code{verb == TRUE}).
#' @details The multimodel SAR curve is constructed using information criterion
#'   weights (see Burnham & Anderson, 2002; Guilhaumon et al. 2010). If
#'   \code{obj} is a vector of n model names the function fits the n models to
#'   the dataset provided. A dataset must have four or more datapoints to fit
#'   the multimodel curve. If any models cannot be fitted they are removed from
#'   the multimodel SAR. If \code{obj} is a fit_collection object, any model
#'   fits in the collection which are NA are removed. In addition, if any other
#'   model checks have been selected (i.e. residual normality and heterogeneity
#'   tests, and checks for negative predicted richness values), these are
#'   undertaken and any model that fails the selected test(s) is removed from
#'   the multimodel SAR. The order of the additional checks inside the function
#'   is: normality of residuals, homogeneity of residuals, and a check for
#'   negative fitted values. Once a model fails one test it is removed and thus
#'   is not available for further tests. Thus, a model may fail multiple tests
#'   but the returned warning will only provide information on a single test.
#'
#'   The resultant models are then used to construct the multimodel SAR curve.
#'   For each model in turn, the model fitted values are multiplied by the
#'   information criterion weight of that model, and the resultant values are
#'   summed across all models (Burnham & Anderson, 2002). Confidence intervals
#'   can be calculated (using \code{confInt}) around the multimodel averaged
#'   curve using the bootstrap procedure outlined in Guilhaumon et al (2010).The
#'   procedure transforms the residuals from the individual model fits and
#'   occasionally NAs / Inf values can be produced - in these cases, the model
#'   is removed from the confidence interval calculation (but not the multimodel
#'   curve itself). When several SAR models are used and the number of bootstraps
#'   (\code{ciN}) is large, generating the confidence intervals can take a long
#'   time.
#'   
#'   The \code{sar_models()} function can be used to bring up a list of the 20
#'   model names. \code{display_sars_models()} generates a table of the 20 models
#'   with model information.
#'
#' @return A list of class "multi" and class "sars" with two elements. The first
#'   element ('mmi') contains the fitted values of the multimodel sar curve. The
#'   second element ('details') is a list with the following components:
#'   \itemize{ \item{mod_names} { Names of the models that were successfully
#'   fitted and passed any model check} \item{fits} { A fit_collection object
#'   containing the successful model fits} \item{ic} { The information criterion
#'   selected} \item{norm_test} { The residual normality test selected}
#'   \item{homo_test} { The residual homogeneity test selected}
#'   \item{alpha_norm_test} { The alpha value used in the residual normality
#'   test} \item{alpha_homo_test} { The alpha value used in the residual
#'   homogeneity test} \item{ics} { The information criterion values (e.g. AIC
#'   values) of the model fits} \item{delta_ics} { The delta information
#'   criterion values} \item{weights_ics} { The information criterion weights of
#'   each model fit} \item{n_points} {  Number of data points} \item{n_mods} {
#'   The number of successfully fitted models} \item{no_fit} { Names of the
#'   models which could not be fitted or did not pass model checks} }
#'
#'   The \code{\link{summary.sars}} function returns a more useful summary of
#'   the model fit results, and the \code{\link{plot.multi}} plots the
#'   multimodel curve.
#' @note Occasionally a model fit will converge and pass the model fitting
#'   checks (e.g. residual normality) but the resulting fit is nonsensical (e.g.
#'   a horizontal line with intercept at zero). Thus, it can be useful to plot
#'   the resultant 'multi' object to check the individual model fits. To re-run
#'   the \code{sar_multi} function without a particular model, simply remove it
#'   from the \code{obj} argument.
#'
#'   The generation of confidence intervals around the multimodel curve (using
#'   \code{confInt == TRUE}), may throw up errors that we have yet to come
#'   across. Please report any issues to the package maintainer.
#'
#' @references Burnham, K. P., & Anderson, D. R. (2002). Model selection and
#'   multi-model inference: a practical information-theoretic approach (2nd
#'   ed.). New-York: Springer.
#'
#'   Guilhaumon, F., Mouillot, D., & Gimenez, O. (2010). mmSAR: an R-package for
#'   multimodel species-area relationship inference. Ecography, 33, 420-424.
#' @import stats
#' @import cli
#' @import crayon
#' @examples
#' data(galap)
#' #attempt to construct a multimodel SAR curve using all twenty sar models
#' fit <- sar_multi(galap)
#' summary(fit)
#' plot(fit)
#'
#' # construct a multimodel SAR curve using a fit_collection object
#' s1 <- sar_power(galap)
#' s2 <- sar_expo(galap)
#' s3 <- sar_koba(galap)
#' ff <- fit_collection(s1, s2, s3)
#' fit2 <- sar_multi(galap, obj = ff)
#' summary(fit2)
#'
#' # construct a multimodel SAR curve without conducting any model checks
#' fit3 <- sar_multi(galap, normaTest = "none", homoTest = "none", neg_check = FALSE)
#'
#' @export


sar_multi <- function(data,
                       obj = c("power", "powerR","epm1","epm2","p1","p2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg", "linear"),
                       crit = "Info",
                       normaTest = "lillie",
                       homoTest = "cor.fitted",
                       neg_check = TRUE,
                       alpha_normtest = 0.05,
                       alpha_homotest = 0.05,
                       confInt = FALSE,
                       ciN = 100,
                       verb = TRUE){

  if (!((is.character(obj))  | (class(obj) == "sars")) ) stop("obj must be of class character or sars")

  if (nrow(data) < 4) stop("Multi SAR needs at least four data points")
  if (nrow(data) == 4 & normaTest == "lillie") stop("The Lilliefors test cannot be performed with less than 5 data points")

  if (is.character(obj) & is.null(data)) stop("if obj is character then data should be provided")

  if (is.character(obj)) {
    if (any(!(obj %in% c("linear","power","powerR","epm1","epm2","p1","p2","expo","koba","mmf","monod","negexpo","chapman","weibull3","asymp","ratio","gompertz","weibull4","betap","heleg")))) stop("provided model names do not match with model functions")
  }

  if (length(obj) < 2) stop("more than 1 fit is required to construct a sar_multi")

  normaTest <- match.arg(normaTest, c("none", "shapiro", "kolmo", "lillie"))
  homoTest <- match.arg(homoTest, c("none","cor.area","cor.fitted"))

  if (normaTest == "none") alpha_normtest <- "none"
  if (homoTest == "none") alpha_homotest <- "none"

  #if (verb) cat_line(rule(left = paste0(cyan(symbol$bullet),bold(" multi_sars")),right="multi-model SAR"))
  if (verb & is.character(obj)) {
    cat("\n", paste("Now attempting to fit the", length(obj), "SAR models:"), "\n\n")
    cat_line(rule(left = bold(" multi_sars"),right="multi-model SAR"))
  }
  #if (verb) cat_line(magenta(symbol$arrow_right)," Data set is: ")
  #if (verb) cat_line(rule(left = paste0(magenta(symbol$bullet))))
  #if (verb) bullet("O | S : model", bullet = blue_arrow())

  #if not yet fitted, fit the models to the datasquare_small_filled
  if (is.character(obj)) {

    mods <- paste0("sar_",obj)
    names(mods) <- obj

    fits <- suppressWarnings(lapply(obj, function(x){

      f <- eval(parse(text = paste0(mods[x],"(data", ", normaTest = ", paste0("'", normaTest, "'"), ", homoTest = ", paste0("'", homoTest, "'"), ")")))

      if (verb) {
        if(is.na(f$value)) {
          cat_line( paste0(red(symbol$arrow_right)," ",col_align(x,max(nchar(obj)))," : ", red(symbol$cross)))
        }else{

          if (!is.matrix(f$sigConf)){
            cat_line( paste0(yellow(symbol$arrow_right)," ",col_align(x,max(nchar(obj)))," : Warning: could not compute parameters statistics"))
          }else{
            cat_line( paste0(cyan(symbol$arrow_right)," ",col_align(x,max(nchar(obj)))," : ",green(symbol$tick)))
          }
        }
      }

      f

    }))#eo suppressWarnings(lapply)

    #remove models with no parameter estimates
  #  sigC <- vapply(fits, function(x) any(is.na(x$sigConf)), FUN.VALUE = logical(1))
    #if(all(sigC)) stop("No model could be fitted, aborting multi_sars\n")

   # if (any(sigC)){
    #  warning("Could not compute parameter statistics for one or more models and these have been excluded from the multi SAR", call. = FALSE)
    #  badNames2 <- vapply(fits[sigC], FUN = function(x){x$model$name}, FUN.VALUE = character(1))
    #  if (badMods != 0) {
    #    badMods <- c(badMods, badNames2)
    #  } else{
     #   badMods <- badNames2
      #}
      #fits <- fits[!sigC]
   #}
  }else{
    if (attributes(obj)$type == "fit_collection"){
      fits <- obj
      if(all(is.na(fits))){
        stop("The fit collection had no fits, aborting multi_sars\n")
      }
    }else{
      stop("an object of class 'sars' is passed but not of type 'fit_collection':(")
    }
  }#eo if else is.character(obj)

  #####BAD MODEL CHECKS#######################

  #NA CHECKS
  f_nas <- unlist(lapply(fits,function(b)b$value))

  if(all(is.na(f_nas))){
    stop("No model could be fitted, aborting multi_sars\n")
  }

  badMods <- vector(length = 0, mode = "character")

  if(any(is.na(f_nas))){
    badNames <- is.na(f_nas)
    message("\n", paste(sum(is.na(f_nas)), "models could not be fitted and have been excluded from the multi SAR"), "\n")
    badMods <- obj[badNames] #extract the bad model names from the obj vector (not from fits, as no model name if NA)
    fits <- fits[!is.na(f_nas)]
  }

  bml <- length(badMods)


  if ((normaTest != "none" | homoTest != "none" | neg_check) & verb){
    if (is.character(obj)){
      if (any(is.na(f_nas))){
       cat("\nModel fitting completed. Now undertaking model validation checks.\nAdditional models will be excluded if necessary:\n")
      } else {
        cat("\nModel fitting completed - all models succesfully fitted. Now undertaking model validation checks.\nAdditional models will be excluded if necessary:\n")
      }
    } else {
      cat("\nNow undertaking model validation checks. Additional models will\nbe excluded if necessary\n")
    }
  }

  #if checks for normality and / or homoscedasticity enabled, then check and remove bad fits from fits

  if (normaTest != "none") {
    np <- vapply(fits, function(x) x$normaTest[[2]]$p.value, FUN.VALUE = numeric(1))
    #sometimes bad models produce calculated values with all same richness values and no correlation
    #can be done. Remove these
    if (anyNA(np)){
      wnn <- is.na(np)
      mn <- vapply(fits, function(x) x$model$name, FUN.VALUE = character(1))#get all names in fit collection
      message("\n", paste(sum(is.na(np)),"models returned NAs in the residuals normality test and have been excluded from the multi SAR:"), "\n",
              paste(mn[wnn], collapse = ", "), "\n")
      badMods <- c(badMods, mn[wnn])#select the model names with NAs
      fits <- fits[!wnn]
      np <- np[!wnn]
    }
    whp <- np < alpha_normtest
    if (any(whp)) {
      mn2 <- vapply(fits, function(x) x$model$name, FUN.VALUE = character(1))#get all names in fit collection
      message("\n", paste(sum(np < alpha_normtest), "models failed the residuals normality test and have been excluded from the multi SAR:"), "\n",
              paste(mn2[whp], collapse = ", "), "\n")
      badMods <- c(badMods, mn2[whp])#select the model names for models with p < 0.05
      fits <- fits[!whp]#then remove these models from the fit collection
    }
  }
  if (homoTest != "none") {
    hp <- vapply(fits, function(x) x$homoTest[[2]]$p.value, FUN.VALUE = numeric(1))

    #sometimes bad models produce calculated values with all same richness values and no correlation
    #can be done. Remove these
    if (anyNA(hp)){
      whn <- is.na(hp)
      mn3 <- vapply(fits, function(x) x$model$name, FUN.VALUE = character(1))#get all names in fit collection
      message("\n", paste(sum(is.na(hp)),"returned NAs in the residuals homogeneity test and have been excluded from the multi SAR:"), "\n",
              paste(mn3[whn], collapse = ", "), "\n")
      badMods <- c(badMods, mn3[whn])#select the model names with NAs
      fits <- fits[!whn]
      hp <- hp[!whn]
    }
    whh <- hp < alpha_homotest
    if (any(whh)) {
      mn4 <- vapply(fits, function(x) x$model$name, FUN.VALUE = character(1))#get all names in fit collection
      message("\n", paste(sum(hp < alpha_homotest),"models failed the residuals homogeneity test and have been excluded from the multi SAR:"), "\n",
              paste(mn4[whh], collapse = ", "), "\n")
      badMods <- c(badMods, mn4[whh])#select the model names for models with p < 0.05
      fits <- fits[!whh]
    }
  }
  #negative values
  if (neg_check){
    nc <- vapply(fits, function(x) any(x$calculated < 0), FUN.VALUE = logical(1))
    if (any(nc)) {
      mn5 <- vapply(fits, function(x) x$model$name, FUN.VALUE = character(1))#get all names in fit collection
      message("\n", paste(sum(nc), "models have negative fitted values and have been excluded from the multi SAR:"), "\n",
              paste(mn5[nc], collapse = ", "), "\n")
      badMods <- c(badMods, mn5[nc])#select the model names for models with p < 0.05
      fits <- fits[!nc]
    }
  }

  if (length(badMods) == bml & verb) cat("\n", "All models passed the model validation checks","\n","\n")

  if (length(badMods) == 0) badMods <- 0
  if (length(fits) < 2) stop("Fewer than two models could be fitted and / or passed the model checks")


  sf <- vapply(fits, function(x) x$model$name, FUN.VALUE = character(1))
  if (verb) cat(length(sf), "remaining models used to construct the multi SAR:\n",
       paste(sf, collapse = ", "), "\n")

  if (is.character(obj)) fits <- fit_collection(fits = fits)

####################################

  #setting variables
  nPoints <- length(fits[[1]]$data$A)
  nMods <- length(fits)

  modNames <- vapply(fits, FUN = function(x){x$model$name}, FUN.VALUE = character(1))

  #choosing an IC criterion (AIC or AICc or BIC)
  IC <- switch(crit,
         Info= if ( (nPoints / 3) < 40 ) { "AICc" } else { "AIC" },
         Bayes= "BIC"
  )

  #get ICs
  ICs <- vapply(X = fits, FUN = function(x){x[[IC]]}, FUN.VALUE = double(1))

  #get delta ICs
  delta_ICs <- ICs - min(ICs)

  #get akaike weights
  akaikesum <- sum(exp( -0.5*(delta_ICs)))
  weights_ICs <- exp(-0.5*delta_ICs) / akaikesum

  #ERROR: produce weight averaged diversity measures
 # mmi <- vapply(fits,FUN=function(x){x$calculated},FUN.VALUE=double(nPoints))
 # mmi <- apply((mmi * weights_ICs), 1 , sum)

  #produce weight averaged diversity measures
  mmi <- vapply(fits,FUN=function(x){x$calculated},FUN.VALUE=double(nPoints))
  wm <- matrix(nrow = nPoints, ncol = length(fits))
  for (i in seq_along(weights_ICs)) {wm[ ,i] <- mmi[ ,i] * weights_ICs[i]}
  mmi <- apply(wm, 1 , sum)

  res <- mmi

    details <- list(
      mod_names = modNames,
      fits = fits,
      crit = crit,
      ic = IC,
      norm_test = normaTest,
      homo_test = homoTest,
      alpha_norm_test = alpha_normtest,
      alpha_homo_test = alpha_homotest,
      ics = ICs,
      delta_ics = delta_ICs,
      weights_ics = weights_ICs,
      n_points = nPoints,
      n_mods = nMods,
      no_fit = as.vector(badMods)
    )

  res <- list(mmi = mmi, details = details)

  class(res) <- c("multi", "sars")
  attr(res, "type") <- "multi"

  #if (verb) cat_line(rule(left = cyan(symbol$bullet)))
  if (verb) cat_line(rule())

  if (confInt){
    cat("\nCalculating sar_multi confidence intervals - this may take some time:\n")
    cis <- sar_conf_int(res, n = ciN, crit = crit, normaTest = normaTest,
                        homoTest = homoTest,
                        neg_check = neg_check,
                        alpha_normtest = alpha_normtest,
                        alpha_homotest = alpha_homotest, verb = verb)
    res$details$confInt <- cis
  } else {
    res$details$confInt <- NA
  }

  invisible(res)

}#end of multi_sars
