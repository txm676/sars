#### INTERNAL FUNCTION(S)

######################################## optimization function

#' @importFrom stats formula optim shapiro.test ks.test cor.test pt confint coef
#' @importFrom nortest lillie.test

rssoptim <- function(model, data, start = NULL, algo = "Nelder-Mead",
                     normaTest = "none", homoTest = "none",
                     homoCor = "spearman", verb = TRUE){

  #initial parameters
  if (is.null(start)) {
    start <- model$init(data)
  }else{
    start <- start
  }

  #if outside ranges : rescaling
  for (i in seq_along(start)) {
    if (model$parLim[i] != "R") {
      if (start[i] <= 0) { start[i] <- 0.1 }
    }
    if (model$parLim[i] == "unif") {
      if (start[i] > 1) { start[i] <- .8 }
    }
  }#eo for

  #starting values on the link function scale
  startMod  <-  transLink(start,model$parLim)
  names(startMod) <- model$parNames

  #RSS function
  rssfun <- model$rss.fun

  #optimization (first result)
  res1 <- tryCatch(optim(startMod, rssfun, hessian = TRUE, data = data,
                         method = algo, control = list(maxit = 50000) ),
                   error = function(e){e}
                   )

  #Backtransformation of parameter values; the rssfun backtransforms
  #the values tried by optim for Rplus parameters before calculating rss
  #and so we need to do this here to get the true value that was used to get rss.
  res1$par  <-  backLink(res1$par,model$parLim)

  #renaming the parameters vector
  names(res1$par) <- model$parNames
  
  #for asymptotic model, if a negative z is returned, just error, as this
  #causes very strange fits. This will then return the model did not fit
  #message to the user.
  if (model$name == "Asymptotic regression"){
    if (res1$par[3] < 0){
      stop ("Asymp negative z")
    }
  }

  #calculating expected richness
  S.calc <- model$mod.fun(data$A,res1$par)

  #residuals (changed order Nov 2020)
  #residu  <-  as.vector(S.calc - data$S)
  residu  <-  as.vector(data$S - S.calc)
  
  #squared residuals
  sq_residu <- residu^2
  
  #second result
  res2  <-  list(startvalues=start,data=data,model=model,
                 calculated=S.calc,residuals=residu)

  #Residual tests

  normaTest <- match.arg(normaTest, c("none", "shapiro", "kolmo", "lillie"))
  homoTest <- match.arg(homoTest, c("none","cor.area","cor.fitted"))

  l <- data[[2]]
  
  if (length(l) < 5 & normaTest == "lillie") {
    
    warning("The Lilliefors test cannot be performed with less than 5",
            " data points, changing to no residual normality test\n")
    
    normaTest <- "none"
    
  }#eo if length
  
  #normality of residuals
  if (normaTest == "shapiro") {
    normaTest <- list("test" = "shapiro", tryCatch(shapiro.test(residu),
                                                   error = function(e)NA))
    } else if (normaTest == "lillie"){
      normaTest <- list("test" = "lillie", tryCatch(lillie.test(residu),
                                                    error = function(e)NA))
    } else if (normaTest == "kolmo"){
      normaTest <- list("test" = "kolmo", tryCatch(ks.test(residu, "pnorm"),
                                                   error = function(e)NA))
      } else{
        normaTest <- list("test" = "none", "none")
        }

  #Homogeneity of variance

  if (homoTest == "cor.area"){
    homoTest  <- list("test" = "cor.area", 
                      tryCatch(suppressWarnings(cor.test(sq_residu,data$A, 
                      method = homoCor)), 
                      error = function(e)list(estimate=NA,p.value=NA)))
  } else if (homoTest == "cor.fitted"){
    homoTest  <- list("test" = "cor.fitted", 
                      tryCatch(suppressWarnings(cor.test(sq_residu,S.calc,
                      method = homoCor)),
                    error = function(e)list(estimate=NA,p.value=NA)))
  } else {
    homoTest <- list("test" = "none", "none")
  }

  #R2, AIC, AICc, BIC

  #common vars
  n <- length(l)
  P <- length(model$parLim) + 1  # + 1 for the estimated variance

  #R2 (Kvaleth, 1985, Am. Statistician)
  R2 <-  1 - ( (res1$value) /  sum((data$S - mean(data$S))^2) )

  #R2a (He & Legendre 1996, p724)
  R2a <-  1 - ( ((n-1)*(res1$value)) /
                  ((n-P)*sum((data$S - mean(data$S))^2)) )

  #old formulas based on rss
  #AIC <- n * log(res1$value / n) + 2 * P
  #AICc <- n * log(res1$value / n) + 2*P*(n / (n - P - 1))
  #BIC <- n *log(res1$value / n) + log(n) * P
  
  #using log likelihood (and AIC etc) based on nls approach
  #log likelihood code from nls
  val <-  -n * (log(2 * pi) + 1 - log(n) + log(sum(residu^2)))/2
  AIC <- (2 * P) - (2 * val)
  AICc <- -2 * val + 2 * P * (n / (n - P - 1))
  BIC <- (-2 * val) + (P * log(n))
  
  res3 <- list(AIC=AIC, AICc=AICc, BIC=BIC, R2=R2, R2a=R2a)

  verge <- ifelse(res1$convergence==0, TRUE, FALSE)
  #Removed Nov 2020
  #(Korvath - negative R2 indicates complete lack of fit)
  #verge <- ifelse(R2 <= 0, FALSE, TRUE)

  res <- c(res1,list(verge=verge),res2,res3)

  #estimates significance and confidence interval (95%) using nls;
  #but using our fitted parameter estimates rather than re-fitting a new
  #model using nls (i.e. the par estimates are identical), we are simply
  #using the nls framework to get se, t-values, p-values etc.

  #constructing a nlsModel object
  formul <- formula(paste("S ~",as.character(model$exp)))
  env <- environment(formul)
  if (is.null(env)){
    env <- parent.frame()
    environment(formul) <- env
  }

  nMod <- tryCatch(stats_nlsModel(formul,data,res1$par),
                   error = function(e)NA)

  if (!inherits(nMod, "nlsModel")){
    if (verb){
    warning(model$name,": singular gradient at parameter estimates:
             no parameters significance and conf. intervals.", call. = FALSE)
    }
    res$sigConf <- NA
  }else{
    #number of parameters
    p <- length(model$parLim)

    #residuals degrees of freedom
    rdf <- n - p

    #residuals variance
    resvar <- res1$value / rdf

    #calculating the inverse of the upper triangular factor
    #of the gradient array at estimated parameter values
    XtXinv <- chol2inv(nMod$Rmat())
    dimnames(XtXinv) <- list(names(start), names(start))

    #formatting the table of estimates, standard error, t value and
    #significance of parameters
    se <- sqrt(diag(XtXinv) * resvar)
    tval <- res1$par/se
    param <- cbind(res1$par, se, tval, 2 * pt(abs(tval), rdf,
                                              lower.tail = FALSE))
    dimnames(param) <- list(model$paramnames, c("Estimate", "Std. Error",
                                                "t value", "Pr(>|t|)"))
    
    #95% confidence interval, simply 2 * the SE: method used
    #in mmSAR
    conf <- matrix(c(param[,"Estimate"] - 2 * param[,"Std. Error"],
                     param[,"Estimate"] + 2 * param[,"Std. Error"]),p,2)
    colnames(conf) <- c("2.5%","97.5%")

    sigConf <- cbind(param, conf)
    
    #If power model is fitted, properly fit the model using nls and get the par
    #estimates & confidence intervals using their approach. Sometimes the model
    #will fit with nls and converge, but the confint function will not work with
    #the resultant object (e.g. states singular gradient), so if this happens we
    #return nothing for the nls CIs
    if (model$name == "Power"){
      yobs <- data$S
      xobs <- data$A
      nls_pow <- tryCatch(nls(yobs ~ c * xobs^z,
                     start = list("c" = res1$par[1], "z" = res1$par[2])),
      error = function(e)NA)
      if (!anyNA(nls_pow)) {
        coef_nls_pow <- as.vector(coef(nls_pow))
        CI_nls_pow <- tryCatch(suppressMessages(confint(nls_pow)),
                               error = function(e)NA)
        if (!anyNA(CI_nls_pow)){
          sigConf <- cbind(sigConf, "nls.Est." = coef_nls_pow, CI_nls_pow)
          colnames(sigConf)[8:9] <- c("nls.2.5%", "nls.97.5%")
        }
      }
    }#eo if power
    res$sigConf <- sigConf
  }#eo if else is.na(nMod)

  res$normaTest <- normaTest
  res$homoTest <- homoTest

  invisible(res)

}#eo rssoptim

######################################## Multiple starting values optimization function

#function to work out if all values in a vector are the same
vec.equal <- function (x, tol = .Machine$double.eps^0.5) 
{
  if (length(x) == 1) 
    return(TRUE)
  x <- range(x)/mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

#' @importFrom stats runif sd

grid_start_fit <- function(model, data, n, type, algo = "Nelder-Mead",
                           normaTest = "none", homoTest = "none",
                           homoCor = "spearman", 
                           verb = TRUE) {
  
  #if type == "partial", just sample the 500 values from the start.list
  #sequence values.
  
  #  if(length(model$parNames)<4){
  ns <- 100
  # }else{
  #   ns <- 10
  # }
  

  start.list <- lapply(model$parLim,function(x){
    res <- switch(x,
                  R = sample(seq(-500,500), ns),
                  Rplus = c(seq(.1,500,length.out = ns)),
                  unif = runif(ns)
    )
    return(res)
  })
  

  names(start.list) <- model$parNames
  
  grid.start <- expand.grid(start.list)
  
  if (n > nrow(grid.start)){
    n <- nrow(grid.start)
    warning(paste0("grid_n larger than possible combinations for ", model$name,
                   ": setting grid_n to max possible\n"))
  }
  
  grid.start <- grid.start[sample.int(dim(grid.start)[1],n),]
  
  #add our default in as possible 
  def.start <- model$init(data)
  names(def.start) <- colnames(grid.start)
  grid.start <- rbind(grid.start, def.start)
  
  #if type == exhaustive, then do a more expansive search of starting par
  #space.
  if (type == "exhaustive"){
  
  #ensure very small values are included as useful for some models
  sm_val <- lapply(model$parNames, function(x){
    c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1)
  })

  sm_grid <- expand.grid(sm_val)

  colnames(sm_grid) <- colnames(grid.start)

  grid.start <- rbind(grid.start, sm_grid)


# the asymptote parameters are all Rplus, and when e.g. PD is used as the
# response, the asymptote can occur much higher than 500, so for Rplus,
# we tag on the 4 largest richness values on the end (and one 75% of largest).
# Values also tagged on
# if max richness < 500, but this doesn't matter (and might help as should
# be closer to true asymptote value)
if (any(model$parLim == "Rplus")){
  RPM <- sort(data$S, decreasing = TRUE)[1:4]
  RPM75 <- max(data$S) * 0.75
  RPM <- c(RPM, RPM75)
  RPM <- c(RPM)
  WPM <- which(model$parLim == "Rplus")
  WPM2 <- which(!model$parLim == "Rplus")

  ZZ <- vector("list", length = length(model$parLim))
  #iterate across model parameters, and if Rplus store the 4 largest values,
  #and if not, store the relevant values. Then create a new expanded grid
  #and add onto grid.start
  for (i in 1:length(model$parLim)){
    if (model$parLim[i] == "Rplus"){
      ZZ[[i]] <- RPM
    } else if(model$parLim[i] == "R"){
      ZZ[[i]] <- c(-500, -250, -50, 0.001, 0.01, 0.1, 1, 50, 250, 500)
    } else{
      ZZ[[i]] <- runif(5)
    }
  }#eo for
  lar_grid <- expand.grid(ZZ)
  colnames(lar_grid) <- colnames(grid.start)
  grid.start <- rbind(grid.start, lar_grid)
}#eo if Rplus

#some more specific z values for Chapman and Gompertz models
if (model$name == "Chapman Richards"){
  zseq <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
            0.01, 0.05, 0.1, 0.5)
  gs2 <- data.frame(rep(def.start[1], 10), zseq, rep(def.start[3], 10))
  colnames(gs2) <- colnames(grid.start)
  grid.start <- rbind(grid.start, gs2)
  #Some extra c values, based on holding z and d constant
  z1 <-  1 / sd(data$A)#stack forum presents this as useful z start
  c1 <- seq(0.01, 25, 0.001)
  c1 <- c(c1, seq(25, 100, 0.01))
  gs3 <- data.frame(rep(def.start[1], length(c1)), 
                    rep(z1, length(c1)), c1)
  colnames(gs3) <- colnames(grid.start)
  grid.start <- rbind(grid.start, gs3)
}
if (model$name == "Gompertz"){
  d2 <- sort(data$S, decreasing = TRUE)[1:3]
  zz2 <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005,
            0.01, 0.05, 0.1, 0.5)
  cc <- def.start[3]
  cc2 <- c(cc + 5, cc + 10, cc + 50, cc + 100, cc + 200, cc + 500, cc + 800,
           cc - 5, cc - 10, cc - 50, cc - 100, cc - 200, cc - 500, cc - 800)
  gs2 <- expand.grid(d2, zz2, cc2)
  colnames(gs2) <- colnames(grid.start)
  grid.start <- rbind(grid.start, gs2)
}
  
}#eo if exhaustive

  #####################################################
  
  fit.list <- suppressWarnings(apply(grid.start, 1, function(x){
    #  if (verb) cat(".")
    tryCatch(rssoptim(model, data, start = x, algo = algo,
                      normaTest = normaTest, homoTest = homoTest,
                      homoCor = homoCor), 
             error = function(e) list(value = NA))
  }))
  
  fit.list <- as.list(fit.list)
  
  values <- unlist(lapply(fit.list,function(x){x$value}))
  
  # note this just returns one value even if there are multiple values with
  # the same lowest rss - and there almost always will be as lots of starting
  # par estimates will converge on same final pars, so this is fine.
  
  #This finds the min RSS, then checks if multiple par estimates return the same
  #rss. If so, it throws a warning and randomly selects a set. If not,
  #it just returns the min set
  min_val <- min(values, na.rm = TRUE)
  w_mult <- as.vector(which(values == min_val))
  if (length(w_mult) > 1){
    mult_pars <- vapply(fit.list[w_mult], function(x) x$par, 
         FUN.VALUE = numeric(length(model$parNames)))
    mult_pars <- round(mult_pars, 2)
    mult_equal <- apply(mult_pars, 2, function(x) vec.equal(x))
    if (any(!mult_equal) & verb){
      warning(paste0(model$name, ": Multiple parameter estimates returned the ", 
                     "same minimum rss; one set have been randomly selected"))
    }
    min <- sample(w_mult, 1)
  } else{
    min <- which.min(values)
  }
  
  if (length(min) != 0) {
    fit.list[[min]]
  } else{
    list(value = NA)
  }
  
}#eo grid_start_fit

######################################## optimization wrapper
get_fit <- function(model = model, data = data, start = NULL,
                    grid_start = "partial", grid_n = NULL, algo = "Nelder-Mead",
                    normaTest = "none", homoTest = "none", 
                    homoCor = "spearman",
                    verb = TRUE){
  
  if (isFALSE(is.null(start)) & (grid_start != "none")){
    stop("You must choose between 'start' and 'grid_start',",
         " but choose wisely\n")
  }
  
  ##found that for some models, grid start often pushes it into weird parameter 
  ##space (e.g. weibull 3 becomes wiggly), but not forbidden space (e.g. pars 
  ##are all positive) so for now, for these models, just fit using our starting
  ##parameter estimates.
  if (grid_start != "none") { #use grid_search
    
    if (!model$name %in% c("Cumulative Weibull 3 par.",
                           "Cumulative Weibull 4 par.")){
      #for grid_start == partial, use n of 500.
      if (grid_start == "partial") grid_n <- 500
    
    fit <- grid_start_fit(model = model, data = data, n = grid_n, 
                          type = grid_start,
                          algo = algo, normaTest = normaTest,
                          homoTest = homoTest, homoCor = homoCor, verb = verb)
    } else{
      fit <- tryCatch(rssoptim(model = model, data = data, algo = algo,
                               normaTest = normaTest, homoTest = homoTest,
                               homoCor = homoCor, verb = verb),
                      error=function(e) list(value = NA))
    }
  } else if (!is.null(start)){#or provided start values
    fit <- tryCatch(rssoptim(model = model, data = data,
                             start = start, algo = algo,
                             normaTest = normaTest, homoTest = homoTest,
                             homoCor = homoCor, verb = verb),
                    error = function(e) list(value = NA))
  } else { #or if neither selected, use default start values from within sars
    fit <- tryCatch(rssoptim(model = model, data = data, algo = algo,
                             normaTest = normaTest, homoTest = homoTest,
                             homoCor = homoCor, verb = verb),
                    error=function(e) list(value = NA))
  }
  
  if(is.na(fit$value) & verb){
    warning("The model could not be fitted :(\n")
    return(list(value = NA))
  }else{
    return(fit)
  }
}#eo get_fit



###################################################################
##########stats:::nlsModel###########################
####################################################

#this is the stats:::nlsModel function. It needs to be manually included
#as CRAN does not allow :::

#' @importFrom stats start numericDeriv


stats_nlsModel <- function (form, data, start, wts, upper = NULL)
{
  thisEnv <- environment()
  env <- new.env(hash = TRUE, parent = environment(form))
  for (i in names(data)) assign(i, data[[i]], envir = env)
  ind <- as.list(start)
  parLength <- 0
  for (i in names(ind)) {
    temp <- start[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp, envir = env)
    ind[[i]] <- parLength + seq_along(start[[i]])
    parLength <- parLength + length(start[[i]])
  }
  getPars.noVarying <- function() unlist(mget(names(ind), env))
  getPars <- getPars.noVarying
  internalPars <- getPars()
  if (!is.null(upper))
    upper <- rep_len(upper, parLength)
  useParams <- rep_len(TRUE, parLength)
  lhs <- eval(form[[2L]], envir = env)
  rhs <- eval(form[[3L]], envir = env)
  .swts <- if (!missing(wts) && length(wts))
    sqrt(wts)
  else rep_len(1, length(rhs))
  assign(".swts", .swts, envir = env)
  resid <- .swts * (lhs - rhs)
  dev <- sum(resid^2)
  if (is.null(attr(rhs, "gradient"))) {
    getRHS.noVarying <- function() {
      if (is.null(upper))
        numericDeriv(form[[3L]], names(ind), env)
      else numericDeriv(form[[3L]], names(ind), env, ifelse(internalPars <
                                                              upper, 1, -1))
    }
    getRHS <- getRHS.noVarying
    rhs <- getRHS()
  }
  else {
    getRHS.noVarying <- function() eval(form[[3L]], envir = env)
    getRHS <- getRHS.noVarying
  }
  dimGrad <- dim(attr(rhs, "gradient"))
  marg <- length(dimGrad)
  if (marg > 0L) {
    gradSetArgs <- vector("list", marg + 1L)
    for (i in 2L:marg) gradSetArgs[[i]] <- rep_len(TRUE,
                                                   dimGrad[i - 1])
    useParams <- rep_len(TRUE, dimGrad[marg])
  }
  else {
    gradSetArgs <- vector("list", 2L)
    useParams <- rep_len(TRUE, length(attr(rhs, "gradient")))
  }
  npar <- length(useParams)
  gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
  gradCall <- switch(length(gradSetArgs) - 1L,
          call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], drop = FALSE),
          call("[", gradSetArgs[[1L]],  gradSetArgs[[2L]], gradSetArgs[[2L]],
          drop = FALSE), call("[", gradSetArgs[[1L]], gradSetArgs[[2L]],
          gradSetArgs[[2L]], gradSetArgs[[3L]], drop = FALSE),
          call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
               gradSetArgs[[3L]], gradSetArgs[[4L]], drop = FALSE))
  getRHS.varying <- function() {
    ans <- getRHS.noVarying()
    attr(ans, "gradient") <- eval(gradCall)
    ans
  }
  if (length(gr <- attr(rhs, "gradient")) == 1L)
    attr(rhs, "gradient") <- gr <- as.vector(gr)
  QR <- qr(.swts * gr)
  qrDim <- min(dim(QR$qr))
  if (QR$rank < qrDim)
    stop("singular gradient matrix at initial parameter estimates")
  getPars.varying <- function() unlist(mget(names(ind), env))[useParams]
  setPars.noVarying <- function(newPars) {
    assign("internalPars", newPars, envir = thisEnv)
    for (i in names(ind)) assign(i, unname(newPars[ind[[i]]]),
                                 envir = env)
  }
  setPars.varying <- function(newPars) {
    internalPars[useParams] <- newPars
    for (i in names(ind)) assign(i, unname(internalPars[ind[[i]]]),
                                 envir = env)
  }
  setPars <- setPars.noVarying
  on.exit(remove(i, data, parLength, start, temp, m))
  m <- list(resid = function() resid, fitted = function() rhs,
            formula = function() form, deviance = function() dev,
            lhs = function() lhs, gradient = function() .swts * attr(rhs,
            "gradient"), conv = function() {
             if (npar == 0) return(0)
             rr <- qr.qty(QR, resid)
             sqrt(sum(rr[1L:npar]^2)/sum(rr[-(1L:npar)]^2))
             }, incr = function() qr.coef(QR, resid),
            setVarying = function(vary = rep_len(TRUE,
            length(useParams))) {
             assign("useParams", if (is.character(vary)) {
             temp <- logical(length(useParams))
             temp[unlist(ind[vary])] <- TRUE
             temp
            } else if (is.logical(vary) && length(vary) != length(useParams))
              stop("setVarying : 'vary' length must match",
                   " length of parameters") else {
              vary
              }, envir = thisEnv)
              gradCall[[length(gradCall) - 1L]] <<- useParams
              if (all(useParams)) {
              assign("setPars", setPars.noVarying, envir = thisEnv)
              assign("getPars", getPars.noVarying, envir = thisEnv)
              assign("getRHS", getRHS.noVarying, envir = thisEnv)
              assign("npar", length(useParams), envir = thisEnv)
              } else {
             assign("setPars", setPars.varying, envir = thisEnv)
             assign("getPars", getPars.varying, envir = thisEnv)
             assign("getRHS", getRHS.varying, envir = thisEnv)
              assign("npar", length(seq_along(useParams)[useParams]),
              envir = thisEnv)
              }
              }, setPars = function(newPars) {
              setPars(newPars)
              assign("resid", .swts * (lhs - assign("rhs", getRHS(),
              envir = thisEnv)), envir = thisEnv)
              assign("dev", sum(resid^2), envir = thisEnv)
              if (length(gr <- attr(rhs, "gradient")) == 1L) gr <- c(gr)
              assign("QR", qr(.swts * gr), envir = thisEnv)
              (QR$rank < min(dim(QR$qr)))
              }, getPars = function() getPars(),
            getAllPars = function() getPars(),
            getEnv = function() env, trace = function() {
              cat(format(dev), ": ", format(getPars()))
              cat("\n")
            }, Rmat = function() qr.R(QR),
            predict = function(newdata = list(), qr = FALSE) eval(form[[3L]],
                                                      as.list(newdata), env))
  class(m) <- "nlsModel"
  m
}
