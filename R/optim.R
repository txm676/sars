#### INTERNAL FUNCTION(S)

######################################## optimization function

#' @importFrom stats formula optim shapiro.test ks.test cor.test pt
#' @importFrom nortest lillie.test

rssoptim <- function(model, data, start = NULL, algo = "Nelder-Mead",
                     normaTest = "lillie", homoTest = "cor.fitted"){

  #initial parameters
  if (is.null(start)) {
    start <- model$init(data)
  } else start <- start

  #if outside ranges : rescaling
  for (i in seq_along(start)) {
    if (model$parLim[i] != "R") {
      if (start[i] <= 0) { start[i] <- 0.1 }
    }
    if (model$parLim[i] == "unif") {
      if (start[i] > 1) { start[i] <- .8 }
    }
  }#eo for

  # starting values on the link function scale
  startMod  <-  transLink(start,model$parLim)
  names(startMod) <- model$parNames

  #RSS function
  rssfun <- model$rss.fun

  #optimization (first result)
  res1 <- tryCatch(optim(startMod, rssfun, hessian = TRUE, data = data,
                         method = algo, control = list(maxit = 50000) ),
                   error = function(e){e}
                   )

  #Backtransformation of parameters values
  res1$par  <-  backLink(res1$par,model$parLim)

  #renaming the parameters vector
  names(res1$par) <- model$parNames

  #calculating expected richness
  S.calc <- model$mod.fun(data$A,res1$par)

  #residuals
  residu  <-  as.vector(S.calc - data$S)

  #second result
  res2  <-  list(startvalues=start,data=data,model=model,
                 calculated=S.calc,residuals=residu)

  #Residuals normality test
  #normaTest  <-  tryCatch(list(shapiro =shapiro.test(residu),
  #kolmo = ks.test(residu, "pnorm") , lillie = list(statistic=NA,p.value=NA)),
  #error = function(e) list(shapiro =list(statistic=NA,p.value=NA),
  #kolmo = list(statistic=NA,p.value=NA) , lillie = list(statistic=NA,
  #p.value=NA) )) #lillie.test(residu)

  #Residuals normality test

  l <- data[[2]]

  if (length(l) < 5) {

      warning("The Lilliefors test cannot be performed with less than 5",
              " data points\n")

  }#eo if length

  if (length(l) < 3) {

      warning("The Shapiro test cannot be performed with less than 3",
              " data points\n")

  }#eo if length

  normaTest <- match.arg(normaTest, c("none", "shapiro", "kolmo", "lillie"))
  homoTest <- match.arg(homoTest, c("none","cor.area","cor.fitted"))

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
        normaTest <- "none"
        }

  #Homogeneity of variance

  if (homoTest == "cor.area"){
    homoTest  <- list("test" = "cor.area", tryCatch(cor.test(residu,data$A),
                            error = function(e)list(estimate=NA,p.value=NA)))
  } else if (homoTest == "cor.fitted"){
    homoTest  <- list("test" = "cor.fitted", tryCatch(cor.test(residu,S.calc),
                            error = function(e)list(estimate=NA,p.value=NA)))
  } else homoTest <- "none"

  #R2, AIC, AICc, BIC

  #common vars
  n <- length(l)
  P <- length(model$parLim) + 1  # + 1 for the estimated variance

  #R2 (Kvaleth, 1985, Am. Statistician)
  R2 <-  1 - ( (res1$value) /  sum((data$S - mean(data$S))^2) )

  #R2a (He & Legendre 1996, p724)
  R2a <-  1 - ( ((n-1)*(res1$value)) /
                  ((n-P)*sum((data$S - mean(data$S))^2)) )

  #AIC
  AIC <- n * log(res1$value / n) + 2 * P

  #AICc
  AICc <- n * log(res1$value / n) + 2*P*(n / (n - P - 1))

  #BIC
  BIC <- n *log(res1$value / n) + log(n) * P

  res3 <- list(AIC=AIC, AICc=AICc, BIC=BIC, R2=R2, R2a=R2a)

  #convergence verif -> 71 is R2<=0
  verge <- ifelse(res1$convergence==0, TRUE, FALSE)
  verge <- ifelse(R2 <= 0, FALSE, TRUE)

  res <- c(res1,list(verge=verge),res2,res3)

  #estimates signifiance and confidence interval (95%)

  #constructing a nlsModel object

  formul <- formula(paste("S ~",as.character(model$exp)))
  env <- environment(formul)
  if (is.null(env)){
    env <- parent.frame()
    environment(formul) <- env
  }

  nMod <- tryCatch(stats_nlsModel(formul,data,res1$par),
                   error = function(e)NA)

  if(class(nMod) != "nlsModel"){
    warning(model$name,": singular gradient at parameter estimates:
   no parameters significance and conf. intervals.", call. = FALSE)
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

    #formating the table of estimates, standard eroor, t value and
    #significance of parameters
    se <- sqrt(diag(XtXinv) * resvar)
    tval <- res1$par/se
    param <- cbind(res1$par, se, tval, 2 * pt(abs(tval), rdf,
                                              lower.tail = FALSE))
    dimnames(param) <- list(model$paramnames, c("Estimate", "Std. Error",
                                                "t value", "Pr(>|t|)"))

    #95% confidence interval
    conf <- matrix(c(param[,"Estimate"] - 2 * param[,"Std. Error"],
                     param[,"Estimate"] + 2 * param[,"Std. Error"]),p,2)
    colnames(conf) <- c("2.5%","97.5%")

    sigConf <- cbind(param,conf)

    res$sigConf <- sigConf

  }#eo if else is.na(nMod)


  res$normaTest <- normaTest
  res$homoTest <- homoTest

  invisible(res)

}#eo rssoptim

######################################## Multiple starting values optimization function

#' @importFrom stats runif

grid_start_fit <- function(model, data, n, algo = "Nelder-Mead",
                           normaTest = "lillie", homoTest = "cor.fitted",
                           verb = TRUE) {

  ns <- ifelse(length(model$parNames) < 4, 100, 10)

  start.list <- lapply(model$parLim,function(x){
    switch(x,
        R = sample(seq(-500, 500), ns),
        Rplus = seq(.1, 500, length.out = ns),
        unif = runif(ns)
    )
  })

  names(start.list) <- model$parNames

  grid.start <- expand.grid(start.list)

  grid.start <- grid.start[sample.int(dim(grid.start)[1],n),]

  if (verb) cat("- running grid optim: ")

  fit.list <- apply(grid.start, 1, function(x){
    if (verb) cat(".")
    tryCatch(rssoptim(model, data , start = x, algo = algo,
                      normaTest = normaTest, homoTest = homoTest)
             , error = function(e) list(value = NA))
  })

  fit.list <- as.list(fit.list)

  values <- unlist(lapply(fit.list,function(x){x$value}))

  min <- which.min(values)

  if (length(min)) {
    fit.list[[min]]
  } else{
    list(value = NA)
  }

}#eo grid_start_fit

######################################## optimization wrapper
get_fit <- function(model = model, data = data, start = NULL,
                    grid_start = NULL, algo = "Nelder-Mead",
                    normaTest = "lillie", homoTest = "cor.fitted",
                    verb = TRUE){

  if (!is.null(start) & !is.null(grid_start)){
    stop("You must choose between 'start' and 'grid_start',",
         " but choose wisely\n")
  }

  if (is.null(start)) {
    fit <- tryCatch(rssoptim(model = model, data = data, algo = algo,
                             normaTest = normaTest, homoTest = homoTest)
                    ,error=function(e) list(value = NA))
    if (is.na(fit$value)) {
      if (!is.null(grid_start)){
        if (grid_start != FALSE){
          n <- min(grid_start, 1000)
          fit <- grid_start_fit(model = model, data = data, n = n,
                                algo = algo, normaTest = normaTest,
                                homoTest = homoTest,verb = verb)
          grid_start <- NULL
        }
      }
    }
  }

  if (!is.null(start)) {
    fit <- tryCatch(rssoptim(model = model, data = data,
                             start = start, algo = algo,
                             normaTest = normaTest, homoTest = homoTest),
                    error = function(e) list(value = NA))
  }

  if (!is.null(grid_start)) {
    if (grid_start != FALSE){
      fit <- grid_start_fit(model = model, data = data, n = grid_start,
                            algo = algo, normaTest = normaTest,
                            homoTest = homoTest, verb = verb)
    }#eo if (grid_start != FALSE)
  }#eo if(!is.null(grid_start))

  if (is.na(fit$value)) {
    warning("The model could not be fitted :(\n")
    return(list(value = NA))
  } else {
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
