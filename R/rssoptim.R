rssoptim <- function(model,data,start=NULL, algo = "Nelder-Mead"){

  #initial parameters
  if (is.null(start)) {
    start <- model$init(data)
  }else{
    start <- start
  }

  #if outside ranges : rescaling
  for (i in 1:length(start)) {
    if (model$parLim[i] != "R") {
      if (start[i] <= 0) { start[i] <- 0.1 }
    }
    if (model$parLim[i] == "unif") {
      if (start[i] > 1) { start[i] <- .8 }
    }
  }#eo for

  #sarting values on the link function scale
  startMod  <-  transLink(start,model$parLim)
  names(startMod) <- model$parNames

  #RSS function
  rssfun <- model$rss.fun

  #optimization (first result)
  res1 <- tryCatch(optim(startMod, rssfun, hessian = TRUE, data = data, method = algo, control = list(maxit = 50000) ),
                   error = function(e){e}
                   )

  #Backtransformation of parameters values
  res1$par  <-  sars:::backLink(res1$par,model$parLim)

  #renaming the parameters vector
  names(res1$par) <- model$parNames

  #calculating expected richness
  S.calc <- model$mod.fun(data$A,res1$par)

  #residuals
  residu  <-  as.vector(S.calc - data$S)

  #second result
  res2  <-  list(startvalues=start,data=data,model=model,calculated=S.calc,residuals=residu)

  #Residuals normality test
  #normaTest  <-  tryCatch(list(shapiro =shapiro.test(residu),kolmo = ks.test(residu, "pnorm") , lillie = list(statistic=NA,p.value=NA) ),  error = function(e) list(shapiro =list(statistic=NA,p.value=NA),kolmo = list(statistic=NA,p.value=NA) , lillie = list(statistic=NA,p.value=NA) )) #lillie.test(residu)

  #Residuals normality test

  l <- data[[2]]

  if(length(l)<5) {

      warning("The Lilliefors test cannot be performed with less than 5 data points \n")

  }#eo if length

  if(length(l)<3) {

      warning("The Shapiro test cannot be performed with less than 3 data points \n")

  }#eo if length

  normaTest <- list(shapiro = tryCatch(shapiro.test(residu), error = function(e)NA),
                     lillie = tryCatch(nortest::lillie.test(residu), error = function(e)NA),
                      kolmo = tryCatch(ks.test(residu, "pnorm"), error = function(e)NA)
                    )
  
  #Homogeneity of variance
  
  homoTest  <- list(cor.area = tryCatch(cor.test(residu,data$A), error = function(e)list(estimate=NA,p.value=NA)),
                    cor.fitted = tryCatch(cor.test(residu,S.calc), error = function(e)list(estimate=NA,p.value=NA))
  )

  #R2, AIC, AICc, BIC

  #common vars
  n <- length(l)
  P <- length(model$parLim) + 1  # + 1 for the estimated variance

  #R2 (Kvaleth, 1985, Am. Statistician)
  R2 <-  1 - ( (res1$value) /  sum((data$S - mean(data$S))^2) )

  #R2a (He & Legendre 1996, p724)
  R2a <-  1 - ( ((n-1)*(res1$value)) /  ((n-P)*sum((data$S - mean(data$S))^2)) )

  #AIC
  AIC <- n * log(res1$value / n) + 2 * P

  #AICc
  AICc <- n * log(res1$value / n) + 2*P*(n / (n - P - 1))

  #BIC
  BIC <- n *log(res1$value / n) + log(n) * P

  res3 = list(AIC=AIC, AICc=AICc, BIC=BIC, R2=R2, R2a=R2a)

  #convergence verif -> 71 is R2<=0
  verge <- ifelse(res1$convergence==0, TRUE, FALSE)
  verge <- ifelse(R2 <= 0, FALSE, TRUE)

  res <- c(res1,list(verge=verge),res2,res3)

  #estimates signifiance and confidence interval (95%)

  #constructing a nlsModel object

  formul <- formula(gsub("==","~",as.character(model$formula)))
  env <- environment(formul)
  if (is.null(env)){
    env <- parent.frame()
    environment(formul) <- env
  }

  nMod <- stats:::nlsModel(formul,data,res1$par)

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

  #formating the table of estimates, standard eroor, t value and significance of parameters
  se <- sqrt(diag(XtXinv) * resvar)
  tval <- res1$par/se
  param <- cbind(res1$par, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
  dimnames(param) <- list(model$paramnames, c("Estimate", "Std. Error",
                                              "t value", "Pr(>|t|)"))

  #95% confidence interval
  conf <- matrix(c(param[,"Estimate"] - 2 * param[,"Std. Error"], param[,"Estimate"] + 2 * param[,"Std. Error"]),p,2)
  colnames(conf) <- c("2.5%","97.5%")

  sigConf <- cbind(param,conf)
    #sigConf <- matrix(NA,(P-1),6)
    #colnames(sigConf) <- c("Estimate", "Std. Error","t value", "Pr(>|t|)","2.5%","97.5%")
  #
  res$sigConf <- sigConf
   
  invisible(res)

}#eo rssoptim
