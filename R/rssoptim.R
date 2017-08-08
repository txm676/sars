rssoptim <- function(model,data,custstart=NULL,normtest,algo="Nelder-Mead"){

  #paramters bounds
  parLim <- model$parLim

  #initial parameters
  if(is.null(custstart)){
    start <- model$init(data)
  }else{
    start <- custstart
  }

  #if outside ranges : rescaling
  for (i in 1:length(start)) {
    if(parLim[i]!="R"){
      if(start[i]<=0){ start[i] <- 0.1 }
    }
    if(parLim[i]=="unif"){
      if(start[i]>1){ start[i] <- .8 }
    }
  }#eo for

  #sarting values on the link function scale
  startMod  <-  transLink(start,parLim)
  names(startMod) <- model$parNames

  #RSS function
  rssfun <- model$rss.fun

  #optimization (first result)
  res1 <- optim(startMod,rssfun,hessian=F,data=data,method=algo,control=list(maxit=50000))

  #Backtransformation of parameters values
  res1$par  <-  backLink(res1$par,parLim)

  #renaming the parameters vector
  names(res1$par) <- model$parNames

  #calculating expected richness
  S.calc <- model$mod.fun(data$A,res1$par)

  #residuals
  residu  <-  as.vector(data$S - S.calc)

  #second result
  res2  <-  list(startvalues=start,data=data,model=model,calculated=S.calc,residuals=residu)

  #Residuals normality test
  #normaTest  <-  tryCatch(list(shapiro =shapiro.test(residu),kolmo = ks.test(residu, "pnorm") , lillie = list(statistic=NA,p.value=NA) ),  error = function(e) list(shapiro =list(statistic=NA,p.value=NA),kolmo = list(statistic=NA,p.value=NA) , lillie = list(statistic=NA,p.value=NA) )) #lillie.test(residu)

  #Residuals normality test

  if(normtest=="lillie"){

    if(length(l)<5) {

      error("The Lilliefors test cannot be used with less than 5 data points -> switch to the Shapiro test (3 data points still required) \n")

      }#eo if length
  }#eo if lillie

  if(normtest=="shapiro"){

    if(length(l)<3) {

      error("The Shapiro test cannot be used with less than 3 data points -> switch to 'kolmo' or 'none' \n")

    }#eo if length
  }#eo if shapiro

  normaTest <- switch(normtest, "shapiro" = shapiro.test(residu) , "lillie" = lillie.test(residu) , "kolmo" = ks.test(residu, "pnorm"), "none" = list(statistic=NA,p.value=NA) )

  #Homogeneity of variance
  homoTest  <-  tryCatch(list(cor.area = cor.test(residu,data$data$A),cor.fitted = cor.test(residu,S.calc)), error = function(e) list(cor.area = list(estimate=NA,p.value=NA),cor.fitted = list(estimate=NA,p.value=NA)))

  #R2, AIC, AICc, BIC

  #common vars
  n <- length(data$A)
  P <- length(model$parNames) + 1  # + 1 for the estimated variance

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
  verge <- ifelse(R2<=0,71,69)

  res <- c(res1,list(verge=verge,normaTest=normaTest,homoTest=homoTest),res2,res3)

  invisible(res)

}#eo rssoptim
