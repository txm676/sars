######################################## optimization function
rssoptim <- function(model, data, start = NULL, algo = "Nelder-Mead"){

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

  nMod <- tryCatch(stats:::nlsModel(formul,data,res1$par), error = function(e)NA)
  
  if(class(nMod) != "nlsModel"){
    warning(model$name,": singular gradient matrix at parameter estimates. Could not compute parameters significance and conf intervals.", call. = FALSE)
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
    
    res$sigConf <- sigConf
    
  }#eo if else is.na(nMod)
  
   
  invisible(res)

}#eo rssoptim

######################################## Multiple starting values optimization function
grid_start_fit <- function(model, data, n, algo = "Nelder-Mead", verb = TRUE) {
  
  if(length(model$parNames)<4){
    ns <- 100
  }else{
    ns <- 10
  }
  
  start.list <- lapply(model$parLim,function(x){
    res = switch(x,
                 R = sample(seq(-500,500),ns),
                 Rplus = seq(.1,500,length.out = ns),
                 unif = runif(ns)
    )
    return(res)
  })
  
  names(start.list) <- model$parNames
  
  grid.start <- expand.grid(start.list)
  
  grid.start <- grid.start[sample.int(dim(grid.start)[1],n),]
  
  if (verb) cat("- running grid optim: ")
  
  fit.list <- apply(grid.start, 1, function(x){
    if (verb) cat(".")
    tryCatch(rssoptim(model, data , start = x, algo = algo), error = function(e) list(value = NA))
  })
  
  fit.list <- as.list(fit.list)
  
  values <- unlist(lapply(fit.list,function(x){x$value}))
  
  min <- which.min(values)
  
  if(length(min) != 0) {
    fit.list[[min]]
  }else{
    list(value = NA)
  }

}#eo grid_start_fit

######################################## optimization wrapper
get_fit <- function(model = model, data = data, start = NULL, grid_start = NULL, algo = "Nelder-Mead", verb = TRUE){
  
  if(!is.null(start) & !is.null(grid_start)){
    stop("You must choose between 'start' and 'grid_start', but choose wisely\n")
  }
  
  if(is.null(start)){
    fit <- tryCatch(rssoptim(model = model, data = data, algo = algo),error=function(e) list(value = NA))
    if(is.na(fit$value)){
      if(!is.null(grid_start)){
        if(grid_start != FALSE){
          n <- min(grid_start,1000)
          fit <- grid_start_fit(model = model, data = data, n = n, algo = algo, verb = verb)
          grid_start <- NULL
        }
      }
    }
  }
  
  if(!is.null(start)){
    fit <- tryCatch(rssoptim(model = model, data = data, start = start, algo = algo), error = function(e) list(value = NA))
  } 
  
  if(!is.null(grid_start)) {
    if (grid_start != FALSE){
      fit <- grid_start_fit(model = model, data = data, n = grid_start, algo = algo, verb = verb)
    }#eo if (grid_start != FALSE)
  }#eo if(!is.null(grid_start))
  
  if(is.na(fit$value)){
    warning("The model could not be fitted :(\n")
    return(list(value = NA))
  }else{
    return(fit)
  }
}#eo get_fit

