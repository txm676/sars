
# Negative Exponential Family (internal)

#NEGATIVE EXPONENTIAL (Holdridge et al. 1971)

negexpo <- list(
  name=c("Negative exponential"),
  formula=expression(s == d*(1 - exp(-z*a) )),
  exp=expression(d*(1-exp(-z*A))),
  shape="convex",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","unif"),
  custStart=function(data)c(max(data$S),.01),
  #initials values function
  init=function(data){d=max(data$S); Z=( -log( (-data$S/(max(data$S)+1))+1))/data$A; z = mean(Z); c(d,z)}
)

#Chapman–Richards 3 S = a [1 − exp(−bA)]c Flather (1996)

chapman <- list(
  name=c("Chapman Richards"),
  formula=expression(S == d * (1 - exp(-z*A)^c )),
  exp=expression(d * (1 - exp(-z*A)^c )),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","R","R"),
  custStart=function(data)c(10,.01,max(data$S)),
  #initials values function
  init=function(data){d=max(data$S); Z=( -log( (-data$S/(max(data$S)+1))+1))/data$A; z = mean(Z); c(d,z,1)}
)

#CUMULATIVE WEIBULL DISTRIBUTION with 3 parameters

weibull3 <- list(
  name=c("Cumulative Weibull 3 par."),
  formula=expression(S == d(1 - exp(-c*A^z)) ),
  exp=expression(d*(1 - exp(-c*A^z)) ),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","Rplus","Rplus"),
  custStart=function(data)c(10,.01,max(data$S)),
  init=function(data){
    data = data[data$S!=0,]
    #c calculation (asymptote)
    d=max(data$S)+max(data$S)/4
    #z calculation
    Z=log(-log((d-data$S)/d)) # = log z + flogX
    Z[][Z == Inf]=NA
    c=exp(min(Z))
    dat=data.frame("A"=log(data$A),"S"=Z)
    c=exp(lm(S~A,dat)$coefficients[[1]])
    #f calculation
    z=lm(S~A,dat)$coefficients[[2]]
    c(d,c,z)
  }
)

#CUMULATIVE WEIBULL DISTRIBUTION with 4 parameters

weibull4 <- list(
  name=c("Cumulative Weibull 4 par."),
  formula=expression(S == d * (1 - exp(-c*A^z))^f ),
  exp=expression(d * (1 - exp(-c*A^z))^f ),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","Rplus","Rplus","Rplus"),
  custStart=function(data)c(10,.01,max(data$S),.01),
  #initial values function
  init=function(data){
    data = data[data$S!=0,]
    #c calculation (asymptote)
    d=max(data$S)+max(data$S)/4
    #z calculation
    Z=log(-log((d-data$S)/d)) # = log z + flogX
    Z[][Z == Inf]=NA
    c=exp(min(Z))
    dat=data.frame("A"=log(data[[1]]),"S"=Z)
    c=exp(lm(S~A,dat)$coefficients[[1]])
    #f calculation
    z=lm(S~A,dat)$coefficients[[2]]
    c(d,c,z,1)
  }
)
