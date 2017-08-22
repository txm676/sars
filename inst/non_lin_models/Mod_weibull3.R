#CUMULATIVE WEIBULL DISTRIBUTION with 3 parameters
model <- list(
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
