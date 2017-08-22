#CUMULATIVE WEIBULL DISTRIBUTION with 4 parameters
model <- list(
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
