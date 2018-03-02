#Asymptotic Regression
model <- list(
  name=c("Asymptotic regression"),
  exp=expression(d - c*z^A),
  shape="convex",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R"),
  #initial values function
  init=function(data){#Ratkowsky 1983 p178
    #d determination (asymptote)
    d=max(data$S)+max(data$S)/4
    #Intermediate variable calculation
    Z=log(d-data$S)
    #we have also Z=log(c)+Xlog(z) -> linear regression
    dat=data.frame("a"=data$A,"Z"=Z)
    zf=stats::lm(Z~a,dat)$coefficients
    c(d,exp(zf[1]),exp(zf[2]))
  }
)
