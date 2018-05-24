#LOGISTIC FUNCTION (HE & LEGENDRE 1996)
model <- list(
  name=c("Heleg(Logistic)"),
  formula=expression(S == c/(f + A^(-z))),
  exp=expression(c/(f + A^(-z))),
  shape="sigmoid",
  asymp=function(pars)pars["c"]/pars["f"],
  parLim = c("Rplus","Rplus","Rplus"),
  custStart=function(data)c(max(data$S),10,.01),
  #initial values function
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    #c calculation (asymptote)
    c=max(data$S)+1
    #Intermediate variable calculation
    #long=length(data[[2]])
    Z=log((c/data$S) -1)
    #We have the z and f init values by linear regression of Z on data[[1]]
    dat=data.frame("a"=data$A,"Z"=Z)
    zf=stats::lm(Z~a,dat)$coefficients
    c(max(data$S)*exp(-zf[1]),exp(-zf[1]),-zf[2])
  }
)
