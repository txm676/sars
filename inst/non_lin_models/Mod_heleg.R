#LOGISTIC FUNCTION (HE & LEGENDRE 1996)
#doesn't make sense to fit negative relationships with this
#model, so pars should be Rplus.
#Tjorve's table in book chapter says asymptote par of Archibald logistic
#is the c parameter in ours, but have checked and we are correct with c/f.
#Table 3 of Tjorve (2009) states its shape is sigmoid but "only the convex
#part is used", and Williams (2009) discusses it being sigmoid but in logA
#space. We've found it can give both shapes in untransformed space, so
#have classified it as convex/sigmoid.
#Found to be equivalent to mmf and the latter deprecated.
#nb. Tjorve (2009) has a mistake in their formula for this model (Archibald
#logistic), and give corrected version in the sar book chapter
model <- list(
  name=c("Heleg(Logistic)"),
  formula=expression(S == c/(f + A^(-z))),
  exp=expression(c/(f + A^(-z))),
  shape="convex/sigmoid",
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
