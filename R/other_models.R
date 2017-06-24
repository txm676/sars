
## Other models (internal)

#Asymptotic Regression

asymp <- list(
  name=c("Asymptotic regression"),
  formula=expression(S == d - c*z^A),
  exp=expression(d - c*z^A),
  shape="convex",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R"),
  custStart=function(data)c(1,1,max(data$S)*2),
  #initial values function
  init=function(data){#Ratkowsky 1983 p178
    #d determination (asymptote)
    d=max(data$S)+max(data$S)/4
    #Intermediate variable calculation
    Z=log(d-data$S)
    #we have also Z=log(c)+Xlog(z) -> linear regression
    dat=data.frame("a"=data$A,"Z"=Z)
    zf=lm(Z~a,dat)$coefficients
    c(d,exp(zf[1]),exp(zf[2]))
  }
)


#RATIONAL FUNCITON ratkowski (1990)

ratio <- list(
  name=c("Rational function"),
  formula=expression(S == over( (c + z*A) , (1+d*A) ) ),
  exp=expression((c + z*A)/(1+d*A)),
  shape="convex",
  asymp=function(pars)pars["z"]/pars["d"],
  parLim = c("R","Rplus","unif"),
  custStart=function(data)c(1,1,0.000001),
  #initial values function
  init=function(data){c(0,0,.5)}
)


#gompertz model

gompertz <- list(
  name=c("Gompertz"),
  formula=expression(S==d*e^(-e^(-z*(A-c)))),
  exp=expression(d*exp(-exp(-z*(A-c)))),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R"),
  custStart=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    #d determination (asymptote)
    d<-max(data$S)+1
    #t.0 determination (min obs. age)
    t.0<-min(data$A)-1
    #Intermediate variable calculation
    Z=log(-log(data$S/d))
    #we have also Z=-kT + kt0 -> linear regression
    dat=data.frame("a"=data$A,"Z"=Z)
    reg=lm(Z~a,dat)$coefficients
    #transformations of coeficients
    k.first<--reg[2]
    k.second<-reg[1]/t.0
    k.final<-mean(c(k.first,k.second))
    #estimates return
    c(d,t.0,k.final)
  },
  #initials values function
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    #d determination (asymptote)
    d<-max(data$S)+1
    #t.0 determination (min obs. age)
    t.0<-min(data$A)-1
    #Intermediate variable calculation
    Z=log(-log(data$S/d))
    #we have also Z=-kT + kt0 -> linear regression
    dat=data.frame("a"=data$A,"Z"=Z)
    reg=lm(Z~a,dat)$coefficients
    #transformations of coeficients
    k.first<--reg[2]
    k.second<-reg[1]/t.0
    k.final<-mean(c(k.first,k.second))
    #estimates return
    c(d,t.0,k.final)
  }
)


#Beta-P function (cumulative)

betap <- list(
  name=c("Beta-P cumulative"),
  formula=expression(S == d*(1-(1+(A/c)^z)^-f) ),
  exp=expression(d*(1-(1+(A/c)^z)^-f)),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R","R"),
  custStart=function(data=data)c(max(data$S),10,.01,.5),
  #initial values function
  init=function(data){c(max(data$S)+1,.5,.5,.5)}
)


#LOGISTIC FUNCTION (HE & LEGENDRE 1996)

heleg <- list(
  name=c("Logistic (He & Legendre)"),
  formula=expression(S == over( c , (f + A^(-z)) ) ),
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
    zf=lm(Z~a,dat)$coefficients
    c(max(data$S)*exp(-zf[1]),exp(-zf[1]),-zf[2])
  }
)
