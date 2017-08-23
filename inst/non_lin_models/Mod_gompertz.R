#gompertz model
model <- list(
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
    reg=stats::lm(Z~a,dat)$coefficients
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
    reg=stats::lm(Z~a,dat)$coefficients
    #transformations of coeficients
    k.first<--reg[2]
    k.second<-reg[1]/t.0
    k.final<-mean(c(k.first,k.second))
    #estimates return
    c(d,t.0,k.final)
  }
)

