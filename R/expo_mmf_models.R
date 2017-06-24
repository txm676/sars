
### Exponential and mmf families (internal)


# EXPONENTIAL MODEL (GLEASON 1922)

expo <- list(
  name=c("Exponential"),
  formula=expression(s==c+z*log(A)),
  exp=expression(c+z*log(A)),
  shape="convex",
  asymp=function(pars)FALSE,
  #limits for parameters
  parLim = c("R","R"),
  custStart=function(data)c(5,.25),
  #initials values function
  init=function(data){
    semilog.data = data.frame(log(data$A),data$S)
    names(semilog.data)=c("A","S")
    par=lm(S~A,semilog.data)$coefficients
    names(par)=c("c","z")
    par
  }
)


# Kobayashi logarithmic (KOBAYASHI 1975), convex upward, no asymptote

koba <- list(
  name=c("Kobayashi"),
  formula=expression(s==c*log(1+ a/z)),
  exp=expression(c*log(1+ A/z)),
  shape="convex",
  asymp=function(pars)FALSE,
  #limits for parameters
  parLim = c("R","Rplus"),
  custStart=function(data)c(5,0.1),
  #initials values function
  init=function(data){
    c(max(data$data$S),30)
  }
)


#"Morgan Mercier Family" curve (Williams et al. 2009 formula)

mmf <- list(
  name=c("MMF"),
  formula=expression(s==over(d,1+c*a^(-z))),
  exp=expression(d/(1+c*A^(-z))),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","Rplus","Rplus"),
  custStart=function(data)c(max(data$S),5,.25),
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    d=as.real(max(data$S)*4)
    newVar = log((d/data$S) - 1)
    reg = lm(newVar~log(data$A))
    c=exp(reg$coefficients[1])
    z=-reg$coefficients[2]
    c(d,c,z)
  }
)


# MONOD CURVE (MONOD 1950, Willimas et al. 2009 formula)

monod <- list(
  name=c("monod"),
  formula=expression(s==over(d,1+c*a^(-1))),
  exp=expression(d/(1+c*A^(-1))),
  shape="convex",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","Rplus"),
  custStart=function(data)c(quantile(data$A,c(0.25)),max(data$S)),
  #initials values function
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    d=as.real(max(data$A)+max(data$S)/4)
    c=data[[1]]*(d/data$S - 1)
    c(d,quantile(c,c(0.25)))
  }
)
