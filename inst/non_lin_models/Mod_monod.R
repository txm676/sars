# MONOD CURVE (MONOD 1950, Willimas et al. 2009 formula)
model <- list(
  name=c("Monod"),
  formula=expression(S==d/(1+c*A^(-1))),
  exp=expression(d/(1+c*A^(-1))),
  shape="convex",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","Rplus"),
  custStart=function(data)c(stats::quantile(data$A,c(0.25)),max(data$S)),
  #initials values function
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    d=as.double(max(data$A)+max(data$S)/4)
    c=data[[1]]*(d/data$S - 1)
    c(d,stats::quantile(c,c(0.25)))
  }
)
