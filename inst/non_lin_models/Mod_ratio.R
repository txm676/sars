#RATIONAL FUNCITON ratkowski (1990)
model <- list(
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
