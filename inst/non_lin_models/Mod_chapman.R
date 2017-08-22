#Chapman–Richards 3 S = a [1 − exp(−bA)]c Flather (1996)
model <- list(
  name=c("Chapman Richards"),
  formula=expression(S == d * (1 - exp(-z*A)^c )),
  exp=expression(d * (1 - exp(-z*A)^c )),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","R","R"),
  custStart=function(data)c(10,.01,max(data$S)),
  #initials values function
  init=function(data){d=max(data$S); Z=( -log( (-data$S/(max(data$S)+1))+1))/data$A; z = mean(Z); c(d,z,1)}
)
