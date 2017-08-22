# EXTENDED POWER MODEL 2 (TJORVE 2009)
model <- list(
  name=c("Extended Power model 2"),
  formula=expression(S==c*A^(z-(d/A))),
  exp=expression(c*A^(z-(d/A))),
  shape="sigmoid",
  asymp=function(pars)FALSE,
  custStart=function(data)c(5,.25,1),
  #limits for parameters
  parLim = c("Rplus","unif","R"),
  #initials values function
  init=function(data){return(c(0,0,0))}
)
