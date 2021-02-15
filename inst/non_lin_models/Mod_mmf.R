#"Morgan Mercier Family" curve (Williams et al. 2009 formula)
#have double checked and the Williams formula is definitely equivalent
#to the Tjorve and Godeau et al formulas.
#Found to be equivalent to mmf and deprecated
model <- list(
  name=c("MMF"),
  formula=expression(S==d/(1+c*A^(-z))),
  exp=expression(d/(1+c*A^(-z))),
  shape="convex/sigmoid",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","Rplus","Rplus"),
  custStart=function(data)c(max(data$S),5,.25),
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    d=(max(data$S)*4)
    newVar = log((d/data$S) - 1)
    reg = stats::lm(newVar~log(data$A))
    c=exp(reg$coefficients[1])
    z=-reg$coefficients[2]
    c(d,c,z)
  }
)
