# PERSISTENCE FUNCTION 1 (TJORVE 2009)
model <- list(
  name=c("Persistence function 1"),
  formula=expression(s == c*a^z * exp(-d*a)),
  exp=expression(c*A^z * exp(-d*A)),
  mod=s~c*a^z * exp(-d*a),
  shape="convex",
  asymp=function(pars)FALSE,
  custStart=function(data)c(5,.25,.15),
  #limits for parameters
  parLim = c("Rplus","Rplus","Rplus"),
  #initials values function
  init=function(data){if(any(data$S==0)){log.data=data.frame(A=log(data$A),S=log(data$S+.5))}else{log.data=log(data)};res=lm(S~A,log.data)$coefficients;res=c(res[1],res[2],1);names(res)=P1$paramnames;return(res)}
)
