# EXTENDED POWER MODEL 1 (TJORVE 2009)
model <- list(
  name=c("Extended Power model 1"),
  formula=expression(S==c*A^(z*A^-d)),
  exp=expression(c*A^(z*A^-d)),
  shape="sigmoid",
  asymp=function(pars)FALSE,
  parLim = c("Rplus","R","R"),
  custStart=function(data)c(5,.25,.15),
  #initials values function
  init=function(data){if(any(data$S==0)){log.data=data.frame(S=log(data$a),S=log(data$S+.5))}else{log.data=log(data)};res=lm(S~A,log.data)$coefficients;res=c(exp(res[1]),res[2],.15);names(res)=model$parNames;return(res)}
)

