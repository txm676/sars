# POWER MODEL BIS (ROSENSWEIG 1995)
model <- list(
  name=c("PowerR"),
  formula=expression(S == f + c*A^z),
  exp=expression(f + c*A^z),
  shape="convex",
  asymp=function(pars)FALSE,
  custStart=function(data)c(5,.25,0),
  #limits for parameters
  parLim = c("R","R","R"),
  #initials values function
  init=function(data){if(any(data$S==0)){log.data=data.frame(A=log(data$A),S=log(data$S+.5))}else{log.data=log(data)};res=lm(S~A,log.data)$coefficients;res = c(exp(res[1]),res[2]);res = c(res,lm(S~A,data)$coefficients[1]);names(res)=c("f","c","z");return(res)}
)
