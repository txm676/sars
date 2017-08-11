# POWER MODEL BIS (ROSENSWEIG 1995)

power_R <- list(
  name=c("Power_R"),
  formula=expression(s==f + c*a^z),
  exp=expression(f + c*A^z),
  shape="convex",
  asymp=function(pars)FALSE,
  custStart=function(data)c(5,.25,0),
  #limits for parameters
  parLim = c("R","R","R"),
  #initials values function
  init=function(data){if(any(data$s==0)){log.data=data.frame(a=log(data$a),s=log(data$s+.5))}else{log.data=log(data)};res=lm(s~a,log.data)$coefficients;res = c(exp(res[1]),res[2]);res = c(res,lm(s~a,data)$coefficients[1]);names(res)=c("f","c","z");return(res)}
)
