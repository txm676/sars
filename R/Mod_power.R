# POWER MODEL (ARRHENIUS 1921)

power <- list(
  name = c("Power"),
  formula = expression(s == c * a ^ z),
  exp = expression(c * A ^ z),
  shape = "convex",
  asymp = function(pars)FALSE,
  custStart = function(data) c(5, .25, 0),
  #limits for parameters
  parLim = c("R", "R"),
  #initials values function
  init=function(data){if(any(data$S==0)){log.data=data.frame(A=log(data$A),S=log(data$S+.5))}else{log.data=log(data)};res=lm(S~A,log.data)$coefficients;res = c(exp(res[1]),res[2]);names(res)=c("c","z");return(res)}
)
