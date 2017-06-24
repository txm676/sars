

###Power family models (internal)

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


# EXTENDED POWER MODEL 1 (TJORVE 2009)

epm1 <- list(
  name=c("Extended Power model 1"),
  formula=expression(s==c*a^(z*a^-d)),
  exp=expression(c*A^(z*A^-d)),
  shape="sigmoid",
  asymp=function(pars)FALSE,
  parLim = c("Rplus","R","R"),
  custStart=function(data)c(5,.25,.15),
  #initials values function
  init=function(data){if(any(data$S==0)){log.data=data.frame(S=log(data$a),S=log(data$S+.5))}else{log.data=log(data)};res=lm(S~A,log.data)$coefficients;res=c(exp(res[1]),res[2],.15);names(res)=model$parNames;return(res)}
)


# EXTENDED POWER MODEL 2 (TJORVE 2009)

epm2 <- list(
  name=c("Extended Power model 2"),
  formula=expression(s==c*a^(z-(d/a))),
  exp=expression(c*A^(z-(d/A))),
  shape="sigmoid",
  asymp=function(pars)FALSE,
  custStart=function(data)c(5,.25,1),
  #limits for parameters
  parLim = c("Rplus","unif","R"),
  #initials values function
  init=function(data){return(c(0,0,0))}
)


# PERSISTENCE FUNCTION 1 (TJORVE 2009)

P1 <- list(
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

# PERSISTENCE FUNCTION 2 (TJORVE 2009)

P2 <- list(
  name=c("Persistence function 2"),
  formula=expression(s == c*A^z * exp(-d/A)),
  exp=expression(c*A^z * exp(-d/A)),
  mod=s~c*a^z * exp(-d/a),
  shape="sigmoid",
  asymp=function(pars)FALSE,
  custStart=function(data)c(5,.25,.15),
  #limits for parameters
  parLim = c("Rplus","Rplus","R"),
  #initials values function
  init=function(data){if(any(data$S==0)){log.data=data.frame(A=log(data$A),S=log(data$S+.5))}else{log.data=log(data)};res=lm(S~A,log.data)$coefficients;c(res[1],res[2],.1)}
)
