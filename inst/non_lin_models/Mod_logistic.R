#Standard Logistic function
#Formula given in Tjorve (2003)
#d confirmed as asymptote in Tjorve (2003)
#Starting parameter method taken from:
#https://bscheng.com/2014/05/07/modeling-logistic-growth-data-in-r/
model <- list(
  name=c("Logistic(Standard)"),
  formula=expression(S==d/(1 + exp(-z*A + c))),
  exp=expression(d/(1 + exp(-z*A + c))),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","Rplus","Rplus"),
  init=function(data){
    p <- data$S/(max(data$S) + 1)
    logitR <- log(p / (1-p)) #logit function (same as car::logit)
    cc <- coef(lm(logitR ~ data$A, data = data))
    d1 <- max(galap$s) + 1
    #use abs, because in the formulation they use in that website, the sign
    #of the c parameter is reversed (seemingly negative most of the time)
    c(d1, cc[2], abs(cc[1]))
  }
)