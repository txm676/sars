#Standard Logistic function
#Formula given in Tjorve (2003)
#d confirmed as asymptote in Tjorve (2003)
#Starting parameter method taken from:
#https://bscheng.com/2014/05/07/modeling-logistic-growth-data-in-r/
#Shape given as sigmoid, but it can have observed shapes of linear, convex up
#and convex down. 
#Pars all set to Rplus, similar to heleg and mmf. d has to be as is the 
#asymptote, and have tested the other two - z should be Rplus but with c it is
#trickier. For some datasets, setting c to R results in a better fit, but then 
#this also means for other datasets it converges on a worse fit. Having it as
#Rplus seems to allow it to better fit the "typical" sigmoid shape so have left
#it as that for now.
model <- list(
  name=c("Logistic(Standard)"),
  formula=expression(S==d/(1 + exp(-z*A + c))),
  exp=expression(d/(1 + exp(-z*A + c))),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  #limits for parameters
  parLim = c("Rplus","Rplus","Rplus"),
  init=function(data){
    if (any(data$S == 0)){
      p <- (data$S + 0.1)/(max(data$S) + 1)
    } else{
      p <- data$S/(max(data$S) + 1)
    }
    logitR <- log(p / (1 - p)) #logit function (same as car::logit)
    cc <- coef(lm(logitR ~ data$A))
    d1 <- max(data$S) + 1
    #use abs, because in the formulation they use in that website, the sign
    #of the c parameter is reversed (seemingly negative most of the time)
    c(d1, cc[2], abs(cc[1]))
  }
)