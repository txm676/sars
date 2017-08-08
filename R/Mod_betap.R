
#Beta-P function (cumulative)

betap <- list(
  name=c("Beta-P cumulative"),
  formula=expression(S == d*(1-(1+(A/c)^z)^-f) ),
  exp=expression(d*(1-(1+(A/c)^z)^-f)),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R","R"),
  custStart=function(data=data)c(max(data$S),10,.01,.5),
  #initial values function
  init=function(data){c(max(data$S)+1,.5,.5,.5)}
)

