#Beta-P function (cumulative)
model <- list(
  name=c("Beta-P cumulative"),
  exp=expression(d*(1-(1+(A/c)^z)^-f)),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R","R"),
  #initial values function
  init=function(data){c(max(data$S)+1,.5,.5,.5)}
)

