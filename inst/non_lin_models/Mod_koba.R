# Kobayashi logarithmic (KOBAYASHI 1975), convex upward, no asymptote
model <- list(
  name=c("Kobayashi"),
  formula=expression(s==c*log(1+ a/z)),
  exp=expression(c*log(1+ A/z)),
  shape="convex",
  asymp=function(pars)FALSE,
  #limits for parameters
  parLim = c("R","Rplus"),
  custStart=function(data)c(5,0.1),
  #initials values function
  init=function(data){
    c(max(data$data$S),30)
  }
)
