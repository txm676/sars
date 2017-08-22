# Kobayashi logarithmic (KOBAYASHI 1975), convex upward, no asymptote
model <- list(
  name=c("Kobayashi"),
  formula=expression(S==c*log(1+ A/z)),
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
