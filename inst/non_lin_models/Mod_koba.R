# Kobayashi logarithmic (KOBAYASHI 1975), convex upward, no asymptote
model <- list(
  name = c("Kobayashi"),
  formula = expression(S == c*log(1+ A/z)),
  exp = expression(c*log(1 + A/z)),
  shape = "convex",
  asymp = function(pars)FALSE,
  #limits for parameters
  parLim = c("R","Rplus"),
  custStart = function(data) c(5,0.1),
  #initials values function
  init = function(data){
    c(max(data$S), max(data$S))
  }
)



################## test init

#from mKobayashi 1975

# S = c if A = (e - 1) * z 
# with e = natural base log ~ 2.718

# p269
# c is the number of species occuring in it's caracteristic area (e - 1) * z




