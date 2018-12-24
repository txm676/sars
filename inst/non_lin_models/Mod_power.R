# POWER MODEL (ARRHENIUS 1921)
model <- list(
    name = c("Power"),
    formula = expression(S == c * A ^ z),
    exp = expression(c * A ^ z),
    shape="convex",
    asymp=function(pars)FALSE,
    parLim = c("R", "R"),
    custStart=function(data)c(5,0.25),
    init = function(data){
      if (any(data$S == 0)){
        log.data = data.frame(A = log(data$A), S = log(data$S + .5))
      } else {
        log.data = log(data)
      }
      res = stats::lm(S ~ A, log.data)$coefficients
      res = c(exp(res[1]), res[2])
      names(res) = c("c", "z")
      return(res)
    }
)
