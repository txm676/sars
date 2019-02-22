#Asymptotic Regression
model_asym <- function(data) {
  info <- sars::Table1[1L,]
  list(
  name = as.character(info$Model),
  formula = expression(S==d - c*z^A),
  exp = expression(d - c*z^A),
  shape = info$Model.shape,
  asymp = function(pars)pars["d"],
  parLim  =  c("Rplus","R","R"),
  #initial values function
  init = function(data){#Ratkowsky 1983 p178
    #d determination (asymptote)
    d = max(data$S) + .25*max(data$S)
    #Intermediate variable calculation
    Z = log(d-data$S)
    #we have also Z=log(c)+Xlog(z) -> linear regression
    dat = data.frame("a"=data$A, "Z"=Z)
    zf = stats::lm(Z~a,dat)$coefficients
    c(d, exp(zf[1L]), exp(zf[2L]))
  }
)
}
