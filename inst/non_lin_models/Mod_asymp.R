#Asymptotic Regression
#Note that Tjorve and Williams et al give the function to the 
#power of -A. Checked and Flater (1996), https://www.jstor.org/stable/3001471?seq=1
#and the original Ratkowsky source give it as power of positive A so we have kept
#our use of this here.
model <- list(
  name = c("Asymptotic regression"),
  formula = expression(S == d - c*z^A),
  exp = expression(d - c*z^A),
  shape = "convex/sigmoid",
  asymp = function(pars)pars["d"],
  parLim  =  c("Rplus","R","Rplus"),
  #initial values function
  init = function(data){#Ratkowsky 1983 p178
    #d determination (asymptote)
    d = max(data$S)+max(data$S)/4
    #Intermediate variable calculation
    Z = log(d-data$S)
    #we have also Z=log(c)+Xlog(z) -> linear regression
    dat = data.frame("a"=data$A,"Z"=Z)
    zf = stats::lm(Z~a,dat)$coefficients
    c(d,exp(zf[1]),exp(zf[2]))
  }
)
