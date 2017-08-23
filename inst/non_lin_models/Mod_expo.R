# EXPONENTIAL MODEL (GLEASON 1922)
model <- list(
    name=c("Exponential"),
    formula=expression(S==c+z*log(A)),
    exp=expression(c+z*log(A)),
    shape="convex",
    asymp=function(pars)FALSE,
    #limits for parameters
    parLim = c("R","R"),
    custStart=function(data)c(5,.25),
    #initials values function
    init=function(data){
      semilog.data = data.frame(log(data$A),data$S)
      names(semilog.data)=c("A","S")
      par=stats::lm(S~A,semilog.data)$coefficients
      names(par)=c("c","z")
      par
    }
)
