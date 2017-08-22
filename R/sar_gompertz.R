#' Fit the Asymptotic regression model

#' @description Fit the Asymptotic regression model to SAR data
#' @usage sar_gompertz(data, custstart = NULL, normtest = 'lillie')
#' @param data A dataset in the form of a dataframe with two columns: 
#'   the first with island/site areas, and the second with the species richness
#'   of each island/site.
#' @return 
#' @examples
#' data(galap)
#' fit <- sar_gompertz(galap)
#' summary(fit)
#' plot(fit)
#' @export

sar_gompertz <- function(data=galap, start = NULL){
if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
if (is.matrix(data)) data <- as.data.frame(data) 
if (anyNA(data)) stop('NAs present in data') 
data <- data[order(data[,1]),] 
colnames(data) <- c('A','S') 
#gompertz model
model <- list(
  name=c("Gompertz"),
  formula=expression(S==d*e^(-e^(-z*(A-c)))),
  exp=expression(d*exp(-exp(-z*(A-c)))),
  shape="sigmoid",
  asymp=function(pars)pars["d"],
  parLim = c("Rplus","R","R"),
  custStart=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    #d determination (asymptote)
    d<-max(data$S)+1
    #t.0 determination (min obs. age)
    t.0<-min(data$A)-1
    #Intermediate variable calculation
    Z=log(-log(data$S/d))
    #we have also Z=-kT + kt0 -> linear regression
    dat=data.frame("a"=data$A,"Z"=Z)
    reg=lm(Z~a,dat)$coefficients
    #transformations of coeficients
    k.first<--reg[2]
    k.second<-reg[1]/t.0
    k.final<-mean(c(k.first,k.second))
    #estimates return
    c(d,t.0,k.final)
  },
  #initials values function
  init=function(data){
    if(any(data$S==0)){data=data[data$S!=0,]}
    #d determination (asymptote)
    d<-max(data$S)+1
    #t.0 determination (min obs. age)
    t.0<-min(data$A)-1
    #Intermediate variable calculation
    Z=log(-log(data$S/d))
    #we have also Z=-kT + kt0 -> linear regression
    dat=data.frame("a"=data$A,"Z"=Z)
    reg=lm(Z~a,dat)$coefficients
    #transformations of coeficients
    k.first<--reg[2]
    k.second<-reg[1]/t.0
    k.final<-mean(c(k.first,k.second))
    #estimates return
    c(d,t.0,k.final)
  }
)


model <- compmod(model) 
fit <- rssoptim(model = model, data = data, custstart = start, algo = 'Nelder-Mead') 
obs <- obs_shape(fit) 
fit$observed_shape <- obs$fitShape 
fit$asymptote <- obs$asymp 
class(fit) <- 'sars' 
attr(fit, 'type') <- 'fit' 
return(fit) 
}#end of sar_gompertz
