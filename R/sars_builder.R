## Try to get a function that build that creates sars.


sars_builder <- function() {

}

# if (!(is.matrix(data) | is.data.frame(data)))
# stop('data must be a matrix or dataframe')
# if (is.matrix(data)) data <- as.data.frame(data)
# if (anyNA(data)) stop('NAs present in data')
# data <- data[order(data[,1]),]
# colnames(data) <- c('A','S')
# #"Morgan Mercier Family" curve (Williams et al. 2009 formula)
# model <- list(
#   name=c("MMF"),
#   formula=expression(S==d/(1+c*A^(-z))),
#   exp=expression(d/(1+c*A^(-z))),
#   shape="sigmoid",
#   asymp=function(pars)pars["d"],
#   #limits for parameters
#   parLim = c("Rplus","Rplus","Rplus"),
#   custStart=function(data)c(max(data$S),5,.25),
#   init=function(data){
#     if(any(data$S==0)){data=data[data$S!=0,]}
#     d=(max(data$S)*4)
#     newVar = log((d/data$S) - 1)
#     reg = stats::lm(newVar~log(data$A))
#     c=exp(reg$coefficients[1])
#     z=-reg$coefficients[2]
#     c(d,c,z)
#   }
# )
#
# model <- compmod(model)
# fit <- get_fit(model = model, data = data, start = start,
# grid_start = grid_start, algo = 'Nelder-Mead',
#
# normaTest =  normaTest, homoTest = homoTest, verb = TRUE)
# if(is.na(fit$value)){
#   return(list(value = NA))
# }else{
#   obs <- obs_shape(fit)
#   fit$observed_shape <- obs$fitShape
#   fit$asymptote <- obs$asymp
#   class(fit) <- 'sars'
#   attr(fit, 'type') <- 'fit'
#   return(fit)
# }

check_data <- function(data) {
  if (!(is.matrix(data) | is.data.frame(data)))
    stop('data must be a matrix or dataframe')
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop('NAs present in data')
  data <- data[order(data[,1L]),]
  #
  colnames(data) <- c('A','S')
  data
}
