########################################################################
###########Dengler Functions################################
############################################################

library(sars)

dat <- ro()
colnames(dat) <- c("A", "S")

pp <- sar_power(dat)

lin_pow(dat)

A = 5
B = 3
z1 = 0.2
z2 = 0.3
c = 40
k = 1

S =	10 ^ (log10(c) + (log10(A) < log10(B)) * (z1 * log10(A)) + 
            (log10(A) >= log10(B)) * ((z1 * (log10(B) + z2)) * (log10(A) - log10(B))))

#get dengler to independently calculate an S value for each function using set A,B, z and c pars

res1 <- tryCatch(optim(startMod, rssfun, hessian = TRUE, data = data, method = algo, control = list(maxit = 50000) ),
                 error = function(e){e})


#starting values

#log10(B) = 5.5 (i.e. about 300000 km²)
#z1 about 0.2
#z2 about 0.7
#log10(c) about 2.5
#k you can set to 1


nls(formula = S ~ 10 ^ (log10(c) +(log10(A) < log10(B)) * (z1 * log10(A)) + 
                          (log10(A) >= log10(B)) * ((z1 * (log10(B) + z2)) * (log10(A) - log10(B)))), 
    data = dat, start = data.frame(c = 100, z1 = 0.2, z2 = 0.7, B = 3000), lower=c(0.1,0.01, 0.01, 1),
    algorithm = "port",  control = list(maxiter = 5000))



#log version
log_S =	log_c+(log_A<log_T)*(z1*log_A)+(log_A>=log_T)*(z1*log_T+z2*(log_A-log_T))

LS = log10(c) + (log10(A) < log10(B)) + (z1 * log10(A)) + 
            (log10(A) >= log10(B)) * ((z1 * (log10(B) + z2)) * (log10(A) - log10(B)))



###smooth (k = 1)

S = 10 ^ (log10(c) + (z2 - z1) * (log(exp(1 * log10(B) - log10(A)) + 1) + log10(A)) +
            z1 * log10(A) - (z2 - z1) * (log(exp(log10(B)) + 1)))
#log version

S = log10(c) + (z2 - z1) * (log(exp(1 * log10(B) - log10(A)) + 1) + log10(A)) +
  z1 * log10(A) - (z2 - z1) * (log(exp(log10(B)) + 1))


###smooth (k varies)

S = 10 ^ (log10(c) + (z2 - z1) * (log(exp(k * log10(B) - k * log10(A)) + 1) / k + log10(A)) +
                                    z1 * log10(A) - (z2 - z1) * (log(exp(k * log10(B)) + 1) /k))


#log version

S = log10(c) + (z2 - z1) * (log(exp(k * log10(B) - k * log10(A)) + 1) / k + log10(A)) +
  z1 * log10(A) - (z2 - z1) * (log(exp(k * log10(B)) + 1) /k)


#uses a brute force approach; if number of combinations is very large this can take a long time

#if given, pars should be a data.frame where each column is a 
#parameter and each row is a set of starting values. Similar to 
#what is returned by expand.grid function. The order of columns shoud be: c, z1, z2, B.

grid_start_deng <- function(dat, pars = NULL, mod = "break"){
  
  if (!(is.data.frame(dat) || is.matrix(dat))) stop("dat must be a dataframe or matrix")
  if (is.matrix(dat)) dat <- as.data.frame(dat)
  if (!all(colnames(dat) == c("A", "S"))){
    colnames(dat) <- c("A", "S")
    message("\n", "Renaming columns as 'A' and 'S', in that order", "\n")
  }
  
  #Dengler's log zero approach
  if (any(dat$S == 0)){
    wh <- which(dat$S == 0)
    dat$S[wh] <- 0.25
    message("\n", "Zeros found in richness values: these have been converted to 0.25", "\n")
  }
  
  if (! mod %in% c("break", "break_log", "smooth4", "smooth4_log", "smooth5", "smooth5_log")) stop("mod not recognised")
  if (is.null(pars)){ #if pars not provided create a df of different starting values
  c <- seq(0.01, 500, 50)
  z1 <- seq(0.01, 1, 0.1)
  z2 <- seq(0.01, 1, 0.1)
  #get range of values for breakpoints
  dr2 <- max(dat$A)[1]; dr1 <- min(dat$A)[1]
  B <- exp(seq(log(dr1), log(dr2))) 
  
  grid_vals <- expand.grid(c, z1, z2, B)
  colnames(grid_vals) <- c("c", "z1", "z2", "B")
  #if smooth 5 parameter need to add starting values for k
  if (mod %in% c("smooth5", "smooth5_log")){
    grid_vals$k <- 1
  }
  } else { #else do some checks and use the user's starting values df
    if (!is.data.frame(pars)) stop("pars argument should be a dataframe or NULL")
    if (mod %in% c("smooth5", "smooth5_log")){
      if (ncol(pars) != 5) stop("not enough columns in pars")
    } else {
      if (ncol(pars) != 4) stop("not enough columns in pars")
    }
    grid_vals <- pars
  }#eo is.null(pars)
  
  message("\n", paste(nrow(grid_vals), "combinations of starting parameter values used"), "\n", "\n")
  
  if (mod == "break"){
    reg <- apply(grid_vals, 1, function(x){

    tryCatch(nls(formula = S ~ 10 ^ (log10(c) +(log10(A) < log10(B)) * (z1 * log10(A)) + 
                            (log10(A) >= log10(B)) * ((z1 * (log10(B) + z2)) * (log10(A) - log10(B)))), 
      data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4]), lower=c(0.01, 0.01, 0.01, 1),
      algorithm = "port"), error = function(e) list(value = NA))
    })#eo apply
  } else if (mod == "break_log"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = log10(S) ~ log10(c) + (log10(A) < log10(B)) + (z1 * log10(A)) + 
                     (log10(A) >= log10(B)) * ((z1 * (log10(B) + z2)) * (log10(A) - log10(B))), 
                   data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4]), lower=c(0.01, 0.01, 0.01, 1),
                   algorithm = "port"), error = function(e) list(value = NA))
    })#eo a
  } else if (mod == "smooth4"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = S ~ 10 ^ (log10(c) + (z2 - z1) * (log(exp(1 * log10(B) - log10(A)) + 1) + log10(A)) +
                                     z1 * log10(A) - (z2 - z1) * (log(exp(log10(B)) + 1))), 
                   data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4]), lower=c(0.01, 0.01, 0.01, 1),
                   algorithm = "port"), error = function(e) list(value = NA))
    })#eo a
  } else if (mod == "smooth4_log"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = log10(S) ~ log10(c) + (z2 - z1) * (log(exp(1 * log10(B) - log10(A)) + 1) + log10(A)) +
                     z1 * log10(A) - (z2 - z1) * (log(exp(log10(B)) + 1)), 
                   data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4]), lower=c(0.01, 0.01, 0.01, 1),
                   algorithm = "port"), error = function(e) list(value = NA))
    })#eo a
  } else if (mod == "smooth5"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = S ~ 10 ^ (log10(c) + (z2 - z1) * (log(exp(k * log10(B) - k * log10(A)) + 1) / k + log10(A)) +
                                                z1 * log10(A) - (z2 - z1) * (log(exp(k * log10(B)) + 1) /k)), 
                   data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4], k = x[5]), 
                   lower=c(0.01, 0.01, 0.01, 1, 0.01),
                   algorithm = "port"), error = function(e) list(value = NA))
    })#eo a
  } else if (mod == "smooth5_log"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = log10(S) ~ log10(c) + (z2 - z1) * (log(exp(k * log10(B) - k * log10(A)) + 1) / k + log10(A)) +
                     z1 * log10(A) - (z2 - z1) * (log(exp(k * log10(B)) + 1) /k), 
                   data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4], k = x[5]), 
                   lower=c(0.01, 0.01, 0.01, 1, 0.01),
                   algorithm = "port"), error = function(e) list(value = NA))
    })#eo a
  }
  reg2 <- reg[sapply(reg, function(x) all(!is.na(x)))] #remove NA elements (i.e. subset a list)
  
  if (length(reg2) == 0) stop("No model fits converged: try providing different starting values")
  
  if (length(reg2) == 1){
    finMod <- reg2[[1]]
  } else {
    regRSS <- vapply(reg2, function(x) x$m$deviance(), FUN.VALUE = numeric(1))
    rssMin <- which.min(regRSS)[1]
    finMod <- reg2[[rssMin]]
  }#eo if
  return(finMod)
}


gg <- grid_start_deng(dat)

gg2 <- grid_start_deng(dat, mod = "break_log")

gg3 <- grid_start_deng(dat, mod = "smooth4")

gg4 <- grid_start_deng(dat, mod = "smooth4_log")

gg5 <- grid_start_deng(dat, mod = "smooth5")

gg6 <- grid_start_deng(dat, mod = "smooth5_log")

######################################
#test plots for Jurgen
#########################################
jpeg("Breakpoint_test.jpeg", width = 25, height = 37, res = 600, units = "cm")

par(mfrow = c(3, 2))
par("mfg")

plot(dat$A, dat$S, col = "black", pch = 16)
points(dat$A, gg$m$fitted(), col = "red")

par("mfg")
plot(log10(dat$A), log10(dat$S), col = "black", pch = 16)
points(log10(dat$A), gg2$m$fitted(), col = "red")

par("mfg")
plot(dat$A, dat$S, col = "black", pch = 16)
points(dat$A, gg3$m$fitted(), col = "red")

par("mfg")
plot(log10(dat$A), log10(dat$S), col = "black", pch = 16)
points(log10(dat$A), gg4$m$fitted(), col = "red")


par("mfg")
plot(dat$A, dat$S, col = "black", pch = 16)
points(dat$A, gg5$m$fitted(), col = "red")

par("mfg")
plot(log10(dat$A), log10(dat$S), col = "black", pch = 16)
points(log10(dat$A), gg6$m$fitted(), col = "red")


dev.off()





deng_NB <- function(data, start = NULL, grid_start = NULL, normaTest =  "lillie",
                    homoTest = "cor.fitted"){
  if (!(is.matrix(data) || is.data.frame(data))) stop('data must be a matrix or dataframe') 
  if (is.matrix(data)) data <- as.data.frame(data) 
  if (anyNA(data)) stop('NAs present in data') 
  data <- data[order(data[,1]),] 
  colnames(data) <- c('A','S') 
  #"Morgan Mercier Family" curve (Williams et al. 2009 formula)
  model <- list(
    name=c("Normal_BP"),
    formula=expression(S==10 ^ (log10(c) +(log10(A) < log10(B)) * (z1 * log10(A)) + 
                                  (log10(A) >= log10(B)) * ((z1 * (log10(B) + z2)) * (log10(A) - log10(B))))),
    exp=expression(10 ^ (log10(c) +(log10(A) < log10(B)) * (z1 * log10(A)) + 
                           (log10(A) >= log10(B)) * ((z1 * (log10(B) + z2)) * (log10(A) - log10(B))))),
    shape="breakpoint",
    asymp=function(pars)FALSE,
    
    ##########################################################
    
    
    #limits for parameters
    parLim = c("Rplus","Rplus","Rplus"),
    custStart=function(data)c(max(data$S),5,.25),
    init=function(data){
      if(any(data$S==0)){data=data[data$S!=0,]}
      d=(max(data$S)*4)
      newVar = log((d/data$S) - 1)
      reg = stats::lm(newVar~log(data$A))
      c=exp(reg$coefficients[1])
      z=-reg$coefficients[2]
      c(d,c,z)
    }
  )
  
  model <- compmod(model) 
  fit <- get_fit(model = model, data = data, start = start, grid_start = grid_start, algo = 'Nelder-Mead', 
                 normaTest =  normaTest, homoTest = homoTest, verb = TRUE) 
  if(is.na(fit$value)){ 
    return(list(value = NA)) 
  }else{ 
    obs <- obs_shape(fit) 
    fit$observed_shape <- obs$fitShape 
    fit$asymptote <- obs$asymp 
    class(fit) <- 'sars' 
    attr(fit, 'type') <- 'fit' 
    return(fit) 
  } 
}#end of sar_mmf






