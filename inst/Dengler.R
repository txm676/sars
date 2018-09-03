########################################################################
###########Dengler Functions################################
############################################################

library(sars)
library(foreach)
library(doParallel)
library(cluster)
library(dplyr)


dat <- read.csv("D:\\documents\\Work\\On-going projects\\Dengler\\test.csv")
colnames(dat) <- c("A", "S")

#####################################
###Second test dataset#############
################################

dat <- read.csv("D:\\documents\\Work\\On-going projects\\Dengler\\test2.csv")#test2
un <- unique(dat$group)
l2 <- lapply(un, function(x){
  f <- dplyr::filter(dat, group == x)
  f <- dplyr::select(f, c(A, S))
  f
})

#NB H has some NAs in area column so removed (check with Jurgen)

#uses a brute force approach; if number of combinations is very large this can take a long time

#if given, pars should be a data.frame where each column is a 
#parameter and each row is a set of starting values. Similar to 
#what is returned by expand.grid function. The order of columns shoud be: c, z1, z2, B.

#nls fit objects bring environment with them so can be very large (maybe just return best startign values)

grid_start_deng <- function(dat, pars = NULL, mod = "break", extensive = FALSE){
  
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
  
  #get range of values for breakpoints
  dr2 <- max(dat$A)[1]; dr1 <- min(dat$A)[1]  
  if (!extensive){
    c <- c(1, 50, 100, 500, 2000)
    z1 <- seq(0.01, 1, 0.15)
    z2 <- seq(0.01, 1, 0.15)
    B <- exp(seq(log(dr1), log(dr2))) 
  } else {
    c <- c(1, 5, 20, 50, 100, 250, 500, 1000, 3000)
    z1 <- seq(0.01, 1, 0.1)
    z2 <- seq(0.01, 1, 0.1)
    B <- c(1, 10, 100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 4000000, 12000000)
  }

  if (mod %in% c("smooth5", "smooth5_log")){ #if smooth 5 parameter need to add starting values for k
    if (!extensive){
      k <- c(1, dr2 * 0.005, dr2 * 0.25)
    } else {
      k <- c(0.5, 1, 3, 5, 10, 25, 50, 100, 500, 1000)
    }
    grid_vals <- expand.grid(c, z1, z2, B, k)
    colnames(grid_vals) <- c("c", "z1", "z2", "B", "k")
  } else {
    grid_vals <- expand.grid(c, z1, z2, B)
    colnames(grid_vals) <- c("c", "z1", "z2", "B")
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
  
  message("\n", paste(nrow(grid_vals), "combinations of starting parameter values used"), "\n")
  
  if (mod == "break"){
    reg <- apply(grid_vals, 1, function(x){

    tryCatch(nls(formula = S ~ 10 ^ (log10(c) +(log10(A) < log10(B)) * (z1 * log10(A)) + 
                            (log10(A) >= log10(B)) * ((z1 * (log10(B) + z2)) * (log10(A) - log10(B)))), 
      data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4]), lower=c(0.01, 0, 0, 1),
      upper = c(max(dat$S), 1, 1, dr2), algorithm = "port"), error = function(e) list(value = NA))
    })#eo apply
  } else if (mod == "break_log"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = log10(S) ~ log10(c) + (log10(A) < log10(B)) * (z1 * log10(A)) + 
                     (log10(A) >= log10(B)) * (z1 * log10(B) + z2 * (log10(A) - log10(B))), 
                   data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4]), lower=c(0.01, 0, 0, 1),
                   upper = c(max(dat$S), 1, 1, dr2), algorithm = "port"), error = function(e) list(value = NA))
    })#eo a
  } else if (mod == "smooth4"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = S ~ 10 ^ (log10(c) + (z2 - z1) * (log(exp(1 * log10(B) - log10(A)) + 1) + log10(A)) +
                                     z1 * log10(A) - (z2 - z1) * (log(exp(log10(B)) + 1))), 
                   data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4]), lower=c(0.01, 0, 0, 1),
                   upper = c(max(dat$S), 1, 1, dr2), algorithm = "port"), error = function(e) list(value = NA))
    })#eo a
  } else if (mod == "smooth4_log"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = log10(S) ~ log10(c) + (z2 - z1) * (log(exp(1 * log10(B) - log10(A)) + 1) + log10(A)) +
                     z1 * log10(A) - (z2 - z1) * (log(exp(log10(B)) + 1)), 
                   data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4]), lower=c(0.01, 0, 0, 1),
                   upper = c(max(dat$S), 1, 1, dr2), algorithm = "port"), error = function(e) list(value = NA))
    })#eo a
  } else if (mod == "smooth5"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = S ~ 10 ^ (log10(c) + (z2 - z1) * (log(exp(k * log10(B) - k * log10(A)) + 1) / k + log10(A)) +
                                                z1 * log10(A) - (z2 - z1) * (log(exp(k * log10(B)) + 1) /k)), 
                   data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4], k = x[5]), 
                   lower=c(0.01, 0, 0, 1, 0.01),
                   upper = c(max(dat$S), 1, 1, dr2), algorithm = "port"), error = function(e) list(value = NA))
    })#eo a
  } else if (mod == "smooth5_log"){
    reg <- apply(grid_vals, 1, function(x){
      
      tryCatch(nls(formula = log10(S) ~ log10(c) + (z2 - z1) * (log(exp(k * log10(B) - k * log10(A)) + 1) / k + log10(A)) +
                   z1 * log10(A) - (z2 - z1) * (log(exp(k * log10(B)) + 1) /k), 
                  data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4], k = x[5]), 
                  lower=c(0.01, 0, 0, 1, 0.01),
                  upper = c(max(dat$S), 1, 1, dr2), algorithm = "port"), error = function(e) list(value = NA))

    })#eo a
  }
  reg2 <- reg[sapply(reg, function(x) all(!is.na(x)))] #remove NA elements (i.e. subset a list)
  
  if (length(reg2) == 0 && extensive == FALSE){
    message("\n","No models converged: re-starting with more starting parameter value combinations","\n")
    newA <- tryCatch(grid_start_deng(dat, pars = NULL, mod = mod, extensive = TRUE), error = function(e) list(value = NA))
    if (length(newA) == 1) stop("No model fits converged: try providing different starting values")
    return(newA)
  }
  
  if (length(reg2) == 0) stop("No models converged: try providing different starting values")
  
  if (length(reg2) == 1){
    finMod <- reg2[[1]]
  } else {
    regRSS <- vapply(reg2, function(x) x$m$deviance(), FUN.VALUE = numeric(1))
    rssMin <- which.min(regRSS)[1]
    finMod <- reg2[[rssMin]]
  }#eo if
  return(finMod)
}


##################################
##power model results
###############################
pow <- function(dat){
  colnames(dat) <- c("A", "S")
  #Dengler's log zero approach
  if (any(dat$S == 0)){
    wh <- which(dat$S == 0)
    dat$S[wh] <- 0.25
    message("\n", "Zeros found in richness values: these have been converted to 0.25", "\n")
  }
  
  nlp <- sar_power(dat)
  
  nlp2 <- tryCatch(nls(formula = S ~ c*A^z, 
               data = dat, start = data.frame(c = nlp$par[1], z =  nlp$par[2])), error = function(e) list(value = NA))
  
  
  log.data =  log.data = data.frame(A = log10(dat$A), S = log10(dat$S))
  lp <- lm(S ~ A, data = log.data)
  
  nlpSum <- c(AICcmodavg::AICc(nlp2), nlp2$m$getPars()[1], nlp2$m$getPars()[2])
  lpSum <- c(AICcmodavg::AICc(lp), 10 ^ lp$coefficients[1], lp$coefficients[2])
  res <- c(nlpSum, lpSum)
  return(res)
}


cores = 9
cl = makeCluster(cores); on.exit(stopCluster(cl))
registerDoParallel(cl)
i = 1 #Dummy line for RStudio warnings


modz <- c("break", "break_log", "smooth4", "smooth4_log", "smooth5", "smooth5_log")

dat <- l2[[5]] 

allList2 = foreach(i=seq(from=1, to=length(modz), by=1))  %dopar% { 
gg2 <- tryCatch(grid_start_deng(dat, mod = modz[i], extensive = FALSE), error = function(e) list(value = NA))
gg2
}



##extract the best models from the ll list object
gg <- allList[[1]]
gg2 <- allList[[2]]
gg3 <- allList[[3]]
gg4 <- allList[[4]]
gg5 <- allList[[5]]
gg6 <- allList[[6]]



x <- c(0.6, 0.3, 0.01, 798600, 1)

nls(formula = S ~ 10 ^ (log10(c) + (z2 - z1) * (log(exp(k * log10(B) - k * log10(A)) + 1) / k + log10(A)) +
                          z1 * log10(A) - (z2 - z1) * (log(exp(k * log10(B)) + 1) /k)), 
    data = dat, start = data.frame(c = x[1], z1 = x[2], z2 = x[3], B = x[4], k = x[5]), 
    lower=c(0.01, 0, 0, 1, 0.01),
    upper = c(max(dat$S), 1, 1, dr2, 10000), control = list(maxiter = 5000), algorithm = "port", trace = TRUE)


allList = foreach(i=seq(from=1, to=length(l2), by=1))  %dopar% { 
  
  dat <- l2[[i]]
  
  gg <- tryCatch(grid_start_deng(dat), error = function(e) list(value = NA))
  
  gg2 <- tryCatch(grid_start_deng(dat, mod = "break_log"), error = function(e) list(value = NA))
  
  gg3 <- tryCatch(grid_start_deng(dat, mod = "smooth4"), error = function(e) list(value = NA))
  
  gg4 <- tryCatch(grid_start_deng(dat, mod = "smooth4_log"), error = function(e) list(value = NA))
  
  gg5 <- tryCatch(grid_start_deng(dat, mod = "smooth5"), error = function(e) list(value = NA))
  
  gg6 <- tryCatch(grid_start_deng(dat, mod = "smooth5_log", extensive = T), error = function(e) list(value = NA))
  
  ll <- list(gg, gg2, gg3, gg4, gg5, gg6)
  ll
}

save(allList, file = "allList.R")



replicate(8, tryCatch(grid_start_deng(dat, mod = "smooth5_log", extensive = TRUE), error = function(e) list(value = NA)))


g##############################################
####extract each element of allList and plot#####
######################################################

#check if any are NAs
sapply(allList, function(x) sapply(x, length))

#get AICc and model parameter values
mod_sum <- function(x){
 # if (class(x) == "sars"){
  #  ac <- summary(x)$AICc
   # pars <- x$par
 # } else {
    ac <- AICcmodavg::AICc(x)
    pars <- x$m$getPars()
 # }
  return(c(ac, pars))
}

allSum <- vector("list", length = length(allList))


for (i in 1:length(allList)){
  
  ll <- allList[[i]]
  gg <- ll[[1]]
  gg2 <- ll[[2]]
  gg3 <- ll[[3]]
  gg4 <- ll[[4]]
  gg5 <- ll[[5]]
  gg6 <- ll[[6]]
  
  nam <- paste(letters[i], "_breakpoint.jpeg", sep="")
  
  dat <- l2[[i]]
  if (any(dat$S == 0)){
    wh <- which(dat$S == 0)
    dat$S[wh] <- 0.25
    message("\n", "Zeros found in richness values: these have been converted to 0.25", "\n")
  }
  
  #dat <- arrange(dat, A)
  
  #get AICc and parameter values for each model object (excluding the one with no fit for i == 4[[5]])
  if (!i == 4){
    allSum[[i]] <- lapply(allList[[i]], mod_sum)
  } else {
    al2 <- allList[[i]]
    al2[[5]] <- NA
    allSum[[i]] <- lapply(allList[[i]], function(x){
      if (length(x) > 1){
        mod_sum(x)
      } else {
        rep(NA, 6)
      }})
  }

  allSum[[i]][[7]] <- pow(dat)
  #########################################
  jpeg(paste(nam), width = 25, height = 49, res = 200, units = "cm")
  
  par(mfrow = c(4, 2))
  
  #power model results
  
  par("mfg")
  dat2 <- sar_power(dat)
  plot(dat2, pcol = "black", lcol = "red")
  
  par("mfg")
  log.data =  log.data = data.frame(A = log10(dat$A), S = log10(dat$S))
  lp <- lm(S ~ A, data = log.data)
  dat2 <- cbind(dat, "Fitted" = lp$fitted.values)
  dat2 <- arrange(dat2, A)
  plot(log10(dat2$A), log10(dat2$S), col = "black", pch = 16)
  lines(log10(dat2$A), dat2$Fitted, col = "red")
  title("Power")
  
  #Dengler model results
  par("mfg")
  dat2 <- cbind(dat, "Fitted" = gg$m$fitted())
  dat2 <- arrange(dat2, A)
  plot(dat2$A, dat2$S, col = "black", pch = 16)
  lines(dat2$A, dat2$Fitted, col = "red")
  title("Breakpoint")
  
  par("mfg")
  dat2 <- cbind(dat, "Fitted" = gg2$m$fitted())
  dat2 <- arrange(dat2, A)
  plot(log10(dat2$A), log10(dat2$S), col = "black", pch = 16)
  lines(log10(dat2$A), dat2$Fitted, col = "red")
  title("Breakpoint")
  
  par("mfg")
  dat2 <- cbind(dat, "Fitted" = gg3$m$fitted())
  dat2 <- arrange(dat2, A)
  plot(dat2$A, dat2$S, col = "black", pch = 16)
  lines(dat2$A, dat2$Fitted, col = "red")
  title("Smooth4")
  
  par("mfg")
  dat2 <- cbind(dat, "Fitted" = gg4$m$fitted())
  dat2 <- arrange(dat2, A)
  plot(log10(dat2$A), log10(dat2$S), col = "black", pch = 16)
  lines(log10(dat2$A), dat2$Fitted, col = "red")
  title("Smooth4")
  
  par("mfg")
  if (!i == 4){ 
  dat2 <- cbind(dat, "Fitted" = gg5$m$fitted())
  dat2 <- arrange(dat2, A)
  plot(dat2$A, dat2$S, col = "black", pch = 16)
  lines(dat2$A, dat2$Fitted, col = "red")
  title("Smooth5")
  } else{
    plot(dat2$A, dat2$S, col = "black", pch = 16)
    title("Smooth5")
  }

  par("mfg")
  dat2 <- cbind(dat, "Fitted" = gg6$m$fitted())
  dat2 <- arrange(dat2, A)
  plot(log10(dat2$A), log10(dat2$S), col = "black", pch = 16)
  lines(log10(dat2$A), dat2$Fitted, col = "red")
  title("Smooth5")
  
  dev.off()

}


##format AICc and parameter estimates from each model
matz <- matrix(nrow = 0, ncol = 6)
colnames(matz) <- c("AICc", "c", "z1", "z2", "B", "k")

for (i in 1:length(allSum)){
  
  dum <- allSum[[i]][1:6]
  pres <- allSum[[i]][[7]]#power model results
  
  #fill in k for models with out it
  for (j in 1:4){
    dum[[j]] <- c(dum[[j]], NA)
  }
  m2 <- matrix(unlist(dum), byrow = TRUE, ncol = 6)
  
  m1A <- c(pres[1:3], rep(NA, 3))
  m1B <- c(pres[4:6], rep(NA, 3))
  
  m3 <- rbind(m1A, m1B, m2)
  matz <- rbind(matz, m3)
}

matz <- as.data.frame(matz)
matz <- apply(matz, 2, round, digits = 2)
matz <- cbind(matz, as.vector(unlist(vapply(1:9, function(x) rep(x, 8), FUN.VALUE = numeric(8)))))

matz <- cbind(matz, rep(c(1,2), (nrow(matz) / 2)))
matz <- as.data.frame(matz)
colnames(matz)[7:8] <- c("Sample", "Type")

matzNL <- filter(matz, Type == 1)
matzNL$Type <- rep(c("Power", "Normal_breakpoint", "Smooth_4", "smooth_5"), 9)
matzL <- filter(matz, Type == 2)
matzL$Type <- rep(c("Power", "Normal_breakpoint", "Smooth_4", "smooth_5"), 9)

rbind(matzNL, matzL) %>% wcs()




#############################################################
###RAW FUNCTIONS###########################
################################################
#normal breakpoint

S ~ 10 ^ (log10(c) +(log10(A) < log10(B)) * (z1 * log10(A)) + 
            (log10(A) >= log10(B)) * ((z1 * (log10(B) + z2)) * (log10(A) - log10(B))))

#log version
log10(S) ~ log10(c) + (log10(A) < log10(B)) * (z1 * log10(A)) + 
  (log10(A) >= log10(B)) * (z1 * log10(B) + z2 * (log10(A) - log10(B)))


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
















###################################################################

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






