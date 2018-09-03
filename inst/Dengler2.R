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