

##Questions:

#In Pereira's script, their nls function calls have lower bounds,
#is this something that is essential?

#does it work with non-nested landscapes / islands?

countryside_optim <- function(mod_nam, startPar = NULL, data){
  

  #to be generic it needs to build based on number of habitats
  #provided by the user
  CNs <- colnames(data)
  Nhab <- length(which(grepl("Area", CNs)))
  
  if (is.null(startPar)){
  
  start.vec <- c(0.005,0.05,0.5,1,5)
  start.list <- rep(list(start.vec), Nhab)
  start.list$z <- c(0.001, 0.01, 0.1, 0.5, 1)
  
  LNs <- sapply(1:Nhab, function(x) paste0("h", x))
  names(start.list) <- c(LNs, "z")
  
  grid.start <- expand.grid(start.list)
  
  } else { #user provided start values
    
    if (length(startPar) != Nhab + 1){
      stop("Length of startPar must equal the number of",
           " habitats plus one")
    }
    
    grid.start <- as.data.frame(matrix(startPar, 
                                      nrow = 1))
    LNs <- sapply(1:Nhab, function(x) paste0("h", x))
    colnames(grid.start) <- c(LNs, "z")
  }
  
 y <- sapply(1:Nhab, function(x) paste0("h", x, "*", "Area", x))
 y <- toString(y)
 y <- gsub(", ", " + ", y)
 y <- paste0("S ~ (", y, ")", "^z")
 mod_nam2 <- formula(y)
  
  
#  mod_nam2 <- switch(mod_nam,
         #            "Country_power" = formula(y),
            #         "jigsaw" = formula(S ~ (c1 * H^d) * ((A / H)^z)))
  
  fit.list <- suppressWarnings(apply(grid.start, 1, function(x){
    tryCatch(minpack.lm::nlsLM(mod_nam2,
                               start = x,
                               control = minpack.lm::nls.lm.control(maxiter = 10000, 
                                                                    maxfev = 100000),
                               data = data),
             error = function(e) NA)
  }))
  
  len.fit.list <- sapply(fit.list, length)
  if (any(len.fit.list > 1)){ #otherwise all NA
    good.fit.list <- which(len.fit.list > 1)
    new.fit.list <- fit.list[good.fit.list]
    AIC.fit.list <- vapply(new.fit.list, AIC, 
                           FUN.VALUE = numeric(1))
    #if multiple min, it just picks the first
    best.fit <- new.fit.list[[which.min(AIC.fit.list)]]
  } else {
    best.fit <- NA 
  }
  return(best.fit)
}


sar_countryside <- function(data, modType = "power_log", 
                        con = NULL, logT = log){
  
  
  
}
