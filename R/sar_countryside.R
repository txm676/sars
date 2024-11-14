

countryside_startPars <- function(dat, modType,
                                  sp_grp,
                                  gridStart, Nhab){
  
  A2 <- rowSums(dat[,1:(ncol(dat) - 1)])
  d2 <- data.frame("A" = A2, "S" = dat[,ncol(dat)])
  
  if (modType == "logarithmic"){
    sp2 <- tryCatch(sar_loga(d2),
                  error = function(e) NA)
  } else {
    sp2 <- tryCatch(sar_power(d2),
                    error = function(e) NA)
  }
  
  if (length(sp2) == 1){
    c1 <- 5
    z1 <- 0.25
  } else {
    if (length(sp2$par) != 2){
      c1 <- 5
      z1 <- 0.25
    } else {
      c1 <- sp2$par[1]
      z1 <- sp2$par[2]
    }
  }
  
  if (modType == "logarithmic"){
    hmax <- exp(c1/z1)
  } else {
    hmax <- c1 ^ (1/z1)
  }
  
  if (gridStart == "partial"){
    start.vec <- c(0.0000000001,
                   0.000001,0.1,
                   5000,100000,
                   10000000, 100000000, 999)
  } else if (gridStart == "exhaustive"){
    start.vec <- c(0.0000000001,0.00000001,
                   0.000001,0.0001,0.01,
                   1,1000, 10000,100000,1000000,
                   10000000, 100000000, 
                   10000000000, 999)
  }
  start.list <- rep(list(start.vec), Nhab)
  #add in the calculated value
  if (!is.null(sp_grp)){
  start.list[[sp_grp]][length(start.vec)] <- hmax 
  } else {
    #for ubiquitous, we don't know which group is
    #hmax, so add it to each element.
    start.list <- lapply(start.list, function(x){
      x[length(start.vec)] <- hmax 
      x
    })
  }
  #include the calculated value
  if (gridStart == "partial"){
    start.list$z <- c(0.01, 0.1, 0.7, z1)
  } else if (gridStart == "exhaustive"){
    start.list$z <- c(0.001, 0.01, 0.1, 0.25,
                      0.5, 1, z1)  
  }
  
  LNs <- sapply(1:Nhab, function(x) paste0("h", x))
  names(start.list) <- c(LNs, "z")
  
  grid.start <- expand.grid(start.list)
  
  return(grid.start)
}

countryside_optim <- function(dat, modType, 
                              gridStart = "partial",
                              startPar, zLower = 0,
                              sp_grp){
  
  #to be generic it needs to build based on number of habitats
  #provided by the user
  CNs <- colnames(dat)
  Nhab <- length(which(grepl("Area", CNs)))
  
  if (is.null(startPar)){
  
    grid.start <- countryside_startPars(dat, modType = modType,
                                        sp_grp,
                                        gridStart, Nhab)
  
  #random for testing
 # grid.start <-  grid.start[sample(1:nrow(grid.start), 150),]
  
  } else { #user provided start values
    
    grid.start <- as.data.frame(matrix(startPar, 
                                      nrow = 1))
    LNs <- sapply(1:Nhab, function(x) paste0("h", x))
    colnames(grid.start) <- c(LNs, "z")
  }
  
 y <- sapply(1:Nhab, function(x) paste0("h", x, "*", "Area", x))
 y <- toString(y)
 y <- gsub(", ", " + ", y)
 y <- switch(modType,
          "logarithmic" = paste0("S ~ z * log(", y, ")"),
          "power" = paste0("S ~ (", y, ")", "^z"))
 mod_nam2 <- formula(y)
 
 #lower bounds (0 for hab variables and z
 x <- grid.start[1,]
 xl <- rep(0, length(x))
 names(xl) <- names(x)
 #can set to -Inf for full search of par space
 if (zLower != 0) xl["z"] <- zLower
 
  fit.list <- suppressWarnings(apply(grid.start, 1, function(x){
    tryCatch(minpack.lm::nlsLM(mod_nam2,
                               start = x,
                               lower = xl,
                               control = minpack.lm::nls.lm.control(maxiter = 1000, 
                                                                    maxfev = 100000),
                               data = dat),
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

countryside_affinity <- function(mods, modType,
                                 habNam){
  
  MN <- names(mods)
  r3 <- lapply(mods, function(x){
  modP <- x$m$getPars()
  wZ <- which(names(modP) == "z")
  modP2 <- modP[-wZ]
  scP <- modP2 / max(modP2)
  names(scP) <- habNam
  if (modType == "logarithmic"){
    scP <- c(scP, (log(max(modP2))*modP[wZ]))
  } else {
    scP <- c(scP, (max(modP2)^modP[wZ]))
  }
  names(scP)[length(scP)] <- "c"
  scP
  })
  return(r3)
}

#gridStart = if startPar not NULL, this is ignored.
#Warning that exhaustive can take a while.

#spNam = optional vector of species-group names (matching 
#the column order, otherwise takes names of sp columns)

#habNam = optional vector of habitat names (matching 
#the column order, otherwise just uses Habitat1 etc)

#startPar = if not null, needs to be a numeric matrix, where no
#of rows = number of habitats, and no of cols = no of species
#groups (including ubiquitous sp, if provided). Row and column
#order needs to match the column order of data (i.e., )

#zLower = the lower bound to be used for the z-parameter in the
#minpack.lm::nlsLM function. Default is set to zero, but can be 
#changed to any numeric value (e.g., -Inf to allow for a full
#search of parameter space)

#@return #in help file need to explain how to 
#extract raw nls fit objects:
#(i) a list of the non-linear regression model fits for each of
#the i species groups; (ii) the habitat affinity values (see Eq.
#*) for each of the models in (i); and (iii) the c-parameter
#values (see Eq.*) for each of the models in (i). The model fits
#in (i) are objects of class ‘nls’, meaning that all the basic
#non-linear regression R methods can be applied (e.g., generating
#model summary tables or plotting the model residuals).

#Example:
#user-provided starting pars: 1st row = Spcs_AG; and
#columns relate to h paramters for AG, SH, and QF, and
#then z.
# 
# 
#' @export
sar_countryside <- function(data,
                            modType = "power",
                            gridStart = "partial",
                            startPar = NULL,
                            zLower = 0,
                            ubiSp = FALSE,
                            spNam = NULL,
                            habNam = NULL){
  
  if (!(is.matrix(data) | is.data.frame(data)))
    stop('data must be a matrix or dataframe')
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop('NAs present in data')
  if (!is.logical(ubiSp) | anyNA(ubiSp)){
      stop("ubiSp should be a logical vector of length 1 (column no.)")
  } else if (isTRUE(ubiSp)){
    if (ncol(data) %% 2 == 0){
      stop("If ubiSp == TRUE, there should be an odd number of columns")
    } 
  }
  if (length(zLower) != 1 | !is.numeric(zLower)){
    stop("zLower should be a numeric vector of length 1")
  }
  if (!modType %in% c("power", "logarithmic")){
    stop("modType should be one of power or logarithmic")
  }
  
  ##Rename columns
  cnD <- colnames(data)
  CN <- floor((ncol(data)) / 2)#if ubiSp, it will be 0.5 over (hence floor) 
  colnames(data)[1:CN] <- sapply(1:CN, function(x) paste0("Area", x))
  colnames(data)[(CN + 1):(CN + CN)] <- sapply(1:CN, 
                                       function(x) paste0("SR", x))
  if (ubiSp) colnames(data)[ncol(data)] <- "SR_UB"

  if (any(rowSums(data[,1:CN]) == 0)){
    if (modType == "logarithmic"){
    stop("Some sites have total area equal to zero - this is ",
         "not possible with the logarithmic model")
    } else {
      warning("Some sites have total area equal to zero")
    }
  }
  
  #if habNam & spNam provided, check in correct format, otherwise take
  #area / species column names
  if (ubiSp){
    CN2 <- CN + 1
  } else{
    CN2 <- CN
  }
  
  if (!is.null(spNam)){
    if (!is.vector(spNam) |
        !is.character(spNam) | (length(spNam)!= CN2)){
      stop("spNam should be character vector of\n sp. group names",
           " of correct length.")
    }
  } else{
    spNam <- cnD[(CN + 1):length(cnD)]
  }
  if (!is.null(habNam)){
    if (!is.vector(habNam) |
        !is.character(habNam) | (length(habNam)!= CN)){
      stop("habNam should be character vector of\n habitat names",
           " of correct length.")
    }
  } else{
    habNam <- sapply(1:CN, function(x) paste0("Habitat", x))
  }
  
  if (!is.null(startPar)){
    ###needs to be a matrix, with each row corresponding to
    #habitat types in matching column order
    if (!is.matrix(startPar)){
      stop("startPar should be a matrix")
    } else { #+1 is for z
      if (!all(dim(startPar) == c(CN2, (CN + 1)))){
        stop("Dimensions of startPar are incorrect")
      }
      if (!is.numeric(startPar) | anyNA(startPar)){
        stop("startPar should contain only numbers and no NAs")
      }
    }
  } else {
    if (!any(c("partial", "exhaustive") %in% gridStart)){
      stop("gridStart should be either 'partial' or 'exhaustive")
    }
  }#eo is.null(startPar)
  
  ##Need to then fit the models for each SR, including for 
  #SR_UB if ubiSp.
  k <- 1
  res <- lapply(((CN + 1):(CN + CN)), function(x){
    dum <- data[,c(1:CN, x)]
    dum <- dum[order(dum[,ncol(dum)]),]
    colnames(dum)[ncol(dum)] <- "S"
    if (!is.null(startPar)){
      startPar2 <- startPar[k,]
    } else {
      startPar2 <- startPar
    }
    #Sp.group number
    sgn <- x - CN
    CO <- countryside_optim(dat = dum,
                      modType = modType,
                      gridStart = gridStart,
                      startPar = startPar2,
                      zLower = zLower,
                      sp_grp = sgn)
    k <<- k + 1
    CO
  })
  
  if (ubiSp){
    dum <- data[,c(1:CN, ncol(data))]
    dum <- dum[order(dum[,ncol(dum)]),]
    colnames(dum)[ncol(dum)] <- "S"
    if (!is.null(startPar)){
      startPar2 <- startPar[k,]
    } else {
      startPar2 <- startPar
    }
    res$UB <- countryside_optim(dat = dum, 
                                modType = modType,
                                startPar = startPar2,
                                zLower = zLower,
                                sp_grp = NULL)
  }
  
  names(res) <- spNam
  
  len <- sapply(res, length)
  if (all(len == 1)){
    FM <- "All"
  } else if (any(len == 1)){
    w1 <- which(len == 1)
    FM <- w1
  } else {
    FM <- "None"
  }

  ##Calculate affinity values
  #First remove NA models,
  if (!("All" %in% FM)){
    res2 <- res
    if (!("None" %in% FM)){
      res2[w1] <- NULL
    }
    aff <- countryside_affinity(res2, modType = modType,
                                habNam = habNam)
    aff_H <- lapply(aff, function(x) x[1:(length(x) - 1)])
    aff_C <- sapply(aff, function(x) x[length(x)])
    
    res <- list(res, aff_H, aff_C)
  } else {
    res <- list(res, "All models NA - no affinity values", NA)
  }
  
  #Calculate total richness for each site: only do
  #if all models have fit
  if (("None" %in% FM)){
    res2 <- res
    class(res2) <- c("habitat", "sars", "list")
    attr(res2, "type") <- "countryside"
    attr(res2, "modType") <- modType
    TR <- apply(data[,1:CN],1, function(x){
    v <- as.vector(x)
    vc <- countryside_extrap(res2, area = v)
    vc$Total
    })
    totA1 <- rowSums(data[,1:CN])
    totR1 <- rowSums(data[(CN + 1): (ncol(data))])
    modF <- ifelse(modType == "logarithmic", sar_loga,
                   sar_power)
    dd_pow1 <- tryCatch(modF(data.frame("A" = totA1, 
                                   "R" = totR1)),
                        error = function(e) NA)
    if (length(dd_pow1) == 1){
      rss <- NA
    } else {
    ##Calculate RSS
      cs_rss <- sum((TR - totR1)^2)
      pow_rss <- dd_pow1$value
      rss <- c("Countryside_RSS" = cs_rss, 
               pow_rss)
      names(rss)[2] <- ifelse(modType == "logarithmic",
                              "Logarithmic_RSS", "Power_RSS")
    }#eo length dd pow
  } else {
    TR <- "No predicted total richness values as some models could not be fitted"
    rss <- NA
    dd_pow1 <- NA
  }
  
  res[[4]] <- TR
  res[[5]] <- rss
  res[[6]] <- data
  res[[7]] <- dd_pow1
  res[[8]] <- ubiSp
  names(res) <- c("fits", "affinity", "c", "Pred.Tot.Rich",
                  "rss", "data", "pow.model", "ubiSp")
  class(res) <- c("habitat", "sars", "list")
  attr(res, "type") <- "countryside"
  attr(res, "modType") <- modType
  attr(res, "failedMods") <- FM
  return(res)
}


#If any model fits are NA, these are removed along
#with the corresponding area values provided by
#the user (arguably if this is the case, the extrapolations
#are not of use)

#the order of values in 'area' must match the order of
#habitat cols in the original data matrix provided to
#sar_countryside
#' @export
countryside_extrap <- function(fits, area){
  
  #order of area values needs to match the order of 
  #the model fits in 'fits'
  
  if (!inherits(fits, "habitat")){
    stop("fits should be an object generated by sar_countryside()")
  }
  
  if (attributes(fits)$type != "countryside"){
    stop("fits should be an object generated by sar_countryside()")
  }
  
  if (sum(area) == 0 & 
      attributes(fits)$modType == "logarithmic"){
      stop("Provided areas sum to zero - this is ",
           "not possible with the logarithmic model")
    } 
  
  fits <- fits[[1]]
  
  #If any fits are NA, we need to remove these and also
  #the corresponding user-provided area value(s)
  len <- sapply(fits, length)
  if (all(len == 1)) stop("All model fits are NA")
  if (any(len == 1)){
    warning("Some elements in 'fits' are NA; ",
           "\nthese have been removed prior to extrapolation")
    w1 <- which(len == 1)
    fits[w1] <- NULL
    mes <- TRUE
  } else {
    mes <- FALSE
  }
  
  #number of habitat is no. of parameters - 1 (as one is z)
  Nhab <- length(fits[[1]]$m$getPars()) - 1
  if (length(area) != Nhab) {
    stop("Length of 'area' does not equal no. of habitats in 'fits'")
  }
  names(area) <- sapply(1:Nhab, function(x) paste0("Area", x))
  area <- as.list(area)
  
  #run predict() for each model in fits
  #For each component model, this predicts the 
  #total number of species in that group in the landscape,
  #i.e., across all habitats in the landscape. In the plot
  #function we set the other habitats to zero, but here they
  #can be non-zero.
  Pred <- vapply(fits, function(x){
    predict(x, area)
  }, FUN.VALUE = numeric(1))
  
  PredTot <- sum(Pred)
  
  resP <- list("Indiv_mods" = Pred, 
               "Total" = PredTot,
               "Failed_mods" = mes)
  return(resP)
}

