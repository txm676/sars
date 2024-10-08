

countryside_optim <- function(dat, mod_nam = NULL, 
                              grid_start = "partial",
                              startPar){
  

  #to be generic it needs to build based on number of habitats
  #provided by the user
  CNs <- colnames(dat)
  Nhab <- length(which(grepl("Area", CNs)))
  
  if (is.null(startPar)){
  
  if (grid_start == "partial"){
    start.vec <- c(0.005,0.05, 1,5,100)
  } else if (grid_start == "exhaustive"){
    start.vec <- c(0.00005,0.005,0.05,0.5,1,5,100, 
                   1000, 10000, 100000)
  }
  start.list <- rep(list(start.vec), Nhab)
  start.list$z <- c(0.001, 0.01, 0.1, 0.5, 1)
  
  LNs <- sapply(1:Nhab, function(x) paste0("h", x))
  names(start.list) <- c(LNs, "z")

  grid.start <- expand.grid(start.list)
  
  #random for testing
  grid.start <-  grid.start[sample(1:nrow(grid.start), 150),]
  
  
  } else { #user provided start values
    
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
 
 #  
 # mod_nam2 <- formula("S ~ c*(h1 * Area1 + h2 * Area2 + 
 #                     h3 * Area3)^z")  
 
 
#  mod_nam2 <- switch(mod_nam,
         #            "Country_power" = formula(y),
            #         "jigsaw" = formula(S ~ (c1 * H^d) * ((A / H)^z)))
  
 #lower bounds (0 for hab variables; -Inf[default] for z;
 #but checking whether z should be zero also)
 x <- grid.start[1,]
 xl <- rep(0, (length(x) - 1))
 names(xl) <- names(x)[1:length(x)-1]
 xl["z"] <- -Inf
 
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

countryside_affinity <- function(mods, habNam){
  
  MN <- names(mods)
  r3 <- lapply(mods, function(x){
  modP <- x$m$getPars()
  wZ <- which(names(modP) == "z")
  modP2 <- modP[-wZ]
  scP <- modP2 / max(modP2)
  names(scP) <- habNam
  scP <- c(scP, (max(modP2)^modP[wZ]))
  names(scP)[length(scP)] <- "c"
  scP
  })
  return(r3)
}

#grid_start = if startPar not NULL, this is ignored.
#Warning that exhaustive can take a while.

#spNam = optional vector of species-group names (matching 
#the column order, otherwise takes names of sp columns)

#habNam = optional vector of habitat names (matching 
#the column order, otherwise just uses Habitat1 etc)

#startPar = if not null, needs to be a numeric matrix, where no
#of rows = number of habitats, and no of cols = no of species
#groups (including ubiquitous sp, if provided). Row and column
#order needs to match the column order of data (i.e., )

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
# SP <- matrix(rep(c(3.061161e+08, 2.104536e-01,
#                    1.074748e+00, 1.223700e-01),4),
#              ncol = 4) %>% t()
# f <- sar_countryside(data, ubiSp = TRUE,
#                      habNam = c("AG", "SH", "F"),
#                      startPar  = SP)


#ubiSp = T

sar_countryside <- function(data, modType = NULL,
                            grid_start = "partial",
                            startPar = NULL, 
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
  
  ##Rename columns
  cnD <- colnames(data)
  CN <- floor((ncol(data)) / 2)#if ubiSp, it will be 0.5 over (hence floor) 
  colnames(data)[1:3] <- sapply(1:CN, function(x) paste0("Area", x))
  colnames(data)[(CN + 1):(CN + CN)] <- sapply(1:CN, 
                                       function(x) paste0("SR", x))
  if (ubiSp) colnames(data)[ncol(data)] <- "SR_UB"

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
    if (!any(c("partial", "exhaustive") %in% grid_start)){
      stop("grid_start should be either 'partial' or 'exhaustive")
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
    CO <- countryside_optim(dat = dum,
                      grid_start = grid_start,
                      startPar = startPar2)
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
                                startPar = startPar2)
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
    aff <- countryside_affinity(res2, habNam = habNam)
    aff_H <- lapply(aff, function(x) x[1:(length(x) - 1)])
    aff_C <- sapply(aff, function(x) x[length(x)])
    
    res <- list(res, aff_H, aff_C)
  } else {
    res <- list(res, "All models NA - no affinity values")
  }
  
  class(res) <- c("habitat", "sars", "list")
  attr(res, "type") <- "countryside"
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

countryside_extrap <- function(fits, area, 
                               modType = NULL){
  
  #order of area values needs to match the order of 
  #the model fits in 'fits'
  
  if (!inherits(fits, "habitat")){
    stop("fits should be an object generated by sar_countryside()")
  }
  
  if (attributes(fits)$type != "countryside"){
    stop("fits should be an object generated by sar_countryside()")
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
  Pred <- vapply(fits, function(x){
    predict(x, area)
  }, FUN.VALUE = numeric(1))
  
  PredTot <- sum(Pred)
  
  resP <- list("Indiv_mods" = Pred, 
               "Total" = PredTot,
               "Failed_mods" = mes)
  return(resP)
}
