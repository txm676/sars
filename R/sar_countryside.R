

countryside_startPars <- function(dat, sp_grp,
                                  grid_start, Nhab){
  
  A2 <- rowSums(dat[,1:(ncol(dat) - 1)])
  d2 <- data.frame("A" = A2, "S" = dat[,ncol(dat)])
  
  sp2 <- tryCatch(sar_power(d2),
                  error = function(e) NA)
  
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
  
  hmax <- c1 ^ (1/z1)
  
  if (grid_start == "partial"){
    start.vec <- c(0.0000000001,
                   0.000001,0.1,
                   5000,100000,
                   10000000, 100000000, 999)
  } else if (grid_start == "exhaustive"){
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
  if (grid_start == "partial"){
    start.list$z <- c(0.01, 0.1, 0.7, z1)
  } else if (grid_start == "exhaustive"){
    start.list$z <- c(0.001, 0.01, 0.1, 0.25,
                      0.5, 1, z1)  
  }
  
  LNs <- sapply(1:Nhab, function(x) paste0("h", x))
  names(start.list) <- c(LNs, "z")
  
  grid.start <- expand.grid(start.list)
  
  return(grid.start)
}

countryside_optim <- function(dat, mod_nam = NULL, 
                              grid_start = "partial",
                              startPar, z_lower = 0,
                              sp_grp){
  
  #to be generic it needs to build based on number of habitats
  #provided by the user
  CNs <- colnames(dat)
  Nhab <- length(which(grepl("Area", CNs)))
  
  if (is.null(startPar)){
  
    grid.start <- countryside_startPars(dat, sp_grp,
                                        grid_start, Nhab)
  
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
 y <- paste0("S ~ (", y, ")", "^z")
 mod_nam2 <- formula(y)
 
 #  
 # mod_nam2 <- formula("S ~ c*(h1 * Area1 + h2 * Area2 + 
 #                     h3 * Area3)^z")  
 
 
#  mod_nam2 <- switch(mod_nam,
         #            "Country_power" = formula(y),
            #         "jigsaw" = formula(S ~ (c1 * H^d) * ((A / H)^z)))
  
 #lower bounds (0 for hab variables and z
 x <- grid.start[1,]
 xl <- rep(0, length(x))
 names(xl) <- names(x)
 #can set to -Inf for full search of par space
 if (z_lower != 0) xl["z"] <- z_lower
 
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

#z_lower = the lower bound to be used for the z-parameter in the
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
# data(countryside)
# f <- sar_countryside(countryside, ubiSp = TRUE,
#                      habNam = c("AG", "SH", "F"))
# 
# dd <- f[[5]]
# 
# dd_Area <- length(which(grepl("Area", colnames(dd))))
# dd_SR <- ncol(dd) - dd_Area
# 
# dd_Ran <- range(dd[,1:dd_Area])
# 
# dd2_Area <- dd[,1:dd_Area]
# dd2_SR <- dd[,(dd_Area + 1):ncol(dd)]
# 
# #if predicted richness values returned
# if (length(f[[4]]) > 1){
# 
# ##total predicted richness for the actual
# #add total area and totR columns in
# dd3_Area <- as.data.frame(dd2_Area)
# dd3_Area$totA <- rowSums(dd3_Area)
# dd3_Area$totR <-  f[[4]]
# #same for main dataframe
# ddTot <- dd
# ddTot$totA <- rowSums(dd2_Area)
# ddTot$totR <- rowSums(dd2_SR)
# 
# if(!identical(ddTot$totA, dd3_Area$totA)){
#   stop("rownames mismatch in plot.countryside,",
#        " contact the package author")
# }
# 
# plot(ddTot$totR, dd3_Area$totR,
#      xlab = "Observed total richness",
#      ylab = "Predicted total richness",
#      ...)
# abline(0,1)
# 
# #extract power model values
# if (length(f[[7]]) > 1){
# points(f[[7]]$data$S, f[[7]]$calculated,
#        col = "red")
# } else {
#   cat("\n\nPower model could not be fitted\n\n")
#   }#eo if f7
# } else {
#   cat("\nNo predicted total richness values as some models could not be fitted\n\n")
# }
# 
# ##predicted curves for each land-use
# Ar_seq <- seq(dd_Ran[1], dd_Ran[2], 
#               length.out = 1000)
# #convert in N tables, where in each you can N columns,
# #where N = number of land-use types. In each all columns,
# #except the focal habitat are zeros
# ar_ls <- vector("list", length = dd_Area)
# totR_i <- matrix(ncol = dd_Area, nrow = length(Ar_sq))
# for (i in 1:dd_Area){
#   m_ls <- matrix(0, ncol = dd_Area,
#                        nrow = length(Ar_seq))
#   colnames(m_ls) <- names(f[[2]][[1]])
#   m_ls[,i] <- Ar_seq 
#   totR_i <- apply(m_ls[,1:dd_Area],1,function(x){
#     v <- as.vector(x)
#     vc <- countryside_extrap(f, area = v)
#     vc$Total
#   })
#   
# }
# 
# #cant plot empirical totals vs model preds,
# #as different landscapes can have same total
# #but different proportions of habitats, and thus
# #predicted richness differs
# 
# #check inside main function in the affininty bit,
# #as if some models not fitted does habNam then work?
# 
# for (i in 1:dd_Area){
#   d3 <- dd[,i]
#   u3 <- unique(d3)
#   m3 <- matrix(nrow = length(u3), ncol = dd_SR)
#   rownames(m3) <- u3
#   colnames(m3) <- colnames(dd)[(dd_Area + 1):ncol(dd)]
#   for (j in 1:length(u3)){
#     d4 <- subset(dd, dd[,i] == u3[j])
#     m3[j,] <- as.vector(colMeans(d4)[(dd_Area + 1):ncol(dd)])
#   }
# }
# 
# if(!all.equal(as.numeric(rownames(m3)), u3)){
#   stop("rownames mismatch in plot.countryside,",
#        " contact the package author")
# }
# 
# m4 <- m3[order(u3),]
# 
# plot(sort(u3), m4[,1], type = "l")

sar_countryside <- function(data, modType = NULL,
                            grid_start = "partial",
                            startPar = NULL,
                            z_lower = 0,
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
  if (length(z_lower) != 1 | !is.numeric(z_lower)){
    stop("z_lower should be a numeric vector of length 1")
  }
  
  ##Rename columns
  cnD <- colnames(data)
  CN <- floor((ncol(data)) / 2)#if ubiSp, it will be 0.5 over (hence floor) 
  colnames(data)[1:CN] <- sapply(1:CN, function(x) paste0("Area", x))
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
    #Sp.group number
    sgn <- x - CN
    CO <- countryside_optim(dat = dum,
                      grid_start = grid_start,
                      startPar = startPar2,
                      z_lower = z_lower,
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
                                startPar = startPar2,
                                z_lower = z_lower,
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
    aff <- countryside_affinity(res2, habNam = habNam)
    aff_H <- lapply(aff, function(x) x[1:(length(x) - 1)])
    aff_C <- sapply(aff, function(x) x[length(x)])
    
    res <- list(res, aff_H, aff_C)
  } else {
    res <- list(res, "All models NA - no affinity values", NA)
  }
  
  #Calculate total richness for each site: only do
  #if all models have fit
  if (("None" %in% FM)){
    TR <- apply(data[,1:CN],1, function(x){
    v <- as.vector(x)
    vc <- countryside_extrap(f, area = v)
    vc$Total
    })
    totA1 <- rowSums(data[,1:CN])
    totR1 <- rowSums(data[(CN + 1): (ncol(data))])
    dd_pow1 <- tryCatch(sar_power(data.frame("A" = totA1, 
                                   "R" = totR1)),
                        error = function(e) NA)
    if (length(dd_pow1) == 1){
      rss <- NA
    } else {
    ##Calculate RSS
      cs_rss <- sum((TR - totR1)^2)
      pow_rss <- dd_pow1$value
      rss <- c("Countryside_RSS" = cs_rss, 
               "Power_RSS" = pow_rss)
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

