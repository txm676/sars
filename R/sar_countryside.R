

##Internal function to generate starting parameter estimates
##for countryside models
#' @noRd
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


##Internal function to fit countryside models using all
#starting parameter estimates and return the best fit
#' @importFrom stats AIC
#' @importFrom minpack.lm nlsLM nls.lm.control 
#' @noRd
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


##Internal function to return affinity and c values
#' @noRd
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


###############################################################
########Main Function to Fit Countryside Models#############
############################################################

#' Fit the countryside SAR model
#'
#' @description Fit the countryside biogeography SAR model in
#'   either power or logarithmic form, including an optional
#'   component model for ubiquitous species.
#' @usage sar_countryside(data, modType = "power",
#' gridStart = "partial", startPar = NULL, zLower = 0, 
#' ubiSp = FALSE, spNam = NULL, habNam = NULL)
#' @param data A dataset in the form of a dataframe – requires 
#' a specific column order (see 'Details' below).
#' @param modType Fit the power (\code{"power"}) or logarithmic
#'   (\code{"logarithmic"}) form of the countryside model.
#' @param gridStart The type of grid search procedure to be
#'   implemented to test multiple starting parameter values: can
#'   be one of \code{partial} (default) or \code{exhaustive}. If
#'   \code{startPar} is provided, this argument is ignored. Note 
#'   that \code{exhaustive} can be quite time consuming to run.
#' @param startPar Optional (default = NULL) starting parameter
#'   estimates for the constituent models. Must be a numeric
#'   matrix (see 'Details' below).
#' @param zLower The lower bound to be used for the z-parameter
#'   in the \link[minpack.lm]{nlsLM} function. Default is set to
#'   zero, but can be changed to any numeric value (e.g., -Inf to
#'   allow for a full search of parameter space).
#' @param ubiSp A logical argument specifying whether a component
#'   model should be fitted for ubiquitous species. If set to
#'   TRUE, a column of ubiquitous species richness must be
#'   included in \code{data}.
#' @param spNam Optional vector of species-group names (matching
#'   the column order in \code{data}, otherwise takes names of
#'   species columns in \code{data}).
#' @param habNam Optional vector of habitat names (matching the
#'   column order in \code{data}, otherwise just uses Habitat1
#'   etc).
#' @details To work, the countryside SAR model requires that all
#'   species in the study system have been classified based on
#'   the habitats present. For example, in a study system with
#'   two habitats (forest and grassland), all species must a
#'   priori have been classified as either forest species or
#'   grassland species; optionally, species can also be
#'   classified as ubiquitous species (i.e., species that do not
#'   have strong affinity for a specific habitat – true habitat
#'   generalists; controlled using the \code{ubiSp} argument).
#'   The provided input dataset (\code{data}) will typically
#'   relate to a series of landscapes with differing areas of the
#'   N habitats (e.g., forest and grassland), and for each
#'   landscape the number of (for example) forest species and
#'   grassland species present will be provided as well as
#'   (optionally) the number of ubiquitous species.
#' 
#' It is important that the column orders in \code{data} are
#' correct. The first set of columns should be the habitat area
#' columns, followed by the habitat species richness columns
#' (these should be in the same order as the area columns). An
#' optional final column of ubiquitous species richness can also
#' be included. For example, in a
#' dataset with two habitats (forest and grassland) and setting \code{ubiSp =
#' TRUE}, the column order in \code{data} could be:
#' forest-a, grassland-a, forest-s, grassland-s, and Ubi-s, where
#' a = area, s = species richness, and Ubi-s = the number of
#' ubiquitous species.
#' 
#' The countryside SAR model works by fitting individual
#' component models of a particular form (e.g., power), one for
#' each of the habitat types (e.g., one model for forest species,
#' one for grassland species, and so on). The predictions from
#' these component models are then combined to generate a total
#' predicted richness for each site / landscape. The
#' output of the model fitting includes the individual component
#' model fits, the total predicted (fitted) richness values for
#' each site, and the habitat affinity values for each species
#' group. The latter vary from zero to one, and equal 1 for a
#' given species group's affinity to its preferred habitat (e.g.,
#' forest species for forest).
#'   
#' For \code{startPar}, if not NULL, it needs to be a numeric
#' matrix, where number of rows = number of species groups
#' (including ubiquitous sp., if provided), and number of columns
#' equals number of habitats + 1. Matrix row order matches the
#' order of species group columns in \code{data}, and matrix
#' column order matches the order of habitat columns in
#' \code{data} + 1 extra final column for the z-parameter
#' estimate.
#' 
#' Two different types of plot can be generated with the output,
#' using \code{\link{plot.habitat}}. The
#' \code{\link{countryside_extrap}} function can be used with the
#' output of \code{sar_countryside} to predict the species
#' richness of landscapes with varying areas of the analysed
#' habitats.
#' 
#' See Matthews et al. (2025) for further details.
#' 
#' @return A list (of class ‘habitat’ and ‘sars’; and with a
#'   ‘type’ attribute of ‘countryside’) with eight elements: 
#'   \itemize{ 
#'    \item \strong{i.}  A list of the non-linear regression model
#'   fits for each of the species groups.
#'   \item \strong{ii.}  The habitat affinity values for each of
#'   the models in (i).
#'   \item \strong{iii.}   The c-parameter values for each of the
#'   models in (i).
#'   \item \strong{iv.}   The predicted total richness values
#'   (calculated by summing the predictions for each constituent
#'   countryside model) for each site in the dataset.
#'   \item \strong{v.}   The residual sum of squares – calculated
#'   using the predicted and observed total richness values – for
#'   both the countryside model and the Arrhenius power SAR model
#'   (or logarithmic model) to enable model comparison.
#'   \item \strong{vi.}  The dataset used for fitting (i.e., \code{data}).
#'   \item \strong{vii.}   The power (or logarithmic) model fit object.
#'   \item \strong{viii.}   The user-provided \code{ubiSp} argument.}
#'   
#' @note The model fits in (i) are objects of class ‘nls’,
#'   meaning that all the basic non-linear regression R methods
#'   can be applied (e.g., generating model summary tables or
#'   plotting the model residuals).
#' @references Matthews et al. (2025) In prep.
#' 
#' Pereira, H.M. & Daily, G.C. (2006) Modelling biodiversity
#' dynamics in countryside landscapes. Ecology, 87, 1877–1885.
#'   
#'   Proença, V. & Pereira, H.M. (2013) Species–area models to
#'   assess biodiversity change in multi-habitat landscapes: the
#'   importance of species habitat affinity. Basic and Applied
#'   Ecology, 14, 102–114.
#' @author Thomas J. Matthews, Inês Santos Martins, Vânia Proença 
#' and Henrique Pereira
#' @examples
#' data(countryside)
#' \dontrun{
#' #Fit the countryside SAR model (power form) to the data.
#' #Include a component model of ubiquitous species, and use the
#' #function’s starting parameter value selection procedure.
#' #Abbreviations: AG = agricultural land, SH = shrubland, F =
#' #oak forest, UB = ubiquitous species.
#' s3 <- sar_countryside(data = countryside, modType = "power",
#' gridStart = "partial", ubiSp = TRUE, habNam = c("AG", "SH",
#' "F"), spNam = c("AG_Sp", "SH_Sp", "F_Sp", "UB_Sp"))
#' 
#' #Predict the richness of a site which comprises 1000 area units
#' #of agricultural land, 1000 of shrubland and 1000 of forest.
#' countryside_extrap(s3, area = c(1000, 1000, 1000))
#' 
#' #Generate a plot of the countryside model’s predicted total
#' #richness vs. the observed total richness, and include the
#' #predictions of the Arrhenius power model
#' 
#' plot(s3, type = 1, powFit = TRUE)
#' 
#' #Plot the fitted individual SAR curves for each species group,
#' #providing set line colours, including a legend and
#' #positioning it outside the main plotting window, and modifying
#' #other aspects of the plot using the standard base R plotting
#' #commands
#' par(mar=c(5.1, 4.1, 4.1, 7.5), xpd=TRUE)
#' plot(s3, type = 2, lcol = c("black", "aquamarine4",
#' "#CC661AB3" , "darkblue"), pLeg = TRUE,
#' legPos ="topright", legInset = c(-0.27,0.3), lwd = 1.5)
#' 
#' #Provide starting parameter estimates for the component models
#' #instead of using gridStart
#' M2 <- matrix(c(3.061e+08, 2.105e-01, 1.075e+00, 1.224e-01,
#' 3.354e-08, 5.770e+05, 1.225e+01, 1.090e-01,
#' 6.848e-01, 1.054e-01, 4.628e+05, 1.378e-01,
#' 0.20747, 0.05259, 0.49393, 0.18725), nrow = 4,
#' byrow = TRUE)
#'
#' s4 <- sar_countryside(data = countryside,
#'                     modType = "power",
#'                    startPar = M2, ubiSp = TRUE)
#' }
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
                                gridStart = gridStart,
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



#' Use a sar_countryside() model object to predict richness
#'
#' @description Use a fitted model object from sar_countryside() 
#' to predict richness, given a set of habitat area values.
#' @usage countryside_extrap(fits, area)
#' @param fits A fitted model object from \code{\link{sar_countryside}}.
#' @param area A vector of area values - the number (and order)
#'   of area values (i.e., the length of the vector) must match
#'   the number (and order) of habitats in the dataset used in
#'   the \code{\link{sar_countryside}} fit.
#' @details Takes a model fit generated using
#'   \code{\link{sar_countryside}} and uses it to predict
#'   richness values for a set of user-provided habitat
#'   \code{area} values. Note this can either be interpolated or
#'   extrapolated predictions, depending on the range of area
#'   values used in the original model fits. A ubiquitous model
#'   prediction is included if a ubiquitous component model is
#'   included in \code{fits}.
#'
#'   The habitat area values provided through \code{area} need to
#'   be in the same order as the habitat columns in the original
#'   dataset used in \code{\link{sar_countryside}}.
#'   
#'   The function does work with failed component model fits (any model fit that
#'   is NA is removed along with the corresponding area values provided by the
#'   user), as long as at least one component model was successfully fitted.
#'   However, arguably it does not make sense to predict richness values unless
#'   all component models were successfully fitted.
#'
#' @return A list with three elements. The first contains the
#'   predicted richness values from the individual component
#'   models. The second contains the predicted total richness of
#'   the site (i.e., the summed component model predictions), and
#'   the third is a logical value highlighting whether there were
#'   any failed models in \code{fits}, i.e., component models
#'   that could not be fitted in \code{\link{sar_countryside}}.
#' @author Thomas J. Matthews
#' @examples
#' \dontrun{
#' data(countryside)
#' #Fit the sar_countryside model (power version)
#' s3 <- sar_countryside(data = countryside, modType = "power",
#' gridStart = "partial", ubiSp = TRUE, habNam = c("AG", "SH",
#' "F"), spNam = c("AG_Sp", "SH_Sp", "F_Sp", "UB_Sp"))
#' #Predict the richness of a site which comprises 1000 area units
#' #of agricultural land, 1000 of shrubland and 1000 of forest.
#' countryside_extrap(s3, area = c(1000, 1000, 1000))
#' }
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

