

##Internal function to generate starting parameter estimates
##for countryside models
#' @noRd
countryside_startPars <- function(dat, modType,
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
  } else if (gridStart == "none"){
    start.vec <- c(0.000001, 1, hmax)
  }
  start.list <- rep(list(start.vec), Nhab)
  #add in the calculated value - we don't know which group is
    #hmax, so add it to each element.
    start.list <- lapply(start.list, function(x){
      x[length(start.vec)] <- hmax 
      x
    })

  #include the calculated value
  if (gridStart == "partial"){
    start.list$z <- c(0.01, 0.1, 0.7, z1)
  } else if (gridStart == "exhaustive"){
    start.list$z <- c(0.001, 0.01, 0.1, 0.25,
                      0.5, 1, z1)  
  } else if (gridStart == "none"){
    start.list$z <- z1
  }
  
  LNs <- sapply(1:Nhab, function(x) paste0("h", x))
  names(start.list) <- c(LNs, "z")
  
  grid.start <- expand.grid(start.list)
  
  #filter out rows where hmax not in or in more
  #than once
  if (gridStart == "none" & hmax > 5){
    GSN1 <- rowSums(grid.start)
    GSN2 <- which(GSN1 > hmax & GSN1 < (hmax*2))
    grid.start <- grid.start[GSN2,]
  }
  
  return(grid.start)
}


##Internal function to fit countryside models using all
#starting parameter estimates and return the best fit
#' @importFrom stats AIC
#' @importFrom minpack.lm nlsLM nls.lm.control 
#' @noRd
countryside_optim <- function(dat, modType, 
                              gridStart = "partial",
                              startPar, zLower = 0){
  
  #to be generic it needs to build based on number of habitats
  #provided by the user
  CNs <- colnames(dat)
  Nhab <- length(which(grepl("Area", CNs)))
  
  if (is.null(startPar)){
  
    grid.start <- countryside_startPars(dat, modType = modType,
                                      #  sp_grp,
                                        gridStart, Nhab)
  
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
#'   either power or logarithmic form.
#' @usage sar_countryside(data, modType = "power", gridStart =
#'   "partial", startPar = NULL, zLower = 0, habNam = NULL, spNam
#'   = NULL)
#' @param data A dataset in the form of a dataframe, with columns
#'   for habitat area values and species richness values –
#'   requires a specific column order (see 'Details' below).
#' @param modType Fit the power (\code{"power"}) or logarithmic
#'   (\code{"logarithmic"}) form of the countryside model.
#' @param gridStart The type of grid search procedure to be
#'   implemented to test multiple starting parameter values: can
#'   be one of \code{partial} (default), \code{exhaustive} or
#'   \code{none}. If \code{startPar} is provided, this argument
#'   is ignored. Note that \code{exhaustive} can be quite time
#'   consuming to run. In contrast, \code{none} is much quicker
#'   but only checks a very small number of starting paramter
#'   values (technically not "none").
#' @param startPar Optional (default = NULL) starting parameter
#'   estimates for the constituent models. Must be a numeric
#'   matrix (see 'Details' below).
#' @param zLower The lower bound to be used for the z-parameter
#'   in the \link[minpack.lm]{nlsLM} function. Default is set to
#'   zero, but can be changed to any numeric value (e.g., -Inf to
#'   allow for a full search of parameter space).
#' @param habNam Either a vector of habitat names (must be the
#'   same length as the number of habitat area columns in
#'   \code{data}, and in the same order as the area columns), or
#'   the habitat area column numbers in \code{data}.
#' @param spNam Either a vector of species group names (must be
#'   the same length as the number of species richness columns in
#'   \code{data}, and in the same order as the richness columns),
#'   or the species richness column numbers in \code{data}.
#' @details The provided input dataset (\code{data}) will
#'   typically relate to a series of landscapes (sites) with
#'   differing areas of N habitats (e.g., forest and grassland),
#'   and for each landscape the number of species in a priori
#'   defined groups.
#'
#'   To work, the countryside SAR model requires that all species
#'   in the study system have been classified into groups. This
#'   is typically done based on the habitats present in the study
#'   system. For example, in a study system with two habitats
#'   (forest and grassland), the species can be a priori
#'   classified as either forest species or grassland species.
#'   Optionally, species can also be classified as ubiquitous
#'   species (i.e., species that do not have strong affinity for
#'   a specific habitat – true habitat generalists). However, the
#'   model is flexible and species can technically be grouped
#'   into any groups. For example, in a study system with three
#'   habitats (forest, grassland, wetlands), species could be
#'   grouped into two groups: forest species and other species.
#'   Note that species must be classified prior to fitting the
#'   model, but the data can still be used to help guide these
#'   classifications.
#'
#'   It is important that the column orders in \code{data} are
#'   correct. The first set of columns should be all the habitat
#'   area columns, followed by all the group species richness
#'   columns. Within these two sets (i.e., area and richness
#'   columns), the order of columns is not important. The
#'   user must make clear which columns are the area columns and
#'   which the richness columns, using the \code{habNam} and
#'   \code{spNam} arguments. These can either provide habitat and
#'   species group names (e.g., \code{habNam = c("Forest",
#'   "Other")}) or the column numbers in \code{data} (e.g.,
#'   \code{spNam = 4:6}). If names are provided, note that these 
#'   names can be different to the column names in \code{data},
#'   but they need to be in the same order as their respective 
#'   columns in \code{data}. 
#'   
#'   No columns should be present in \code{data} before the area
#'   columns (i.e., the first column must be an area column) and
#'   all columns after the last species richness column are
#'   excluded by the function. And do not use the arguments to
#'   re-order columns as this will not be undertaken, e.g. use
#'   4:6 or c(4,5,6), and not c(4,6,5). If \code{habNam} and
#'   \code{spNam} are numeric (i.e., column numbers), the
#'   habitats and species groups are named Habitat1, Habitat2,
#'   and Sp_grp1, Sp_grp2, and so on, in the output.
#'
#'   The countryside SAR model works by fitting individual
#'   component models of a particular form (e.g., power), one for
#'   each of the species groups (e.g., one model for forest
#'   species, one for grassland species, and so on). The
#'   predictions from these component models are then combined to
#'   generate a total predicted richness for each site /
#'   landscape. The output of the model fitting includes the
#'   individual component model fits, the total predicted
#'   (fitted) richness values for each site, and the habitat
#'   affinity values for each species group. The latter vary from
#'   zero to one, and equal 1 for a given species group's
#'   affinity to its preferred habitat (e.g., forest species for
#'   forest).
#'   
#'   Note that the logarithmic model can generate negative fitted
#'   richness values for small areas in some cases.
#'   
#'   If you find some or all of your component models are not
#'   fitting / converging, you can try using \code{gridStart =
#'   "exhaustive"} to undertake a wider search of parameter space.
#'   If that still doesn't work you will need to provide a wide
#'   range of starting parameter values manually using the
#'   \code{startPar} argument. To speed up, you can try
#'   \code{gridStart = "none"}, which typically runs in seconds,
#'   but does not provide much of a search of starting parameter
#'   values.
#'
#'   For \code{startPar}, if not NULL, it needs to be a numeric
#'   matrix, where number of rows = number of species groups, and
#'   number of columns equals number of habitats + 1. Matrix row
#'   order matches the order of species group columns in
#'   \code{data}, and matrix column order matches the order of
#'   habitat columns in \code{data} + 1 extra final column for
#'   the z-parameter estimate.
#'   
#'   Three different types of plot can be generated with the
#'   output, using \code{\link{plot.habitat}}. The
#'   \code{\link{countryside_extrap}} function can be used with
#'   the output of \code{sar_countryside} to predict the species
#'   richness of landscapes with varying areas of the analysed
#'   habitats.
#'
#'   See Matthews et al. (2025) for further details.
#' 
#' @return A list (of class ‘habitat’ and ‘sars’; and with a
#'   ‘type’ attribute of ‘countryside’) with eight elements: 
#'   \itemize{ 
#'   \item \strong{i.}  A list of the non-linear regression model
#'   fits for each of the species groups. In the model output,
#'   the h coefficients follow the order of the habitat area
#'   columns in \code{data} (e.g., h1 = column 1). \item
#'   \strong{ii.}  The habitat affinity values for each of the
#'   models in (i).
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
#'   \item \strong{viii.}   The \code{habNam} and \code{spNam} vectors.}
#'   
#' @note The model fits in (i) are objects of class ‘nls’,
#'   meaning that all the basic non-linear regression R methods
#'   can be applied (e.g., generating model summary tables or
#'   plotting the model residuals). This also means that
#'   information criteria values can be returned for each
#'   component model, simply by using, for example,
#'   \code{\link[stats]{AIC}}. This can then be compared with
#'   equivalent values from, for example, the power model (see
#'   Examples, below). However, importantly note that while the
#'   values returned from \code{\link[stats]{AIC}} and
#'   \code{\link{sar_power}} are comparable, these values are not
#'   comparable with the AIC / AICc values presented in Proença &
#'   Pereira (2013) and related studies, due to the different
#'   information criteria equations used (although the delta
#'   values (calculated using a given equation) are comparable
#'   across equations). For more information, see the package
#'   vignette.
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
#' #Fit the countryside SAR model (power form) to the data, 
#' #which contrains 3 habitat types and 4 species groups.
#' #Use the function’s starting parameter value selection procedure.
#' #Abbreviations: AG = agricultural land, SH = shrubland, F =
#' #oak forest, UB = ubiquitous species.
#'  s3 <- sar_countryside(data = countryside, modType = "power",
#'  gridStart = "partial", habNam = c("AG", "SH",
#'  "F"), spNam = c("AG_Sp", "SH_Sp", "F_Sp", "UB_Sp"))
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
#' #Generate Type 2 & 3 plots providing set line colours, plot
#' #titles, and modifying other aspects of the plot using the
#' #standard #ase R plotting commands. See ?plot.habitat for more
#' #info
#' 
#'  plot(s3, type = 2, lcol = c("black", "aquamarine4",
#' "#CC661AB3" , "darkblue"), pLeg = TRUE, lwd = 1.5, 
#'  ModTitle = c("Agricultural land", "Shrubland", "Forest"))
#'  
#'  plot(s3, type = 3, lcol = c("black", "aquamarine4",
#' "#CC661AB3" , "darkblue"), pLeg = TRUE, lwd = 1.5, 
#'  ModTitle = c("Agricultural land", "Shrubland", "Forest"))
#'  
#' #Calculate AIC for a component model and compare with the 
#' #power model
#'  AIC(s3$fits$AG_Sp)
#'  SA <- rowSums(countryside[,1:3])#total site area
#'  SR <- countryside[,4] #agriculture column
#'  SP <- sar_power(data.frame(SA, SR))
#'  SP$AIC
#'
#' #Provide starting parameter estimates for the component models
#' #instead of using gridStart.
#' M2 <- matrix(c(3.061e+08, 2.105e-01, 1.075e+00, 1.224e-01,
#' 3.354e-08, 5.770e+05, 1.225e+01, 1.090e-01,
#' 6.848e-01, 1.054e-01, 4.628e+05, 1.378e-01,
#' 0.20747, 0.05259, 0.49393, 0.18725), nrow = 4,
#' byrow = TRUE)
#'
#' #Provide column numbers rather than names
#' s4 <- sar_countryside(data = countryside,
#'                     modType = "power",
#'                    startPar = M2,
#'                    habNam = 1:3, spNam = 4:7)
#'                    
#' #Speed up by trying gridStart = "none"
#'  s5 <- sar_countryside(data = countryside, modType = "power",
#'  gridStart = "none", habNam = c("AG", "SH",
#'  "F"), spNam = c("AG_Sp", "SH_Sp", "F_Sp", "UB_Sp"))
#' }
#' @export

sar_countryside <- function(data,
                            modType = "power",
                            gridStart = "partial",
                            startPar = NULL,
                            zLower = 0,
                            habNam = NULL,
                            spNam = NULL){
  
  if (!(is.matrix(data) | is.data.frame(data)))
    stop('data must be a matrix or dataframe')
  if (is.matrix(data)) data <- as.data.frame(data)
  if (anyNA(data)) stop('NAs present in data')
  if (length(zLower) != 1 | !is.numeric(zLower)){
    stop("zLower should be a numeric vector of length 1")
  }
  if (!modType %in% c("power", "logarithmic")){
    stop("modType should be one of power or logarithmic")
  }
  
  if (is.null(habNam) | is.null(spNam) | 
      length(c(habNam, spNam)) != ncol(data) | 
      !(is.numeric(habNam) | is.character(habNam)) |
      !(is.numeric(spNam) | is.character(spNam))){
    stop("habNam & spNam should be either character vectors\n",
         " of habitat / species group names, or numeric vectors\n",
         " of column numbers, of correct length")
  } 
  
  if (is.numeric(habNam)){
    habNam <- sapply(1:length(habNam), function(x) paste0("Habitat", x))
  }
  if (is.numeric(spNam)){
    if (min(spNam) <= length(habNam)) stop("spNam columns must come after area columns")
    spNam <- sapply(1:length(spNam), function(x) paste0("Sp_grp", x))
  }
  #remove any excess columns after richness cols
  data <- data[,1:(length(c(habNam, spNam)))]
  
  ##Rename columns
  cnD <- colnames(data)
  CN <- length(habNam)
  CN2 <- length(spNam)
  colnames(data)[1:CN] <- sapply(1:CN, function(x) paste0("Area", x))
  colnames(data)[(CN + 1):ncol(data)] <- sapply(1:CN2, 
                                       function(x) paste0("SR", x))
  if (any(rowSums(data[,1:CN]) == 0)){
    if (modType == "logarithmic"){
    stop("Some sites have total area equal to zero - this is ",
         "not possible with the logarithmic model")
    } else {
      warning("Some sites have total area equal to zero")
    }
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
    if (!any(c("partial", "exhaustive", "none") %in% gridStart)){
      stop("gridStart should be one of 'none', 'partial' or 'exhaustive")
    }
  }#eo is.null(startPar)
  
  ##Need to then fit the models for each SR
  k <- 1
  res <- lapply((CN + 1):ncol(data), function(x){
    dum <- data[,c(1:CN, x)]
    dum <- dum[order(dum[,ncol(dum)]),]
    colnames(dum)[ncol(dum)] <- "S"
    if (!is.null(startPar)){
      startPar2 <- startPar[k,]
    } else {
      startPar2 <- startPar
    }
    
    CO <- countryside_optim(dat = dum,
                      modType = modType,
                      gridStart = gridStart,
                      startPar = startPar2,
                      zLower = zLower)
    k <<- k + 1
    CO
  })
  
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
  res[[8]] <- list(habNam, spNam)
  names(res) <- c("fits", "affinity", "c", "Pred.Tot.Rich",
                  "rss", "data", "pow.model", "Group.Names")
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
#'   values used in the original model fits.
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
#' gridStart = "partial", habNam = c("AG", "SH",
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
