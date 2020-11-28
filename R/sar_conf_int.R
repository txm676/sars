#### INTERNAL FUNCTION(S)

#' @importFrom numDeriv jacobian
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats rmultinom

sar_conf_int <- function(fit, n, crit, obj_all, normaTest,
                         homoTest,
                         homoCor,
                         neg_check,
                         alpha_normtest,
                         alpha_homotest, 
                         grid_start,
                         grid_n,
                         verb){

  if (!"multi" %in% class(fit)) stop ("class of 'fit' should be 'multi'")
  if (length(fit$details$mod_names) < 2) 
    stop ("less than two models in the sar average object")

  #model names for matching
  
  #NB lots of the hashed code relates to changes to do with grid_start being more
  #fully incorporated into the package, and having to chat how this works

#  x1 <- c("Power", "PowerR", "Extended_Power_model_1", 
          #"Extended_Power_model_2", "Persistence_function_1",
         # "Persistence_function_2", "Logarithmic", "Kobayashi", "MMF", 
         # "Monod", "Negative_exponential", "Chapman_Richards", 
        #  "Cumulative_Weibull_3_par.", "Asymptotic_regression", 
        #  "Rational_function","Gompertz", "Cumulative_Weibull_4_par.", 
        #  "Beta-P_cumulative", "Heleg(Logistic)", "Linear_model")

 # x2 <- c("sar_power(", "sar_powerR(", "sar_epm1(", "sar_epm2(", "sar_p1(",
      #    "sar_p2(", "sar_loga(", "sar_koba(", "sar_mmf(", "sar_monod(",
      #    "sar_negexpo(", "sar_chapman(", "sar_weibull3(", "sar_asymp(", 
      #    "sar_ratio(", "sar_gompertz(", "sar_weibull4(", "sar_betap(", 
       #   "sar_heleg(", "sar_linear(")

#  x3 <-   c("power", "powerR","epm1","epm2","p1","p2","loga","koba","mmf",
  #          "monod","negexpo","chapman",
   #         "weibull3","asymp","ratio","gompertz","weibull4","betap","heleg",
    #        "linear")


  
  #observed data
  dat <- fit$details$fits[[1L]]$data

  #weights and model names
  wei <- fit$details$weights_ics
  nams <- as.vector(names(wei))
#  ns1 <- which(x3 %in% nams)
 # nams_short <- x3[ns1]#for use below
  
  #loop over all models in sar_average fit and fill matrices of fitted values
  #and, resids & transformed residuals
  calculated <- residuals <- transResiduals <- matrix(nrow = length(nams),
                                                      ncol = nrow(dat))

  for (i in seq_along(nams))
  {

  #select the expression of the selected model
  #wn <- which(x3 %in% nams[i])
 # w2 <- x2[wn]

  #fit the ith model to observed data; extract fitted values and residuals
    
  #me <- suppressWarnings(eval(parse(text = paste(w2, "dat)", sep = ""))))
  #original code did not account for grid_start, when this is used you can get 
  #different parameter values from fitting the model without this, so instead 
  #of re-fitting, just take the fits from the fit object
    
  me <- fit$details$fits[nams[i]][[1]]
  meF <- me$calculated
  meR <- me$residuals

  if (nams[i] == "linear"){
  #based on code from mmSAR in Rforge
   jacob <- suppressWarnings(numDeriv::jacobian(me$model$rss.fun, me$par, 
                              data = me$data[1, ],model = me$model,  
                              opt = FALSE))

   for (k in 2:nrow(dat)) {
     jacob <- rbind(jacob, suppressWarnings(numDeriv::jacobian(me$model$rss.fun,
                                             me$par, data = me$data[k, ],
                                             model = me$model, opt = FALSE)))
  }
  } else {
  #based on code from mmSAR in Rforge
  jacob <- suppressWarnings(numDeriv::jacobian(me$model$rss.fun, me$par, 
                                    data = me$data[1, ], opt = FALSE))

  for (k in 2:nrow(dat)) {
    jacob <- rbind(jacob, suppressWarnings(numDeriv::jacobian(me$model$rss.fun, 
                                           me$par, data = me$data[k, ], 
                                           opt = FALSE)))
  }
  }#eo if linear

  ##occasionally the jacobian function returns NA for certain models: 
  #need to remove these,
  #looks like it happens when the slope parameter (e.g. z) is very close to
  #0, even if still > 0 (e.g. 1.1e-17) as looking in the source code jacobian 
  #seems to have a zero sensitivity of 1e-4
  #so just fill all NAs for this model and it will be removed below
  if (anyNA(jacob)){
    calculated[i, ] <- NA
    residuals[i, ] <- NA
    transResiduals[i, ] <- NA
    next
  }

  jacobbis <- t(jacob) %*% jacob

  ##occasionally the jaccobis function returns Inf, which errors below,
  #so need to removed this model
  if (anyNA(jacobbis)  | any(jacobbis == "Inf")){
    calculated[i, ] <- NA
    residuals[i, ] <- NA
    transResiduals[i, ] <- NA
    next
  }

  s <- svd(jacobbis)
  jacobbismun <- s$v %*% (diag(1 / s$d)) %*% (t(s$u))
  hatMat <- jacob %*% jacobbismun %*% t(jacob)
  matList <- list(jacob = jacob, hatMat = hatMat)

  #Residuals transformation from Davidson and Hinkley, 1997 "Bootstrap
  #methods and their applications" p 259 eq (6.9)
  diagHatMat <- diag(hatMat)
  tr <- meR - mean(meR)
  #checked and this gives same values as mmSAR
  tr <- suppressWarnings(tr / sqrt(1 - diagHatMat))

  #fill in the matrices for this model
  calculated[i, ] <- meF
  residuals[i, ] <- meR
  transResiduals[i, ] <- tr
  }#eo for

 rownames(transResiduals) <- nams
 rownames(residuals) <- rownames(transResiduals)
 rownames(calculated) <- rownames(residuals)

  #some models have NANs in the transResiduals calculation; so remove them 
 #from all the matrices (and names vectors)
 if (anyNA(transResiduals)){

  wna <- which(apply(transResiduals, 1, anyNA))
  warning(paste("The following model(s) has been removed from the", 
                 " confidence interval calculations as NAs/Infs are
                present in the transformed residuals:", 
                paste(names(wna), collapse = ", ")))
  calculated <- calculated[-wna, ]
  residuals <- residuals[-wna, ]
  transResiduals <- transResiduals[-wna, ]
  wei <- wei[-wna]
  nams <- nams[-wna]
#  nams_short <- nams_short[-wna]
  #identical only works with two objects, so use sapply to compare all objects
  #in a pairwise fashion
  if (!all(sapply(list(nrow(residuals),
                    nrow(transResiduals), length(wei),
                    length(nams)), 
               FUN = identical, nrow(calculated))))
    stop ("Problem with removed models following transResiduals checks")
 }

  #run the bootstrapping (based on code in mmSAR)

  if (verb) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    z <- 1
  }
  nBoot <- n

  pointsNames <- paste("S",c(seq_len(nrow(dat))))

  #Matrix of boot Samples
  bootMatrix <- matrix(0, nBoot, nrow(dat))

  #list of optimisation results
  optimBootResult <- vector("list", length = nBoot)

  #choosing an IC criterion (AIC or AICc or BIC): same code as within
  #sar_average as needs to be identical
  IC <- crit

  #listof calculated values
  bootCalculated <- vector("list", length = nBoot)


  #test variable
  nGoodBoot <- 1

  filtModelList <- nams


  while (nGoodBoot < nBoot+1) {

    test <- 1
    chousModel <- filtModelList[rmultinom(1, 1, wei)==1]

    while (test != 0 ) {

for (l in seq_len(nrow(dat))) {
  positives <- transResiduals[chousModel, ][transResiduals[chousModel, ] > 0]
  negatives <- transResiduals[chousModel, ][transResiduals[chousModel, ] < 0]
  #I think this line ensures that there can be no negative S values, but then
  #I am not sure why the test is needed? Check with FG
  vtci <- negatives[abs(negatives) <= calculated[chousModel, l] ]
  vtci <- c(vtci, positives)
  value <- sample(vtci, 1)
  bootMatrix[nGoodBoot, l] <- calculated[chousModel, l] + value
}#end of for

      #test if one species richness is negative
      test <- length( which(bootMatrix[nGoodBoot, ] < 0) )

    }#end of while

    #fit multiSAR
    badBoot <- FALSE

    df <- data.frame("A" = dat$A, "S" = bootMatrix[nGoodBoot, ])
    
    #in here we use for model names obj_all which means that all the models
    #that the user originally tried to fit are used, even if some of them were
    #removed from the original mmi curve due to failing residuals tests etc.
    #It also now uses the same arguments for checks and grid_start etc that
    #the user originally uses with their mmi curve
    optimres <- tryCatch(suppressMessages(sar_average(obj = obj_all,
                                                data = df,
                                                crit = crit, 
                                                normaTest = normaTest,
                                                homoTest = homoTest,
                                                homoCor = homoCor,
                                                neg_check = neg_check,
                                                alpha_normtest = alpha_normtest,
                                                alpha_homotest = alpha_homotest, 
                                                grid_start = grid_start,
                                                grid_n = grid_n,
                                      verb = FALSE)), error = function(e) NA)

    if (length(optimres) == 1) {
      badBoot <- TRUE
    } else {
      #progress bar
      if (verb) {
        setTxtProgressBar(pb, z)
        z <- z + 1
      }
      #matrix for the fitted values
      bootCalculated[[nGoodBoot]] <- 
        matrix(nrow = length(optimres$details$fits), ncol = nrow(df))
      rownames(bootCalculated[[nGoodBoot]]) <- 
        as.vector(optimres$details$mod_names)
      #matrix for AIC etc
      optimBootResult[[nGoodBoot]] <- 
        matrix(nrow = length(optimres$details$fits), ncol = 3)
      rownames( optimBootResult[[nGoodBoot]]) <- 
        as.vector(optimres$details$mod_names)
      colnames(optimBootResult[[nGoodBoot]]) <- 
        c(paste(IC), "Delta", "Weights")

      for (k in seq_along(optimres$details$mod_names)){
      bootCalculated[[nGoodBoot]][k, ] <- 
        optimres$details$fits[[k]]$calculated
      optimBootResult[[nGoodBoot]][k, 1] <- 
        as.vector(unlist(optimres$details$fits[[k]][IC]))
      optimBootResult[[nGoodBoot]][k, 2:3] <- 0
      }#eo k
    }#eo if

    if (badBoot == FALSE) {
      nGoodBoot <- nGoodBoot + 1
    }#end of if
  }#eo top while



  bootHat <- matrix(0, nBoot, nrow(dat))
  f <- 0

  for (k in 1:nBoot) {

    DeltaICvect <- vector()
    akaikeweightvect <- vector()
    filtNlig <- dim(optimBootResult[[k]])[1]

    if (filtNlig != 0) {

      for (i in 1:filtNlig){

        #Delta IC = ICi - ICmin
      DeltaIC <- optimBootResult[[k]][i, 1] - min(optimBootResult[[k]][, 1])
      DeltaICvect <- c(DeltaICvect, DeltaIC)
      }#end of for i

      for (i in 1:filtNlig){
        #Akaike Weigths
        akaikesum <- sum(exp( -0.5*(DeltaICvect)))
        akaikeweight <- exp(-0.5*DeltaICvect[i]) / akaikesum
        akaikeweightvect <- c(akaikeweightvect, akaikeweight)
      }#end of for i

      optimBootResult[[k]][, 2] <- DeltaICvect
      optimBootResult[[k]][, 3] <- akaikeweightvect

      #Averaging
      for (i in seq_len(nrow(dat))) {
        bootHat[k,i] <- sum(akaikeweightvect * bootCalculated[[k]][,i])
      }#end of for i

    } else { bootHat[k,] <- rep(0, nrow(dat)) }

  } #end of for k

  if (verb) cat("\n")#needed to drop warnings below progress bar

  #sort and return the CIs
  bootSort <- apply(bootHat, 2, sort)
  alp <- 0.025
  c1 <- ceiling(n  * alp) #maybe a better way
  c2 <- ceiling(n  * (1-alp))
  CI <- data.frame("L" = bootSort[c1, ], "U" = bootSort[c2, ])
  return(CI)
}
