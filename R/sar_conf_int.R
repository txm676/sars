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
                         verb,
                         display){
  
  if (!"multi" %in% class(fit)) stop ("class of 'fit' should be 'multi'")
  if (length(fit$details$mod_names) < 2) 
    stop ("less than two models in the sar average object")

  #observed data
  dat <- fit$details$fits[[1L]]$data
  
  #weights and model names
  wei <- fit$details$weights_ics
  nams <- as.vector(names(wei))
  
  #loop over all models in sar_average fit and fill matrices of fitted values
  #and, resids & transformed residuals
  calculated <- residuals <- transResiduals <- matrix(nrow = length(nams),
                                                      ncol = nrow(dat))
  
  for (i in seq_along(nams))
  {

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
  
  if (display) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
    z <- 1
  }
  nBoot <- n
  
  #choosing an IC criterion (AIC or AICc or BIC): same code as within
  #sar_average as needs to be identical
  IC <- crit
  
  #test variable
  nGoodBoot <- 1
  
  #matrix to put the final mmi values in
  bootHat <- matrix(0, nBoot, nrow(dat))
  
  filtModelList <- nams
  
  while (nGoodBoot < nBoot+1) {
    
    #create vector to fill with the bootstrapped values to use ad richness
    bootVec <- vector(length = nrow(dat))
    
    test <- 1
    chousModel <- filtModelList[rmultinom(1, 1, wei)==1]
    
    while (test != 0 ) {
      
      for (l in seq_len(nrow(dat))) {
        positives <- transResiduals[chousModel, ][transResiduals[chousModel, ] > 0]
        negatives <- transResiduals[chousModel, ][transResiduals[chousModel, ] < 0]
        #this line ensures that there can be no negative S values resulting from 
        #adding negative trans residuals to them
        vtci <- negatives[abs(negatives) <= calculated[chousModel, l] ]
        vtci <- c(vtci, positives)
        value <- sample(vtci, 1)
        bootVec[l] <- calculated[chousModel, l] + value
      }#end of for
      
      #test if one species richness is negative. Although we filter out any
      #negative trans residuals that would make the calculated value < 0 when
      #added to it, you can still get negative data values if the original
      #calculated value is < 0 and the selected positive trans residual is <
      #abs(calculated value). mmsAR originally dealt with this by only selecting
      #positive trans residuals which were > abs(calculated vale), but if no
      #positive trans residuals are > it just errors. This way here the while
      #loop just restarts and should (eventually) select a different model.
      test <- length(which(bootVec < 0))
      
    }#end of while
    
    #fit multiSAR
    df <- data.frame("A" = dat$A, "S" = bootVec)
    
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
                                                      verb = FALSE,
                                                      display = FALSE)), 
                         error = function(e) NA)
    
    #Nov 2020: replaced lots of code that was just calculating the mmi curve,
    #with the mmi curve values from sar_average. Checked and this produces
    #the same results.
    if (length(optimres) == 1) {
      next
    } else {
      #progress bar
      if (display) {
        setTxtProgressBar(pb, z)
        z <- z + 1
      }
      bootHat[nGoodBoot,] <- optimres$mmi
      
      nGoodBoot <- nGoodBoot + 1
    } 
  }#eo while
  
  if (display) cat("\n")#needed to drop warnings below progress bar
  
  #sort and return the CIs
  bootSort <- apply(bootHat, 2, sort)
  alp <- 0.025
  c1 <- ceiling(n  * alp)#ceiling avoids any cases of 0
  c2 <- ceiling(n  * (1-alp))
  CI <- data.frame("L" = bootSort[c1, ], "U" = bootSort[c2, ])
  return(CI)
}