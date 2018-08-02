#' Generate a non parametric confidence interval for a SAR

#' @description Generate a non parametric confidence interval for a SAR
#' @usage sar_conf_int(obj, n_boot = 99)
#' @param obj An object of class 'sar' or 'multi.sar'.
#' @param n_boot the number of bootstrap sample used to construct the confidence interval.
#' @return 
#' @examples
#' #multi_sar fit
#' data(galap)
#' fit <- sar_multi(data = galap)
#' conf_fit <- sar_conf_int(fit)
#' #simple sar fit
#' data(galap)
#' fit <- sar_power(data = galap)
#' conf_fit <- sar_conf_int(fit)

#note: with many SAR models in the multi fit, it can take a long time
#models are removed from confint (but not sar_multi) if NAs are produced in Jacobian or in the transResiduals

#' @export

#fit <- sar_multi(galap)
#s <- sar_conf_int(fit, n = 100)

#n = number of iterations

sar_conf_int <- function(fit, n, crit = "Info", normaTest = "lillie",
                         homoTest = "cor.fitted",
                         neg_check = TRUE,
                         alpha_normtest = 0.05,
                         alpha_homotest = 0.05){
  
  if (!"multi" %in% class(fit)) stop ("class of 'fit' should be 'multi'")
  if (length(fit$details$mod_names) < 2) stop ("less than two models in the sar multi object")
  
  #model names for matching
  
  x1 <- c("Power", "PowerR", "Extended_Power_model_1", "Extended_Power_model_2", "Persistence_function_1", 
          "Persistence_function_2", "Exponential", "Kobayashi", "MMF", "Monod", "Negative_exponential", 
          "Chapman_Richards", "Cumulative_Weibull_3_par.", "Asymptotic_regression", "Rational_function", 
          "Gompertz", "Cumulative_Weibull_4_par.", "Beta-P_cumulative", "Heleg(Logistic)", "Linear_model")
  
  x2 <- c("sar_power(", "sar_powerR(", "sar_epm1(", "sar_epm2(", "sar_p1(", "sar_p2(", 
          "sar_expo(", "sar_koba(", "sar_mmf(", "sar_monod(", "sar_negexpo(", 
          "sar_chapman(", "sar_weibull3(", "sar_asymp(", "sar_ratio(", "sar_gompertz(", 
          "sar_weibull4(", "sar_betap(", "sar_heleg(", "sar_linear(")
  
  x3 <-   c("power", "powerR","epm1","epm2","p1","p2","expo","koba","mmf","monod","negexpo","chapman",
            "weibull3","asymp","ratio","gompertz","weibull4","betap","heleg", "linear")
  
  
  
  #observed data
  dat <- fit$details$fits[[1]]$data
  
  #weights and model names
  wei <- fit$details$weights_ics
  nams <- as.vector(names(wei))
  ns1 <- which(x1 %in% nams)
  nams_short <- x3[ns1]#for use below
  
  #loop over all models in sar_multi fit and fill matrices of fitted values and, resids & transformed residuals 
  calculated <- residuals <- transResiduals <- matrix(nrow = length(nams), ncol = nrow(dat))
  
  for (i in seq_along(nams))
  {
  
  #select the expression of the selected model
  wn <- which(x1 %in% nams[i])
  w2 <- x2[wn]

    
  #fit the best model to observed data; extract fitted values and residuals
  me <- eval(parse(text = paste(w2, "dat)", sep = "")))
  meF <- me$calculated
  meR <- me$residuals
  
  if (nams[i] == "Linear_model"){
  #based on code from mmSAR in Rforge
   jacob <- numDeriv::jacobian(me$model$rss.fun, me$par, data = me$data[1, ], model = me$model,  opt = FALSE)
  
   for (k in 2:nrow(dat)) {
     jacob <- rbind(jacob, numDeriv::jacobian(me$model$rss.fun, me$par, data = me$data[k, ],model = me$model, opt = FALSE)) 
  }
  } else {
  #based on code from mmSAR in Rforge
  jacob <- numDeriv::jacobian(me$model$rss.fun, me$par, data = me$data[1, ], opt = FALSE)
  
  for (k in 2:nrow(dat)) {
    jacob <- rbind(jacob, numDeriv::jacobian(me$model$rss.fun, me$par, data = me$data[k, ], opt = FALSE)) 
  }
  }#eo if linear
  
  ##occasionally the jacobian function returns NA for certain models: need to remove these,
  #so just fill all NAs for this model and it will be removed below
  if (anyNA(jacob)){
    calculated[i, ] <- NA
    residuals[i, ] <- NA
    transResiduals[i, ] <- NA
    next
  }
  
  jacobbis <- t(jacob) %*% jacob
  s <- svd(jacobbis)
  jacobbismun <- s$v %*% (diag(1 / s$d)) %*% (t(s$u))
  hatMat <- jacob %*% jacobbismun %*% t(jacob)
  matList <- list(jacob = jacob, hatMat = hatMat)
  
  #Residuals transformation from Davidson and Hinkley, 1997 "Bootstrap methods and their applications" p 259 eq (6.9)
  diagHatMat <- diag(hatMat)
  tr <- meR - mean(meR)
  tr <- tr / sqrt(1 - diagHatMat)#checked and this gives same values as mmSAR
  
  #fill in the matrices for this model
  calculated[i, ] <- meF
  residuals[i, ] <- meR
  transResiduals[i, ] <- tr

  }#eo for
  
  rownames(calculated) <- rownames(residuals) <- rownames(transResiduals) <- nams
  
  
  #some models have NANs in the transResiduals calculation; so remove them from all the matrices (and names vectors)
 if (anyNA(transResiduals)){
   
  wna <- which(apply(transResiduals, 1, anyNA))
  warning(paste("The following model(s) has been removed from the confidence interval calculations as NAs are 
                present in the transformed residuals:", names(wna)))
  calculated <- calculated[-wna, ]
  residuals <- residuals[-wna, ]
  transResiduals <- transResiduals[-wna, ]
  wei <- wei[-wna]
  nams <- nams[-wna]
  nams_short <- nams_short[-wna]
  if (identical(!nrow(calculated), nrow(residuals), 
                nrow(transResiduals), length(wei), length(nams), length(nams_short))) stop ("Problem with removed models 
                                                                                            following transResiduals checks")
 }

  #run the boostrapping (based on code in mmSAR)
  
  z <- 1
  nBoot <- n
  
  pointsNames <- paste("S",c(1:nrow(dat)))
  
  #Matrix of boot Samples
  bootMatrix=matrix(0, nBoot, nrow(dat))
  
  #list of optimisation results
  optimBootResult = vector("list", length = nBoot)
  
  #choosing an IC criterion (AIC or AICc or BIC): same code as within sar_multi as needs to be identical
  IC <- switch(crit,
               Info= if ( (nrow(dat) / 3) < 40 ) { "AICc" } else { "AIC" },
               Bayes= "BIC")
  

  #listof calculated values
  bootCalculated <- vector("list", length = nBoot)
  
  
  #test variable
  nGoodBoot = 1
  
  filtModelList <- nams
  

  while (nGoodBoot < nBoot+1) {
    
    test <- 1
    chousModel = filtModelList[rmultinom(1, 1, wei)==1]

    while (test != 0 ) {
      
      for (l in 1:nrow(dat)) {
        positives = transResiduals[chousModel,][transResiduals[chousModel,] > 0]
        negatives = transResiduals[chousModel,][transResiduals[chousModel,] < 0]
        
        vtci = negatives[abs(negatives) <= calculated[chousModel,l] ]
        vtci = c(vtci,positives)
        value = sample(vtci,1)
        
        bootMatrix[nGoodBoot,l] <- calculated[chousModel,l] + value
      }#end of for
      
      #test if one species richness is negative
      test = length( which(bootMatrix[nGoodBoot,]<0) )
      
    }#end of while

    #fit multiSAR
    badBoot=FALSE
  
    df <- data.frame("A" = dat$A, "S" = bootMatrix[nGoodBoot, ])
    optimres <- tryCatch(suppressWarnings(sar_multi(df, obj = nams_short, verb = FALSE)), error = function(e) NA)
 
    if (length(optimres) == 1) {
      badBoot=TRUE
    } else {
      if ((z %% 10) == 0) print(z)
      z <- z + 1
      #matrix for the fitted values
      bootCalculated[[nGoodBoot]] <- matrix(nrow = length(optimres$details$fits), ncol = nrow(df))
      rownames(bootCalculated[[nGoodBoot]]) <- as.vector(optimres$details$mod_names)
      #matrix for AIC etc
      optimBootResult[[nGoodBoot]] <- matrix(nrow = length(optimres$details$fits), ncol = 3)
      rownames( optimBootResult[[nGoodBoot]]) <- as.vector(optimres$details$mod_names)
      colnames(optimBootResult[[nGoodBoot]]) <- c(paste(IC), "Delta", "Weights")
      for (k in 1:length(optimres$details$mod_names)){
      bootCalculated[[nGoodBoot]][k, ] <- optimres$details$fits[[k]]$calculated
      optimBootResult[[nGoodBoot]][k, 1] <- as.vector(unlist(optimres$details$fits[[k]][IC]))
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
      for (i in 1:nrow(dat)) {
        bootHat[k,i] <- sum(akaikeweightvect * bootCalculated[[k]][,i])
      }#end of for i
      
    } else { bootHat[k,] <- rep(0, nrow(dat)) }
    
  } #end of for k
  
  #sort and return the CIs
  bootSort <- apply(bootHat, 2, sort)
  alp <- 0.025
  c1 <- ceiling(n  * alp) #maybe a better way
  c2 <- floor(n  * (1-alp))
  CI <- data.frame("L" = bootSort[c1, ], "U" = bootSort[c2, ])
  return(CI)
}