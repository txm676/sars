#### INTERNAL FUNCTION(S)

########observed shape#######################

#' @importFrom stats uniroot

obs_shape <- function(x){

    minMaxVal <- NA
    asymp <- FALSE
    model <- x$model
    data <- x$data
    pars <- x$par
    Areas <- seq(range(data$A)[1], range(data$A)[2], length.out = 9999)

    
    fun_d1 <- function(A, par, model) {
      nam <- names(par)
      for (i in seq_along(par)) {
        assign(nam[i], par[i])
      }
      eval(model$d1.exp)
    }
    

    fun_d2 <- function(A, par, model) {
      nam <- names(par)
      for (i in seq_along(par)) {
        assign(nam[i], par[i])
      }
      eval(model$d2.exp)
    }
    
    
    
    #function to detect sign changes and provide roots
    getRoots <- function(fun, d1 = TRUE, Areas, pars, model) {
        #values and sign of the function evaluated
        values <- fun(Areas, pars, model)
        signs <- sign(values)
        minMax <- NA
        sigCh <- vector()
        
        #check if the sign of a value is the same as the next value
        for (i in 1:(length(signs)-1)) {
          if (sign(signs[i]) != sign(signs[i + 1])) sigCh <- c(sigCh, i)
        }
          nMinMax <- length(sigCh)
          if (nMinMax != 0){
            #check whether is first derivative function??
            #dc <- as.list(match.call())
            #if (as.character(dc$fun)[3] == "d1.fun") {
            if (d1) {
              if (nMinMax > 1){
             warning("more than one minimum and/or maximum in the derivative,
                      check the model plot to asses whether the model fit
                      looks sensible")
                for (i in 1:nMinMax) {
                  sigBef <- signs[sigCh[i]]
                  sigAft <- signs[sigCh[i] + 1]
                  if (sigBef == -1 & sigAft == 1){
                     minMax[i] <- "minima"
                  }
                  if (sigBef == 1 & sigAft == -1){
                     minMax[i] <- "maxima"
                  }
              }#eo for i
            } else {
              sigBef <- signs[sigCh]
              sigAft <- signs[sigCh + 1]
              if (sigBef == -1 & sigAft == 1){
                minMax <- "minima"
                }
              if (sigBef == 1 & sigAft == -1){
                minMax <- "maxima"
                }
            }#eo if/else

          }#eo if

          #roots
          roots <- vapply(seq_along(sigCh), FUN = function(x){
                            uniroot(fun, c(Areas[sigCh[x]],
                                        Areas[sigCh[x] + 1]),
                                    par = pars, model = model)$root},
                            FUN.VALUE = double(1))
          res <- list(sigCh = sigCh, roots = roots, minMax = minMax)
          return(res)
        } else {
            res <- list(sigCh = NA, roots = NA, minMax = minMax)
          return(res)
        }
      }#eo getRoots

      #Possible fits
      possFits <- rep(0, 4)
      names(possFits) <- c("linear", "convex up", "convex do", "sigmoid")

      #if the fit is linear
      min.area <- min(data$A)
      max.area <- max(data$A)
      ts <- seq(0,1, .01)
      #dif <- vector()

      dif <- vapply(seq_along(ts), FUN = function(x){
        model$mod.fun((ts[x]  * min.area + (1 - ts[x]) * max.area), pars) -
        (ts[x] * model$mod.fun(min.area, pars) +
           (1 - ts[x]) * model$mod.fun(max.area, pars))
      }, FUN.VALUE = double(1))
      
      if (anyNA(dif)){
          fitShape <- paste0(model$shape, " - observed shape algorithm failed:", 
                           " observed shape set to theoretical shape", " (", model$shape, ")")
          asymptote <- model$asymp(pars)
          if (asymptote){
            if (asymptote > range(data$S)[1] & asymptote < range(data$S)[2])
              asymp <- TRUE
          }#eo if
          res <- list(asymp = asymp, fitShape = fitShape)
          return(res)
      } 

      #is linear?
      if (sum(abs(dif) < 0.001) == length(ts)) possFits <- c(1, 0, 0, 0)

      #is it convex upward/downward?
      if ((sum(abs(dif) < 0.001) != length(ts)) &
          (sum(dif <= 0) == length(dif))) possFits <- c(0,0,1,0)
      if ((sum(abs(dif) < 0.001) !=length(ts))  &
          (sum(dif >= 0) == length(dif))) possFits <- c(0,1,0,0)

      #is the asymptote reached?
      asymptote <- model$asymp(pars)
      if (asymptote){
        if (asymptote > range(data$S)[1] & asymptote < range(data$S)[2])
          asymp <- TRUE
      }#eo if

      roots.d1 <- tryCatch(getRoots(fun_d1, d1 = TRUE, Areas, pars, model),
                           error = function(e) list(sigCh = NA, roots = NA,
                                                    minMax = NA))
      if (model$shape %in% c("sigmoid", "convex/sigmoid")) {
        roots.d2 <- tryCatch(getRoots(fun_d2, d1 = FALSE, Areas, pars, model),
                             error = function(e) list(sigCh = NA,
                                                    roots = NA, minMax = NA))
      } else {
          roots.d2 <- list(sigCh = NA, roots = NA, minMax = NA)
          }

      #getting the min/max value from the first derivative analysis
      minMax <- roots.d1$minMax
      if (!is.na(roots.d1$minMax[1])) minMaxVal <-
        model$mod.fun(roots.d1$roots, pars)

      #is it sigmoid?
      if (!is.na(roots.d2$roots[1])) possFits <- c(0, 0, 0, 1)

      #if the model is sigmoid but the shape study failed then
      #we will said the fit to be sigmoid
      fm <- NULL
      if (model$shape %in% c("sigmoid", "convex/sigmoid") & sum(possFits) == 0) {
        fm <- 1
        possFits <- c(0, 0, 0, 1)
      }
      names(possFits) <- c("linear", "convex up", "convex down", "sigmoid")
      if (is.null(fm)){
        fitShape <- names(possFits[possFits == 1])
      } else {
        fitShape <- paste("sigmoid - observed shape algorithm failed: observed",
                     "shape set to theoretical shape (sigmoid)")
      }
      if (all(possFits == 0)){
        fitShape <- paste0(model$shape, " - observed shape algorithm failed:", 
                           " observed shape set to theoretical shape", " (", model$shape, ")")
      }
      res <- list(asymp = asymp, fitShape = fitShape)
      return(res)
}
