###########Fit the GDM#############################

#' @importFrom stringr str_split
#' @export

##decide on whether to use dredge or not, or both
##work out how to add in linear, power and expo models
##allow user to fit and compare power with linear and expo
##print and plot functions


data <- data.frame("A" = c(10,40,80,160,160), "S" = c(1,3,5,8,10), Ti = c(1,2,3,4,5))



gdm <- function(data, model = "lin_pow", mod_sel = F, A = 1, S = 2, Ti = 3){
  if (anyNA(data)) stop("NAs present in data")
  if (!(is.matrix(data) || is.data.frame(data))) stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)
  if (ncol(data) < 3) stop("Not enough columns/variables to fit GDM")
  if (ncol(data) > 3) {
    warning("More than three columns in dataframe: using the first three")
    data <- data[, 1:3]
  }
  if (!(model %in% c("lin_pow", "power", "expo"))) {
    stop("provided model name not available")
  }
  
  if (A == 1 && S == 2 && Ti == 3){
    colnames(data) <- c("A", "S", "Time")
  } else{
    data <- data[, c(A, S, Ti)]
    colnames(data) <- c("A", "S", "Time")
  }
  
  data$Time2 <- data$Time ^ 2
  
  if (model == "lin_pow"){
    cat("\n","Fitting the GDM using the linear (log-log) power model", "\n")
    if (any(data$S == 0)) data$S <- data$S + 0.1
    data$A <- log(data$A)
    data$S <- log(data$S)
    fit <- lm(S ~ A + Time + Time2, data = data)
    
    if (mod_sel == F) class(fit) <- c("sars.gdm", "lm")

    if (mod_sel == T){
      #maybe just do MuMIn dredge instead (issue with na.omit though)
      #MuMIn::dredge(fit)
      fitL <- list()
      fitL[[1]] <- fit
      fitL[[2]] <- lm(S ~ A + Time, data = data)
      fitL[[3]] <- lm(S ~ Time, data = data)
      fitL[[4]] <- lm(S ~ A, data = data)
      fitL[[5]] <- lm(S ~ 1, data = data)
      fit <- fitL
    }
    
  }
  
  if (model == "power"){
    cat("\n","Fitting the GDM using the non-linear power model", "\n")

    ss <- sar_power(data = data[ ,1:2])

    #build power funct.
    ex <- ss$model$exp
    uex <-  unlist(str_split(ex, "expression"))
    newF <- paste("S ~ ", uex, " + j * Time + k * Time2", sep = "")
    nsp <- names(ss$par)

    sdf <- as.data.frame(matrix(ncol = length(nsp)))
    colnames(sdf) <- nsp
    sdf$c <- 1
    sdf$z <- 1
    sdf$j <- 1
    sdf$k <- 0
    
     fit <- tryCatch(nls(formula = S ~ c * A^z + j * Time + k * Time2, data = data, start = sdf),
                         error = function(e){e})
  
     try(nls(S ~ exp(c + z*log(A) + x*Time + y*Time2),data = data, start = data.frame(c=0, z=1, x=1, y=0)))

  }
  
  
  return(fit)
}
  
