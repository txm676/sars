###########Fit the GDM#############################

#data <- data.frame("A" = c(10,40,80,160,160), "S" = c(1,3,5,8,10), Ti = c(1,2,3,4,5))

#' @export

#no R2 provided for non-linear models as only for linear models, residual standard
#error provided

#AIC calculated using AIC function (log likelihood) not the 
#RSS like rest of sars: check with FG this is OK

#all model comparison does not include log-log power as can't use AIC

#3d plotting will be provided in a future version of the package

gdm <- function(data, model = "lin_pow", mod_sel = FALSE, A = 1, S = 2, Ti = 3){
  if (anyNA(data)) stop("NAs present in data")
  if (!(is.matrix(data) || is.data.frame(data))) stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)
  if (ncol(data) < 3) stop("Not enough columns/variables to fit GDM")
  if (ncol(data) > 3) {
    warning("More than three columns in dataframe: using the first three")
    data <- data[, 1:3]
  }
  if (!(model %in% c("lin_pow", "power", "expo", "linear", "all"))) {
    stop("provided model name not available")
  }
  if (!is.logical(mod_sel)) stop("mod_sel argument should be TRUE or FALSE")
  
  if (A == 1 && S == 2 && Ti == 3){
    colnames(data) <-c("Area", "SR", "Time")
  } else{
    data <- data[, c(A, S, Ti)]
    colnames(data) <- c("Area", "SR", "Time")
  }
  
  if (model == "all") allMods <- vector("list", length = 2)

  if (model == "lin_pow"){
    #cat("\n","Fitting the GDM using the linear (log-log) power model", "\n")
    if (any(data$S == 0)) data$S <- data$S + 0.1
    data$Area <- log(data$Area)
    data$SR <- log(data$SR)
    data$Time2 <- data$Time ^ 2
    fit <- lm(SR ~ Area + Time + Time2, data = data)
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 4)
      fitL[[1]] <- fit
      fitL[[2]] <- lm(SR ~ Area + Time, data = data)
      fitL[[3]] <- lm(SR ~ Area, data = data)
      fitL[[4]] <- lm(SR ~ 1, data = data)
      fit <- fitL
    }
    class(fit) <- c("gdm", "lm")
    attr(fit, "Type") <- "lin_pow"
    attr(fit, "mod_sel") <- mod_sel
    
  } 
  if (model == "expo" || model == "all"){
      
    #cat("\n","Fitting the GDM using the exponential model", "\n")

     fit <- nls(formula = SR ~ Int + A * log(Area) + Ti * Time + Ti2 * Time ^ 2, 
                         data = data, start = data.frame(Int = 0, A = 1, Ti = 1, Ti2 = 0))
     if (mod_sel == TRUE){
       fitL <- vector("list", length = 3)
       fitL[[1]] <- fit
       fitL[[2]] <- nls(formula = SR ~ Int + A * log(Area) + Ti * Time, 
                        data = data, start = data.frame(Int = 0, A = 1, Ti = 1))
       fitL[[3]] <- nls(formula = SR ~ Int + A * log(Area), 
                        data = data, start = data.frame(Int = 0, A = 1))
      # fitL[[4]] <- nls(formula = SR, data = data)
       fit <- fitL
     }
     class(fit) <- c("gdm", "nls")
     attr(fit, "Type") <- "expo"
     attr(fit, "mod_sel") <- mod_sel
     if (model == "all") allMods[[1]] <- fit
     
  } 
  if (model == "linear" || model == "all"){
   # cat("\n","Fitting the GDM using the linear model", "\n")
    fit <- nls(SR ~ Int + A * Area + Ti * Time + Ti2 * Time ^ 2, 
               data = data, start = data.frame(Int = 0, A = 1, Ti = 1, Ti2 = 0))
    
    if (mod_sel == TRUE){
      fitL <- vector("list", length = 3)
      fitL[[1]] <- fit
      fitL[[2]] <- nls(formula = SR ~ Int + A * Area + Ti * Time, 
                       data = data, start = data.frame(Int = 0, A = 1, Ti = 1))
      fitL[[3]] <- nls(formula = SR ~ Int + A * Area, 
                       data = data, start = data.frame(Int = 0, A = 1))
      # fitL[[4]] <- nls(formula = SR, data = data)
      fit <- fitL
    }
    class(fit) <- c("gdm", "nls")
    attr(fit, "Type") <- "linear"
    attr(fit, "mod_sel") <- mod_sel
    if (model == "all") allMods[[2]] <- fit
    }
  if (model == "all") {
    class(allMods) <- "gdm"
    attr(allMods, "Type") <- "allMods"
    attr(allMods, "mod_sel") <- mod_sel
    return(allMods)
  } else {
    return(fit)
  }
}
  

