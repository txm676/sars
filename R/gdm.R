###########Fit the GDM#############################

#' @export

##decide on whether to use dredge or not, or both
##work out how to add in linear, power and expo models
##allow user to fit and compare power with linear and expo
##print and plot functions



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
    cat("\n","Fitting the GDM using the linear power model", "\n")
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
  
 # if (model == "power"){
   # cat("\n","Fitting the GDM using the non-linear power model", "\n")
     # fit <- nls(S ~ c * A ^ z + j * Time + k * Time2, data = data, 
              #   start = data.frame(c = 5, z = 0.25, j = 1, k = 0))
  #}
  return(fit)
}
  
