#First, the two functions which are derived from Coleman et al. (1982), which fit the model #and determine the standard deviation around the model’s predicted values:

sa <- function(x, a){
  sa <- (1 - x) ^ a
  return(sa)
}

sa2 <- function(x, a){
  sa2 <- (1 - x)^(2 * a)
  return(sa2)
}
