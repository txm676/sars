###########Fit the GDM#############################

#' @export


gdm <- function(data){
  if (anyNA(data)) stop("NAs present in data")
  if (!(is.matrix(data) || is.data.frame(data))) stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)
  
}
  
  
  