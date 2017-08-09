


#The main function, which takes a species-site abundance matrix and produces the diagnostic #plot:

#' @export


coleman <- function(data, area){
  if (any(area <= 0)) stop("Area value <=0")
  if (anyNA(data) || anyNA(area)) stop("NAs present in data")
  if (!(is.matrix(data) || is.data.frame(data))) stop("data must be a matrix or dataframe")
  if (is.matrix(data)) data <- as.data.frame(data)

  ##check each species has 1 or more individuals/each sites as at least one species present
  if (any(colSums(data) == 0)) warning("Matrix contains sites with no species")
  if (any(rowSums(data) == 0)) warning("Matrix contains species which were not sampled")

  tot_sp <- nrow(data) #Number of species
  tot_area <- sum(area) #Total area

  ####Derive the observed species richness for each column (site)

  ob_sp <- colSums(apply(data, 2, function(x) ifelse(x > 0, 1, 0)))

  ra <- area/tot_area #relative areas
  abun <- rowSums(data) #abundance of each species

  ##obtain predicted values from random placement model
  s_alp <- c()
  s_hat <- c()
  for (j in seq_along(ra)){
    for (i in seq_along(abun)){
      s_alp[i] <- sa(ra[j],abun[i])
    } #eo i
    s_hat[j] <- tot_sp - sum(s_alp)
  } #eo j

  ##Obtain the model variance
  varz <- c()
  for (j in seq_along(ra)){
    s_alp2 <- c()
    s_alp3 <- c()
    for (i in seq_along(abun)){
      s_alp2[i] <- sa(ra[j],abun[i])
      s_alp3[i] <- sa2(ra[j],abun[i])
    } #eo i
    varz[j] <-(sum(s_alp2)) - (sum(s_alp3))
  } #eo j

  ##Standard deviation
  sd <- sqrt (varz)
  plus_sd <- s_hat + sd
  min_sd <- s_hat - sd

res <- list("Predicted_values" = s_hat, "Standard_deviation" = sd,
            "Relative_areas" = ra, "Species_richness" = ob_sp)
class(res) <- "coleman"
return(res)
}


