#First, the two functions which are derived from Coleman et al. (1982), which fit the model
#and determine the standard deviation around the model’s predicted values:

sa <- function(x, a){
  sa <- (1 - x) ^ a
  return(sa)
}

sa2 <- function(x, a){
  sa2 <- (1 - x)^(2 * a)
  return(sa2)
}


#' Fit Coleman's Random Placement Model
#'
#' @description Fit Coleman's (1981) random placement model to a species-site
#'   abundance matrix: rows are species and columns are sites. Note that the
#'   data must be abundance data and not presence-absence data. According to
#'   this model, the number of species occurring on an island depends on the
#'   relative area of the island and the regional relative species abundances.
#'   The fit of the random placement model can be determined through use of a
#'   diagnostic plot (see \code{\link{plot.coleman}}) of island area (log
#'   transformed) against species richness, alongside the model’s predicted
#'   values (see Wang et al., 2010). Following Wang et al. (2010), the model is
#'   rejected if more than a third of the observed data points fall beyond one
#'   standard deviation from the expected curve.
#' @usage coleman(data, area)
#' @param data A dataframe or matrix in which rows are species and columns are
#'   sites. Each element/value in the matrix is the abundance of a given species
#'   in a given site.
#' @param area A vector of site (island) area values. The order of the vector
#'   must match the order of the columns in \code{data}.
#' @import stats
#' @return A list of class "coleman" with four elements. The first element
#'   contains the fitted values of the model. The second element contains the
#'   standard deviations of the fitted values, and the third and fourth contain
#'   the relative island areas and observed richness values, respectively.
#'   \code{\link{plot.coleman}} plots the model.
#' @references Coleman, B. D. (1981). On random placement and species-area
#'   relations. Mathematical Biosciences, 54, 191-215.
#'
#'   Matthews, T. J., Cottee-Jones, H. E. W., & Whittaker, R. J. (2015).
#'   Quantifying and interpreting nestedness in habitat islands: a synthetic
#'   analysis of multiple datasets. Diversity and Distributions, 21, 392-404.
#'
#'   Wang, Y., Bao, Y., Yu, M., Xu, G., & Ding, P. (2010). Nestedness for
#'   different reasons: the distributions of birds, lizards and small mammals on
#'   islands of an inundated lake. Diversity and Distributions, 16, 862-873.
#' @examples
#' data(cole_sim)
#' fit <- coleman(cole_sim[[1]], cole_sim[[2]])
#' plot(fit, ModTitle = "Hetfield")
#' @export


coleman <- function(data, area){
  if (any(area <= 0)) stop("Area value <=0")
  if (anyNA(data) | anyNA(area)) stop("NAs present in data")
  if (!(is.matrix(data) | is.data.frame(data))) stop("data must be a matrix or dataframe")
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
  s_alp <- vector(length = length(ra))
  s_hat <- vector(length = length(ra))
  for (j in seq_along(ra)){
    for (i in seq_along(abun)){
      s_alp[i] <- sa(ra[j],abun[i])
    } #eo i
    s_hat[j] <- tot_sp - sum(s_alp)
  } #eo j

  ##Obtain the model variance
  varz <- vector(length = length(ra))
  for (j in seq_along(ra)){
    s_alp2 <- vector(length = length(abun))
    s_alp3 <- vector(length = length(abun))
    for (i in seq_along(abun)){
      s_alp2[i] <- sa(ra[j],abun[i])
      s_alp3[i] <- sa2(ra[j],abun[i])
    } #eo i
    varz[j] <-(sum(s_alp2)) - (sum(s_alp3))
  } #eo j

  ##Standard deviation
  sd <- sqrt(varz)
  plus_sd <- s_hat + sd
  min_sd <- s_hat - sd

res <- list("Predicted_values" = s_hat, "Standard_deviation" = sd,
            "Relative_areas" = ra, "Species_richness" = ob_sp)
class(res) <- "coleman"
return(res)
}
