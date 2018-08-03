#' A SAR dataset describing the plants of the Galapagos Islands
#'
#' A sample dataset in the correct mmSAR2 format: contains the areas of a number
#' of islands in the Galapagos, and the number of plant species recorded on each
#' island.
#'
#' @usage data(galap)
#' @format A list with two elements. The first element contains the name of the
#'   dataset. The second element contains a data frame with 2 columns and 16
#'   rows. Each row contains the area of an island (km2) in the Galapagos (1st
#'   column) and the number of plants on that island (2nd column).Preston (1962)
#'   also includes the island of Albemarle, but we have excluded this as it is
#'   almost six times larger than the second largest island.
#' @source Preston FW 1962. The Canonical Distribution of Commonness and Rarity:
#'   Part I. – Ecology 43:185-215.
#' @examples
#' data(galap)
"galap"


