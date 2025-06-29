
##filter_points_in_expanding_circles
#A function that picks a random point and starts sampling increasing circles
#from those points, returning the points inside those circles.

#' @importFrom sf st_geometfry st_buffer st_intersection st_area st_intersects

filter_points_in_expanding_circles <- function(points_sf,
                                               radius_vector,
                                               convex_hull) {
  # Randomly select one point from points_sf
  selected_point <- points_sf[sample(1:nrow(points_sf), 1), ]

  # Initialize a list to store points within each circle
  points_within_circles <- list()

  # Loop through the radius_vector
  for (radius in radius_vector) {
    # Create a buffer (circle) around the selected point with the current radius
    circle <- st_geometry(st_buffer(selected_point, dist = radius))

    # Calculate the intersection of the circle with the convex_hull
    intersection <- st_intersection(circle, convex_hull)

    # Check the area of the intersection compared to the area of the circle
    circle_area <- as.numeric(st_area(circle)) #as.numeric to remove units

    intersection_area <- as.numeric(st_area(intersection)) #as.numeric to remove units

    # Stop if the majority of the circle is outside the convex_hull
    if (intersection_area / circle_area < 0.5) {
      break
    }

    # Find points_sf that intersect the current circle
    points_in_circle <- points_sf[st_intersects(points_sf,
                                                circle,
                                                sparse = FALSE), ]

    # Add the points to the list
    points_within_circles[[paste0("radius_", radius)]] <-
      list(points=points_in_circle,circle=circle)
  }

  # Return the list of points within circles
  return(points_within_circles)
}


##summarize_samples
#A function that takes a list of groups of sites as a sample feature collections
#which has a dataframe with the species occurences, the polygons corresponding
#to the geometric delimitation of the groups of sites, a raster of habitat
#types, a vector of habitat names, the list of species in each group (each group
#is a list of numerical values of the columns corresponding to the species in
#the dataframe from the sample fature), and a vector with the species group
#names. The first column in the sample feature collection should be the name of
#the individual sample site. It produces a table with the area of each habitat
#and the number of species in each group and the total number of species.

#' @importFrom terra crop vect mask freq merge
#' @importFrom sf st_drop_geometry

summarize_samples <- function(samples, polygons, 
                             habitat_raster,
                             species_groups, 
                             species_group_names)
{
  # Initialize an empty data frame for the results
  habitat_names <-  terra::unique(habitat_raster, 
                                  na.rm = F)[,1]
  
  results_df <- data.frame(matrix(ncol = length(habitat_names)+
                                    length(species_group_names)+2,
                                  nrow = 0))
  colnames(results_df) <- c(habitat_names,"Area_Total",
                            species_group_names,"Sp_Total")

  # Iterate over each area sample (i.e. group of sites)
  for (i in seq_along(samples)) {
    sample <- samples[[i]]
    polygon <- polygons[[i]]

    # Crop and mask the raster to the polygon extent
    habitat_cropped <- crop(habitat_raster, vect(polygon))
    habitat_masked <- mask(habitat_cropped, vect(polygon))

    # Calculate the area of each habitat type
    habitat_df <- freq(habitat_masked, bylayer=FALSE)
    #if NA is included in land-use raster, add it in
    if (anyNA(habitat_names)){
    habitat_df <- rbind(habitat_df,freq(habitat_masked,
                                        bylayer=FALSE,
                                        value=NA))
    }

    if (nrow(habitat_df) != length(habitat_names)){
    
    # Ensure all possible values are included
    all_values_df <- data.frame(value = habitat_names,
                                count = 0)

    # Merge with actual frequency data, replacing 0 where missing
    habitat_df <- merge(all_values_df, habitat_df,
                        by = "value", all.x = TRUE)

    # Fill NA counts with 0
    habitat_df$count <- ifelse(is.na(habitat_df$count.y),
                               0, habitat_df$count.y)

    # Drop unnecessary column
    habitat_df <- habitat_df[, c("value", "count")]

    }#if nrow habitat_df
    
    habitat_df$area <- habitat_df$count *
      res(habitat_raster)[1] * res(habitat_raster)[2]

    # Store the results
    results_df[i, seq_along(habitat_names)] <- habitat_df$area # store area of each habitat
    results_df[i, length(habitat_names)+1] <-
      sum(results_df[i, seq_along(habitat_names)]) # store total area

    # Subset species occurrences for each group
    total_species <- 0
    for (k in seq_along(species_groups))
    {
      group_species <- st_drop_geometry(sample)[, species_groups[[k]]+1 ]
      species_present <- colSums(group_species > 0)  # Count species occurrences
      num_species <- sum(species_present > 0)       # Number of species in the group
      results_df[i, length(habitat_names)+1+k] <- num_species #store species number in results
      total_species <- total_species + num_species
    }

    #store the total number of species by summing across 
    #species groups
    results_df[i, length(habitat_names) +
                 length(species_groups)+2] <- total_species
  }
  
  #Move Area_total to end
  results_dfa <- subset(results_df, 
                        select = -c(Area_Total,Sp_Total))
  results_dfb <- cbind(results_dfa, 
                       subset(results_df, 
                              select = c(Area_Total,Sp_Total)))
    
  return(results_dfb)
}


##create_squares
#Create an sf where each sampling point is associated with a sampling square
#centred in the point and with a width dx.

#' @importFrom sf st_as_sf st_buffer

create_squares <- function(points_sf, width) {
  # Ensure the input is a point sf object
  if (!inherits(points_sf, "sf") || !inherits(st_geometry(points_sf),
                                              "sfc_POINT")) {
    stop("Input must be an sf object with point geometries.")
  }
  # Calculate half-width (to shift the square corners)
  half_width <- width / 2

  squares_sf <- st_as_sf(st_buffer(points_sf,
                                   dist = half_width,
                                   endCapStyle = "SQUARE"))

  return(squares_sf)
}


##filter_points_in_clusters
#Filter points for type II/III SAR (fractals sampling scheme or sub-divisions),
#Points should have an associated sampling square as geometry feature instead of
#a point geometry.
#' @importFrom terra split
#' @importFrom sf st_geometry st_coordinates
#' @importFrom stats kmeans

filter_points_in_clusters <- function(points_sf, squares_sf,
                                      cluster_size_vector)
{
  npoints <- nrow(points_sf)
  n_clusters_vector = npoints%/%cluster_size_vector
  n_clusters_vector[n_clusters_vector==0] <- 1  # when whole landscape, npoints<cluster_size
  points_within_clusters<- list()
  for (n_clusters in n_clusters_vector)
  {
    if(n_clusters==npoints)
    {
      points_in_clusters <- split(points_sf, 1:npoints)
      clusters_convex_hulls <- split(st_geometry(squares_sf), 1:npoints)
    }
    else
    {
      # Extract coordinates for k-means
      coords <- st_coordinates(points_sf)

      # Perform k-means clustering
      kmeans_result <- kmeans(coords, centers = n_clusters)

      # split the points sf into a list of clusters of points
      points_in_clusters <- split(points_sf, kmeans_result$cluster)
      squares_in_clusters <- split(squares_sf, kmeans_result$cluster)

      # Merge all squares into a single geometry
      merged_geom_within_clusters <- lapply(squares_in_clusters, st_union)

      # Compute convex hull
      clusters_convex_hulls <- lapply(merged_geom_within_clusters,
                                      st_convex_hull)

      #clusters_convex_hulls <- lapply(points_in_clusters,
      #                                   function(points) {st_convex_hull(st_union(points))})
    }
    # Add the points to the list
    points_within_clusters[[paste0("size_",
            cluster_size_vector[n_clusters_vector==n_clusters])]] <-
      list(points=points_in_clusters,chulls=clusters_convex_hulls)
  }
  # Return the list of points within clusters
  return(points_within_clusters)
}


#' @importFrom stats kmeans
##sumarize the samples for SAR type III
extract_species_positions <- function(species_habitat_matrix,
                                      species_site_matrix) {
  # Get habitat names
  habitat_names <- colnames(species_habitat_matrix[,-1])

  # Initialize a list to store species positions for each habitat
  habitat_positions <- list()

  # Loop through each habitat
  for (habitat in habitat_names) {
    # Get species associated with this habitat
    species_in_habitat <-
      rownames(species_habitat_matrix)[species_habitat_matrix[, habitat] == 1]

    # Find positions (indices) of these species in the species-site matrix
    species_positions <- which(rownames(species_site_matrix) %in%
                                 species_in_habitat)

    # Store results in the list
    habitat_positions[[habitat]] <- species_positions
  }

  return(habitat_positions)
}


#####################################################
######MAIN FUNCTION FOR BUILDING THE SAR####################
################################################
#' Build a SAR...
#'
#' @description Build a SAR .....
#' @usage sar_build(data, sp, lurast, crs = 4326,
#' habNam, habVals, sq_width = 2000, cluster_size_vector)
#' @param data A presence-absence matrix in the form of a
#'   dataframe, where rows are sites. The first three columns
#'   are: site names, then the longitude and latitude of the
#'   sites (note these can be any x and y co-ordinates depending
#'   on the \code{crs} argument provided??). The remaining
#'   columns relate to species in the dataset and take the form
#'   of a standard presence-absence matrix (i.e., 1 for presence
#'   in a given site, 0 for absence).
#' @param sp Dataframe of species habitat / specialisation
#'   classifications (e.g., forest specialist, generalist). The
#'   first column must be the species names (which must match the
#'   species column names in \code{data}). The remaining columns
#'   are then the different habitat habitat / specialisation
#'   classifications. Note that these classifications are binary
#'   (i.e., 1 or 0) and a species must only be placed in one
#'   classification, and must have a classification (i.e., the
#'   rows sums of the data frame, minus the first column, must
#'   all equal 1).
#' @param lurast A land use raster of class \code{SpatRaster}.
#'   See \code{\link[terra]{SpatRaster}} for more details of this
#'   class.
#' @param crs Coordinate reference system, something suitable as
#'   input to \code{\link[sf]{st_crs}}. Default is 'NULL',
#'   whereby the coordinate reference system from \code{lurast}
#'   is used. It is up to the user the ensure the correct
#'   coordinate reference system is employed.
#' @param habNam 
#' @param habVals 
#' @param sq_width ...
#' @param cluster_size_vector ...
#' @details ***How it deals with NAs in land-use raster, considers them
#' #to be a habitat type (and all NA = same habitat type)
#' @return A dataframe where rows are sites and the first N columns
#'   are the areas of the different habitats in \code{lurast},
#'   and then the next n columns are the species richness values
#'   associated with each species habitat / specialisation
#'   classification in \code{sp}. These area and richness columns
#'   can be fed directly into {\link{sar_countryside}}. The last
#'   two columns provide the total (i.e., summed) area and
#'   species richness for each site. Note, that these two final
#'   columns need dropping before using {\link{sar_countryside}}.
#' @importFrom sf st_convex_hull st_as_sf
#' @importFrom purrr list_flatten
#' @examples
#' ##To run the example, the sarsData package is required,
#' ##which is hosted on github due to the data object file
#' ##sizes.
#' #install.packages("devtools")
#' #devtools::install_github("txm676/sarsData")
#' #library(sarsData)
#' #data(IberianBirds)
#' #PA <- IberianBirds[[1]]
#' #classif <- IberianBirds[[2]]
#' ##Load in the land-use raster, stored in extadata drive 
#' ##of sarsData
#' #syf <- system.file("extdata", "coa95_raster_20.tif",
#' #package = "sarsData")
#' #lu1995 <- terra::rast(syf)
#' 
#' ##Run the sar_build function
#' #res <- sar_build(data = PA, sp = classif, 
#' #lurast = lu1995, crs = NULL, 
#' #habNam =  c("Forest", "Agriculture", "Shrubland","Other"),
#' #habVals = c(1,2,3,NA), sq_width = 2000,
#' #cluster_size_vector = c(1,4,16,64,256))
#' 
#' ##Log-log power model using summed area and 
#' ##richness values across habitats
#' #lin_pow(data = data.frame(res$Area_Total, 
#' #res$Sp_Total))
#' 
#' ##Fit a countryside SAR model to the output of
#' ##sar_build()
#' #datacsar <- subset(res, select = -c(Area_Total, Sp_Total))
#' #s3 <- sar_countryside(data = datacsar, modType = "power",
#' #gridStart = "none", habNam = habitat_names, 
#' #spNam = species_group_names)
#' @export

sar_build <- function(data, sp, lurast,
                      crs = NULL,
                      sq_width = 2000,
                      cluster_size_vector){

  if (!is.data.frame(data) & !is.matrix(data)){
    stop("data should be a dataframe or matrix")
  }
  data <- as.data.frame(data)

  if ((ncol(data) - 3) != nrow(sp)){
    stop("number of species in data and sp does not match")
  }

  colnames(sp)[1] <- "species"
  if (!all(colnames(data[,-c(1:3)]) %in% sp$species)){
    stop("species names in data and sp do not match")
  }

  if (any(rowSums(sp[,-1]) != 1)){
    stop("Row sums in sp are not all equal to 1")
  }

  #Convert abundance to presence-absence
  if (any(data[,-c(1:3)] > 1)){
    data[,-c(1:3)] <- replace(data[,-c(1:3)], 
                              data[,-c(1:3)]>1, 1)
    warning("Some values in data not 1 - all such values converted to 1")
  }

  if (is.null(crs)){
    crs2 <- crs(lurast)
  } else {
    crs2 <- crs
  }
  
  points <- st_as_sf(data, coords=c("long","lat"),
                     crs=crs2)

  convex_hull <- st_convex_hull(st_union(points))

  squares_sf <- create_squares(points, width = sq_width)

  pt_in_clusters <- filter_points_in_clusters(points,
                                              squares_sf,
                                              cluster_size_vector)

  species_groups <- extract_species_positions(classif, points[,-1])

  samples <- lapply(pt_in_clusters, function(x) x$points)
  polygons <- lapply(pt_in_clusters, function(x) x$chulls)

  res <- summarize_samples(samples = list_flatten(samples), 
                           polygons = list_flatten(polygons),
                           habitat_raster = lurast, 
                           species_groups = species_groups,
                           species_group_names = names(species_groups))

  return(res)
}
