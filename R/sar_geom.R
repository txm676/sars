#' 
#' ##filter_points_in_expanding_circles
#' #A function that picks a random point and starts sampling increasing circles
#' #from those points, returning the points inside those circles.
#' 
#' #' @importFrom sf st_geometry st_buffer st_intersection st_area st_intersects
#' 
#' filter_points_in_expanding_circles <- function(points_sf, 
#'                                                radius_vector, 
#'                                                convex_hull) {
#'   # Randomly select one point from points_sf
#'   selected_point <- points_sf[sample(1:nrow(points_sf), 1), ]
#'   
#'   # Initialize a list to store points within each circle
#'   points_within_circles <- list()
#'   
#'   # Loop through the radius_vector
#'   for (radius in radius_vector) {
#'     # Create a buffer (circle) around the selected point with the current radius
#'     circle <- st_geometry(st_buffer(selected_point, dist = radius))
#'     
#'     # Calculate the intersection of the circle with the convex_hull
#'     intersection <- st_intersection(circle, convex_hull) 
#'     
#'     # Check the area of the intersection compared to the area of the circle
#'     circle_area <- as.numeric(st_area(circle)) #as.numeric to remove units
#'     
#'     intersection_area <- as.numeric(st_area(intersection)) #as.numeric to remove units
#'     
#'     # Stop if the majority of the circle is outside the convex_hull
#'     if (intersection_area / circle_area < 0.5) {
#'       break
#'     }
#'     
#'     # Find points_sf that intersect the current circle
#'     points_in_circle <- points_sf[st_intersects(points_sf, 
#'                                                 circle, 
#'                                                 sparse = FALSE), ]
#'     
#'     # Add the points to the list
#'     points_within_circles[[paste0("radius_", radius)]] <-
#'       list(points=points_in_circle,circle=circle)
#'   }
#'   
#'   # Return the list of points within circles
#'   return(points_within_circles)
#' }
#' 
#' 
#' ##summarize_samples
#' #A function that takes a list of groups of sites as a sample feature collections
#' #which has a dataframe with the species occurences, the polygons corresponding
#' #to the geometric delimitation of the groups of sites, a raster of habitat
#' #types, a vector of habitat names, the list of species in each group (each group
#' #is a list of numerical values of the columns corresponding to the species in
#' #the dataframe from the sample fature), and a vector with the species group
#' #names. The first column in the sample feature collection should be the name of
#' #the individual sample site. It produces a table with the area of each habitat
#' #and the number of species in each group and the total number of species.
#' 
#' #' @importFrom terra crop vect mask freq merge
#' #' @importFrom sf st_drop_geometry
#' 
#' summarize_samples <- function(samples, polygons, habitat_raster, habitat_names,
#'                              habitat_values, species_groups, species_group_names)
#' {
#'   # Initialize an empty data frame for the results
#'   results_df <- data.frame(matrix(ncol = length(habitat_names)+
#'                                     length(species_group_names)+2, 
#'                                   nrow = 0))
#'   colnames(results_df) <- c(habitat_names,"Area_Total",
#'                             species_group_names,"Sp_Total")
#'   
#'   # Iterate over each area sample (i.e. group of sites)
#'   for (i in seq_along(samples)) {
#'     sample <- samples[[i]]
#'     polygon <- polygons[[i]]
#'     
#'     # Crop and mask the raster to the polygon extent
#'     habitat_cropped <- crop(habitat_raster, vect(polygon))
#'     habitat_masked <- mask(habitat_cropped, vect(polygon))
#'     
#'     # Calculate the area of each habitat type
#'     habitat_df <- freq(habitat_masked, bylayer=FALSE)
#'     habitat_df <- rbind(habitat_df,freq(habitat_masked, 
#'                                         bylayer=FALSE,value=NA))
#'     
#'     # Ensure all possible values are included
#'     all_values_df <- data.frame(value = habitat_values, count = 0)
#'     
#'     # Merge with actual frequency data, replacing 0 where missing
#'     habitat_df <- merge(all_values_df, habitat_df,
#'                         by = "value", all.x = TRUE)
#'     
#'     # Fill NA counts with 0
#'     habitat_df$count <- ifelse(is.na(habitat_df$count.y), 
#'                                0, habitat_df$count.y)
#'     
#'     # Drop unnecessary column
#'     habitat_df <- habitat_df[, c("value", "count")]
#'     
#'     habitat_df$area <- habitat_df$count * 
#'       res(habitat_raster)[1] * res(habitat_raster)[2]
#'     
#'     # Store the results
#'     results_df[i, seq_along(habitat_names)] <- habitat_df$area # store area of each habitat
#'     results_df[i, length(habitat_names)+1] <- 
#'       sum(results_df[i, seq_along(habitat_names)]) # store total area
#'     
#'     # Subset species occurrences for each group
#'     total_species <- 0
#'     for (k in seq_along(species_groups))
#'     {
#'       group_species <- st_drop_geometry(sample)[, species_groups[[k]]+1 ]
#'       species_present <- colSums(group_species > 0)  # Count species occurrences
#'       num_species <- sum(species_present > 0)       # Number of species in the group
#'       results_df[i, length(habitat_names)+1+k] <- num_species #store species number in results
#'       total_species <- total_species + num_species
#'     }
#'     
#'     #store the total number of species by summing across species groups
#'     results_df[i, length(habitat_names) + 
#'                  length(species_groups)+2] <- total_species
#'   }
#'   
#'   return(results_df)
#' }
#' 
#' 
#' ##create_squares
#' #Create an sf where each sampling point is associated with a sampling square
#' #centred in the point and with a width dx.
#' 
#' #' @importFrom sf st_as_sf st_buffer
#' 
#' create_squares <- function(points_sf, width) {
#'   # Ensure the input is a point sf object
#'   if (!inherits(points_sf, "sf") || !inherits(st_geometry(points_sf), 
#'                                               "sfc_POINT")) {
#'     stop("Input must be an sf object with point geometries.")
#'   }
#'   # Calculate half-width (to shift the square corners)
#'   half_width <- width / 2
#'   
#'   squares_sf <- st_as_sf(st_buffer(points_sf, 
#'                                    dist = half_width, 
#'                                    endCapStyle = "SQUARE"))
#'   
#'   return(squares_sf)
#' }
#' 
#' 
#' ##filter_points_in_clusters
#' #Filter points for type II/III SAR (fractals sampling scheme or sub-divisions),
#' #Points should have an associated sampling square as geometry feature instead of
#' #a point geometry.
#' #' @importFrom terra split
#' #' @importFrom sf st_geometry st_coordinates
#' #' @importFrom stats kmeans
#' 
#' filter_points_in_clusters <- function(points_sf, squares_sf, 
#'                                       cluster_size_vector) 
#' {
#'   npoints <- nrow(points_sf)
#'   n_clusters_vector = npoints%/%cluster_size_vector
#'   n_clusters_vector[n_clusters_vector==0] <- 1  # when whole landscape, npoints<cluster_size 
#'   points_within_clusters<- list()
#'   for (n_clusters in n_clusters_vector) 
#'   {
#'     if(n_clusters==npoints)
#'     {
#'       points_in_clusters <- split(points_sf, 1:npoints)
#'       clusters_convex_hulls <- split(st_geometry(squares_sf), 1:npoints)
#'     }
#'     else
#'     {
#'       # Extract coordinates for k-means
#'       coords <- st_coordinates(points_sf)
#'       
#'       # Perform k-means clustering
#'       kmeans_result <- kmeans(coords, centers = n_clusters)
#'       
#'       # split the points sf into a list of clusters of points
#'       points_in_clusters <- split(points_sf, kmeans_result$cluster)
#'       squares_in_clusters <- split(squares_sf, kmeans_result$cluster)
#'       
#'       # Merge all squares into a single geometry
#'       merged_geom_within_clusters <- lapply(squares_in_clusters, st_union)
#'       
#'       # Compute convex hull
#'       clusters_convex_hulls <- lapply(merged_geom_within_clusters, 
#'                                       st_convex_hull)
#'       
#'       #clusters_convex_hulls <- lapply(points_in_clusters, 
#'       #                                   function(points) {st_convex_hull(st_union(points))})
#'     }
#'     # Add the points to the list
#'     points_within_clusters[[paste0("size_", 
#'             cluster_size_vector[n_clusters_vector==n_clusters])]] <-
#'       list(points=points_in_clusters,chulls=clusters_convex_hulls)
#'   }
#'   # Return the list of points within clusters
#'   return(points_within_clusters)
#' }
#' 
#' 
#' 
#' ##sumarize the samples for SAR type III
#' extract_species_positions <- function(species_habitat_matrix, 
#'                                       species_site_matrix) {
#'   # Get habitat names
#'   habitat_names <- colnames(species_habitat_matrix[,-1])
#'   
#'   # Initialize a list to store species positions for each habitat
#'   habitat_positions <- list()
#'   
#'   # Loop through each habitat
#'   for (habitat in habitat_names) {
#'     # Get species associated with this habitat
#'     species_in_habitat <- 
#'       rownames(species_habitat_matrix)[species_habitat_matrix[, habitat] == 1]
#'     
#'     # Find positions (indices) of these species in the species-site matrix
#'     species_positions <- which(rownames(species_site_matrix) %in% 
#'                                  species_in_habitat)
#'     
#'     # Store results in the list
#'     habitat_positions[[habitat]] <- species_positions
#'   }
#'   
#'   return(habitat_positions)
#' }
