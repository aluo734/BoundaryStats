# number of overlapping boundary elements ----
#' @name n_overlap_boundaries
#' @title Direct overlap between boundary elements.
#' @description Statistical test for the number of directly overlapping boundary elements of two traits.
#'
#' @param x A SpatRaster object with boundary elements.
#' @param y A SpatRaster object with boundary elements.
#' @param null_distrib A list of probability functions output from overlap_null_distrib().
#'
#' @return The number of directly overlapping boundary elements and a p-value.
#' @examples \donttest{
#' data(T.cristatus)
#' T.cristatus <- terra::rast(T.cristatus_matrix, crs = T.cristatus_crs)
#' terra::ext(T.cristatus) <- T.cristatus_ext
#' 
#' data(grassland)
#' grassland <- terra::rast(grassland_matrix, crs = grassland_crs)
#' terra::ext(grassland) <- grassland_ext
#' 
#' Tcrist_ovlp_null <- overlap_null_distrib(T.cristatus, grassland, rand_both = FALSE,
#'   x_cat = TRUE, n_iterations = 100, x_model = 'random_cluster')
#' Tcrist_boundaries <- define_boundary(T.cristatus, cat = TRUE)
#' grassland_boundaries <- define_boundary(grassland, cat = FALSE, threshold = 0.1)
#' 
#' n_overlap_boundaries(Tcrist_boundaries, grassland_boundaries, Tcrist_ovlp_null)
#' }
#'
#' @author Amy Luo
#' @references
#' Jacquez, G.M., Maruca,I S. & Fortin, M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' Fortin, M.-J., Drapeau, P. & Jacquez, G.M. (1996) Quantification of the Spatial Co-Occurrences of Ecological Boundaries. Oikos, 77, 51-60.
#' @export
n_overlap_boundaries <- function(x, y, null_distrib) {
  xx <- terra::cells(x, 1)[[1]]
  yy <- terra::cells(y, 1)[[1]]
  n_overlapping <- length(intersect(xx, yy))

  p <- null_distrib$n_overlapping(n_overlapping) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)
  names(p) <-'p-value'
  names(n_overlapping) <- 'n_overlapping'

  return(c(n_overlapping, p))
}

# average minimum distance from boundaries in x to boundaries in y ----
#' @name average_min_x_to_y
#' @title Average minimum distance from x boundary elements to nearest y boundary element.
#' @description
#' Statistical test for the average minimum distance between each boundary element in raster x
#' and the nearest boundary element in raster y. Uses Euclidean distance. The boundaries of
#' trait x depend on the boundaries of trait y.
#'
#' @param x A SpatRaster object with boundary elements.
#' @param y A SpatRaster object with boundary elements.
#' @param null_distrib A list of probability functions output from overlap_null_distrib().
#'
#' @return The average minimum distance and a p-value.
#' 
#' @examples \donttest{
#' data(T.cristatus)
#' T.cristatus <- terra::rast(T.cristatus_matrix, crs = T.cristatus_crs)
#' terra::ext(T.cristatus) <- T.cristatus_ext
#' 
#' data(grassland)
#' grassland <- terra::rast(grassland_matrix, crs = grassland_crs)
#' terra::ext(grassland) <- grassland_ext
#' 
#' Tcrist_ovlp_null <- overlap_null_distrib(T.cristatus, grassland, rand_both = FALSE,
#'   x_cat = TRUE, n_iterations = 100, x_model = 'random_cluster')
#' Tcrist_boundaries <- define_boundary(T.cristatus, cat = TRUE)
#' grassland_boundaries <- define_boundary(grassland, cat = FALSE, threshold = 0.1)
#' 
#' average_min_x_to_y(Tcrist_boundaries, grassland_boundaries, Tcrist_ovlp_null)
#' }
#'
#' @author Amy Luo
#' @references
#' Jacquez, G.M., Maruca,I S. & Fortin,M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' Fortin, M.-J., Drapeau, P. & Jacquez, G.M. (1996) Quantification of the Spatial Co-Occurrences of Ecological Boundaries. Oikos, 77, 51-60.
#' @export
average_min_x_to_y <- function(x, y, null_distrib) {
  x_min_distances <- c()

  x_bound_cells <- terra::xyFromCell(x, terra::cells(x, 1)[[1]])
  y_bound_cells <- terra::xyFromCell(y, terra::cells(y, 1)[[1]])
  dists <- terra::distance(x_bound_cells, y_bound_cells, lonlat = T)
  for (i in sequence(nrow(dists))) {x_min_distances <- append(x_min_distances, min(dists[i,]))} # for each x boundary cell, the minimum distance to a y boundary cell

  avg_min_x_to_y <- mean(x_min_distances)
  names(avg_min_x_to_y) <- 'avg_min_x_to_y'

  p <- null_distrib$avg_min_x_to_y(avg_min_x_to_y) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)
  names(p) <- 'p-value'
  return(c(avg_min_x_to_y, p))
}

# average minimum distance between boundaries (reciprocal) ----
#' @name average_min_distance
#' @title Average minimum distance between boundary elements of two variables
#' @description
#' Statistical test for the average minimum distance between boundary elements in two raster layers.
#' Uses Euclidean distance. Boundaries for each trait affect one another reciprocally (x affects y
#' and y affects x).
#'
#' @param x A SpatRaster object with boundary elements.
#' @param y A SpatRaster object with boundary elements.
#' @param null_distrib A list of probability functions output from overlap_null_distrib().
#'
#' @return p-value
#' @examples \donttest{
#' data(T.cristatus)
#' T.cristatus <- terra::rast(T.cristatus_matrix, crs = T.cristatus_crs)
#' terra::ext(T.cristatus) <- T.cristatus_ext
#' 
#' data(grassland)
#' grassland <- terra::rast(grassland_matrix, crs = grassland_crs)
#' terra::ext(grassland) <- grassland_ext
#' 
#' Tcrist_ovlp_null <- overlap_null_distrib(T.cristatus, grassland, rand_both = FALSE,
#'   x_cat = TRUE, n_iterations = 100, x_model = 'random_cluster')
#' Tcrist_boundaries <- define_boundary(T.cristatus, cat = TRUE)
#' grassland_boundaries <- define_boundary(grassland, cat = FALSE, threshold = 0.1)
#' 
#' average_min_distance(Tcrist_boundaries, grassland_boundaries, Tcrist_ovlp_null)
#' }
#'
#' @author Amy Luo
#' @references
#' Jacquez, G.M., Maruca,I S. & Fortin, M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' Fortin, M.-J., Drapeau, P. & Jacquez, G.M. (1996) Quantification of the Spatial Co-Occurrences of Ecological Boundaries. Oikos, 77, 51-60.
#' @export
average_min_distance <- function(x, y, null_distrib) {
  min_distances <- c()
  
  x_bound_cells <- terra::xyFromCell(x, terra::cells(x, 1)[[1]])
  y_bound_cells <- terra::xyFromCell(y, terra::cells(y, 1)[[1]])
  dists <- terra::distance(x_bound_cells, y_bound_cells, lonlat = T)
  for (i in sequence(nrow(dists))) {min_distances <- append(min_distances, min(dists[i,]))} # for each x boundary cell, the minimum distance to a y boundary cell
  for (i in sequence(ncol(dists))) {min_distances <- append(min_distances, min(dists[,i]))} # for each x boundary cell, the minimum distance to a y boundary cell
  
  avg_min_dist <- mean(min_distances)
  names(avg_min_dist) <- 'avg_min_dist'

  p <- null_distrib$avg_min_dist(avg_min_dist) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)
  names(p) <- 'p-value'
  return(c(avg_min_dist, p))
}
