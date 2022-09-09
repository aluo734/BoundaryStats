# Odirect ----
#' @name Odirect
#' @title Direct overlap between boundary elements
#' @description Statistical test for the number of directly overlapping boundary elements of two traits.
#'
#' @param x A RasterLayer object with boundary elements.
#' @param y A RasterLayer object with boundary elements.
#' @param null_distrib A list of functions representing probability distributions. Output from overlap_null_distrib().
#'
#' @return p-value
#' @examples
#' soil_raster <- raster('soil_types.asc')
#' genetic_raster <- raster('genetic_assignment_probabilities.asc')
#'
#' overlap_null_distribs <- overlap_null_distrib(soil_raster, genetic_raster, n_iterations = 1000,
#'                                               y_convert = T, threshold = 0.1, projection = 4326)
#' soil_boundaries <- define_boundary(soil_raster)
#' genetic_boundaries <- define_boundary(genetic_raster)
#'
#' Odirect(soil_boundary, genetic_boundary, overlap_null_distribs)
#'
#' @author Amy Luo
#' @references
#' Jacquez, G.M., Maruca,I S. & Fortin, M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' Fortin, M.-J., Drapeau, P. & Jacquez, G.M. (1996) Quantification of the Spatial Co-Occurrences of Ecological Boundaries. Oikos, 77, 51-60.
#' @export
Odirect <- function(x, y, null_distrib) {
  x_mat <- raster::as.matrix(x)
  y_mat <- raster::as.matrix(y)

  count = 0
  for (row in 1:nrow(x_mat)) {
    for (col in 1:ncol(x_mat)) {
      if (!is.na(x_mat[row, col]) & !is.na(y_mat[row, col])) {
        if (x_mat[row, col] == 1 & y_mat[row, col] == 1) {count = count + 1}
      }
    }
  }

  p <- null_distrib$Odirect(count) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)
  names(p) <-'p-value'
  names(count) <- 'n overlapping boundary elements'

  return(c(count, p))
}

# Ox ----
#' @name Ox
#' @title Average minimum distance from x boundary elements to nearest y boundary element
#' @description
#' Statistical test for the average minimum distance between each boundary element in raster x
#' and the nearest boundary element in raster y. Uses Euclidean distance. The boundaries of
#' trait x depend on the boundaries of trait y.
#'
#' @param x A RasterLayer object with boundary elements.
#' @param y A RasterLayer object with boundary elements.
#' @param null_distrib A list of functions representing probability distributions. Output from overlap_null_distrib().
#'
#' @return p-value
#' @examples
#' soil_raster <- raster('soil_types.asc')
#' genetic_raster <- raster('genetic_assignment_probabilities.asc')
#'
#' overlap_null_distribs <- overlap_null_distrib(soil_raster, genetic_raster, n_iterations = 1000,
#'                                               y_convert = T, threshold = 0.1, projection = 4326)
#' soil_boundaries <- define_boundary(soil_raster)
#' genetic_boundaries <- define_boundary(genetic_raster)
#'
#' Ox(soil_boundary, genetic_boundary, overlap_null_distribs)
#'
#' @author Amy Luo
#' @references
#' Jacquez, G.M., Maruca,I S. & Fortin,M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' Fortin, M.-J., Drapeau, P. & Jacquez, G.M. (1996) Quantification of the Spatial Co-Occurrences of Ecological Boundaries. Oikos, 77, 51-60.
#' @export
Ox <- function(x, y, null_distrib) {
  min_distances <- c()

  x_mat <- raster::as.matrix(x)
  y_mat <- raster::as.matrix(y)
  y_present <- which(y_mat == 1, arr.ind = T)

  for (row in 1:nrow(x_mat)) {
    for (col in 1:ncol(x_mat)) {
      # for each boundary element in x
      if (!is.na(x_mat[row, col])) {
        if (x_mat[row, col] == 1) {
          # calculate the Euclidean distances to each boundary element in y
          distances <- matrix(nrow = nrow(y_present))
          for (i in 1:length(distances)) {distances[i] <- sqrt((y_present[i, 1] - row)^2 + (y_present[i, 2] - col)^2)}
          # and keep the distance to the closest boundary element in y
          min_distances <- append(min_distances, min(distances))
        }
      }
    }
  }

  stat <- mean(min_distances) # average the minimum distances
  names(stat) <- 'average minimum distance (x depends on y)'

  p <- null_distrib$Ox(stat) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)
  names(p) <- 'p-value'
  return(c(stat, p))
}

# Oxy ----
#' @name Oxy
#' @title Average minimum distance between boundary elements of two variables
#' @description
#' Statistical test for the average minimum distance between boundary elements in two raster layers.
#' Uses Euclidean distance. Boundaries for each trait affect one another reciprocally (x affects y
#' and y affects x).
#'
#' @param x A RasterLayer object with boundary elements.
#' @param y A RasterLayer object with boundary elements.
#' @param null_distrib A list of functions representing probability distributions. Output from overlap_null_distrib().
#'
#' @return p-value
#' @examples
#' soil_raster <- raster('soil_types.asc')
#' genetic_raster <- raster('genetic_assignment_probabilities.asc')
#'
#' overlap_null_distribs <- overlap_null_distrib(soil_raster, genetic_raster, n_iterations = 1000,
#'                                               y_convert = T, threshold = 0.1, projection = 4326)
#' genetic_boundaries <- define_boundary(genetic_raster)
#'
#' Oxy(soil_boundary, genetic_boundary, overlap_null_distribs)
#'
#' @author Amy Luo
#' @references
#' Jacquez, G.M., Maruca,I S. & Fortin, M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' Fortin, M.-J., Drapeau, P. & Jacquez, G.M. (1996) Quantification of the Spatial Co-Occurrences of Ecological Boundaries. Oikos, 77, 51-60.
#' @export
Oxy <- function(x, y, null_distrib) {
  min_distances <- c()

  x_mat <- raster::as.matrix(x)
  y_mat <- raster::as.matrix(y)
  x_present <- which(x_mat == 1, arr.ind = T)
  y_present <- which(y_mat == 1, arr.ind = T)

  for (row in 1:nrow(x_mat)) {
    for (col in 1:ncol(x_mat)) {
      if (!is.na(x_mat[row, col])) {
        if (x_mat[row, col] == 1) {
          distances <- matrix(nrow = nrow(y_present))
          for (i in 1:length(distances)) {distances[i] <- sqrt((y_present[i, 1] - row)^2 + (y_present[i, 2] - col)^2)}
          min_distances <- append(min_distances, min(distances))
        }
      }
    }
  }

  for (row in 1:nrow(y_mat)) {
    for (col in 1:ncol(y_mat)) {
      if (!is.na(y_mat[row, col])) {
        if (y_mat[row, col] == 1) {
          distances <- matrix(nrow = nrow(x_present))
          for (i in 1:length(distances)) {distances[i] <- sqrt((x_present[i, 1] - row)^2 + (x_present[i, 2] - col)^2)}
          min_distances <- append(min_distances, min(distances))
        }
      }
    }
  }

  stat <- mean(min_distances)
  names(stat) <- 'average minimum distance'

  p <- null_distrib$Oxy(stat) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)
  names(p) <- 'p-value'
  return(c(stat, p))
}
