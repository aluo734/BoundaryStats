# number of subgraphs ----
#' @name n_subgraph
#' @title Number of subgraphs
#' @description Statistical test the for number of subgraphs, or sets of contiguous boundary elements, in the data.
#'
#' @param x A RasterLayer object with boundary elements.
#' @param null_distrib A list of probability functions output from boundary_null_distrib().
#'
#' @return The number of subgraphs in the RasterLayer and a p-value.
#' @examples
#' song_raster <- raster('song_dialect_boundaries.asc')
#' dialect_boundaries <- define_boundary(song_raster, 0.1)
#' null_dialect_boundary <- boundary_null_distrib(song_raster, threshold = 0.2,
#'                          n_iterations = 1000, model = 'fractional_brownian',
#'                          projection = 4326, fract_dim = 1.5)
#'
#' n_subgraphs(dialect_boundaries, dialect_boundary_null)
#'
#' @author Amy Luo
#' @references Jacquez, G.M., Maruca,I S. & Fortin M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' @export
n_subgraph <- function(x, null_distrib) {
  # clump continuous cells (including cells touching at corners) together; each clump has a different value, 1:n_clumps
  x_subgraphs <- raster::clump(x)

  # number of subgraphs
  n_subgraph = max(x_subgraphs@data@values, na.rm = T)

  # statistical test
  p <- null_distrib$n_subgraph(n_subgraph) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)

  names(n_subgraph) <- 'n subgraphs'
  names(p) <-'p-value'
  return(c(n_subgraph, p))
}

#length of the longest subgraph ----
#' @name max_subgraph
#' @title Length of the longest subgraph
#' @description Statistical test for the length of the longest subgraph, or set of contiguous boundary elements.
#'
#' @param x A RasterLayer object with boundary elements.
#' @param null_distrib A list of probability functions output from boundary_null_distrib().
#' @param projection Numeric. EPSG code of input raster layer.
#'
#' @return The length of the longest subgraph and a p-value.
#' @examples
#' song_raster <- raster('song_dialect_boundaries.asc')
#' dialect_boundaries <- define_boundary(song_raster, 0.1)
#' null_dialect_boundary <- boundary_null_distrib(song_raster, threshold = 0.2,
#'                          n_iterations = 1000, model = 'fractional_brownian',
#'                          projection = 4326, fract_dim = 1.5)
#'
#' max_subgraph(dialect_boundaries, dialect_boundary_null, projection = 4326)
#'
#' @author Amy Luo
#' @references Jacquez, G.M., Maruca,I S. & Fortin M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' @export
max_subgraph <- function(x, null_distrib, projection) {
  # binarize input data and set non-boundary cells to NA
  x_mat <- raster::as.matrix(x)
  for (row in 1:nrow(x_mat)) {
    for (col in 1:ncol(x_mat)) {
      if(!is.na(x_mat[row, col])) {
        if (x_mat[row, col] == 0) {x_mat[row, col] = NA}
      }
    }
  }

  # matrix -> raster -> polygons (cell = boundary element, polygon with contiguous boundary elements = subgraph)
  x_boundary <- raster::raster(x_mat)
  raster::extent(x_boundary) <- raster::extent(x)
  terra::crs(x_boundary) <- paste0('+init=epsg:', projection)
  xpolygon <- raster::rasterToPolygons(x_boundary, na.rm = T, dissolve = T) %>%
    sp::disaggregate(.) %>%
    terra::buffer(., width = 0.001, dissolve = T) %>%
    sp::disaggregate(.) %>%
    sf::st_as_sf(.)

  # calculate the lengths of the subgraphs and keep the longest length
  lengths <- c()
  for (i in 1:nrow(xpolygon)) {
    lengths <- sf::st_geometry(xpolygon[i,]) %>%
      sf::st_cast(., "POINT") %>%
      sf::st_distance(.) %>%
      max(.) %>%
      append(lengths, .) %>%
      as.numeric(.)
  }
  max_length <- max(lengths)

  # statistical test
  p <- null_distrib$longest_subgraph(max_length) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)

  names(p) <-'p-value'
  names(max_length) <- 'length of longest boundary'
  return(c(max_length, p))
}
