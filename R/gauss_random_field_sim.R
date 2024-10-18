#' @name gauss_random_field_sim
#' @title Gaussian random field neutral model
#' @description Simulates a gaussian random field as a neutral landscape of the same extent and resolution as the
#' input raster, using the same spatial autocorrelation range as the input
#'
#' @param x A SpatRaster object.
#' @return A SpatRaster object with boundary elements.
#'
#' @examples \donttest{
#' #' data(grassland)
#' grassland <- terra::rast(grassland_matrix, crs = grassland_crs)
#' terra::ext(grassland) <- grassland_ext
#'
#' simulation <- gauss_random_field_sim(grassland)
#' terra::plot(simulation)
#' }
#'
#' @author Amy Luo
#' @references
#' James, P. M. A., Fleming, R.A., & Fortin, M.-J. (2010) Identifying significant scale-specific spatial boundaries using wavelets and null models: Spruce budworm defoliation in Ontario, Canada as a case study. Landscape Ecology, 6, 873-887.
#' @export
gauss_random_field_sim <- function (x) {
  # estimate autocorrelation range using local Moran's I + LISA clustering
  polys <- terra::as.polygons(x, dissolve = F) %>%
    sf::st_as_sf() %>%
    sf::st_buffer(0.01)
  names(polys)[1] <- 'values'
  
  x_cluster <- polys %>%
    dplyr::mutate(nb = sfdep::st_contiguity(geometry),
                  wt = sfdep::st_weights(nb, allow_zero = T),
                  local_moran = sfdep::local_moran(values, nb, wt)) %>%
    tidyr::unnest(local_moran) %>%
    dplyr::mutate(cluster = ifelse(p_ii <= 0.05, 1, NA)) %>%
    dplyr::select(cluster, geometry) %>%
    terra::rasterize(x, field = 'cluster') %>%
    terra::as.polygons(na.rm = T) %>%
    terra::buffer(., 0.01) %>%
    terra::disagg() %>%
    sf::st_as_sf()

  lengths <- c()
  for (i in 1:nrow(x_cluster)) {
    lengths <- sf::st_geometry(x_cluster[i,]) %>%
      sf::st_cast(., "POINT") %>%
      sf::st_distance() %>%
      max() %>%
      append(lengths, .) %>%
      as.numeric()
  }
  median_length <- median(lengths)

  cell_size <- terra::cellSize(x, transform = T) %>%
    terra::values() %>%
    mean
  corr_range <- median_length/sqrt(cell_size)

  # simulate raster
  repeat {
    invisible(capture.output(
      x_sim <- try(list(1:terra::nrow(x), 1:terra::ncol(x)) %>%
                      fields::circulantEmbeddingSetup(., cov.args = list(p = 2, aRange = corr_range)) %>%
                      fields::circulantEmbedding(.) %>%
                      terra::rast(.),
                    silent = TRUE)
               ))

    if(is(x_sim) == 'try-error') {corr_range = corr_range * 0.9} else {break}
  }

  # make extent and projection match input data
  terra::crs(x_sim) <- terra::crs(x)
  terra::ext(x_sim) <- terra::ext(x)
  x_sim <- terra::mask(x_sim, x)

  # transform value range of simulated raster
  x_sim <- terra::values(x) %>%
    na.omit(.) %>%
    range(.) %>%
    scales::rescale(terra::as.matrix(x_sim, wide = TRUE), to = .) %>%
    matrix(., nrow = nrow(x_sim), ncol = ncol(x_sim)) %>%
    terra::rast(.)

  terra::ext(x_sim) <- terra::ext(x)
  terra::crs(x_sim) <- terra::crs(x)

  # if input values are all integers, make all simulated values integers
  if (all(na.omit(terra::values(x)) %% 1 == 0)) {terra::values(x_sim) <- round(terra::values(x_sim))}

  # output
  return(x_sim)

}
