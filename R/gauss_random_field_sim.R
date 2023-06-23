#' @name gauss_random_field_sim
#' @title Gaussian random field neutral model
#' @description Simulates a gaussian random field as a neutral landscape of the same extent and resolution as the
#' input raster, using the same spatial autocorrelation range as the input
#'
#' @param x A SpatRaster object.
#' @return A SpatRaster object with boundary elements.
#' 
#' @examples
#' data(grassland)
#' grassland <- terra::rast(grassland_matrix, crs = grassland_crs)
#' terra::ext(grassland) <- grassland_ext
#' 
#' simulation <- gauss_random_field_sim(grassland)
#' terra::plot(simulation)
#'
#' @author Amy Luo
#' @references
#' James, P. M. A., Fleming, R.A., & Fortin, M.-J. (2010) Identifying significant scale-specific spatial boundaries using wavelets and null models: Spruce budworm defoliation in Ontario, Canada as a case study. Landscape Ecology, 6, 873-887.
#' @export
gauss_random_field_sim <- function (x) {
  # estimate autocorrelation range using local Moran's I + LISA clustering
  cell_vals <- terra::cells(x) %>%
    terra::values(x)[.] %>%
    as.data.frame(.)

  lisa_clusters <- terra::as.polygons(x, dissolve = F) %>%
    sf::st_as_sf(.) %>%
    rgeoda::queen_weights(.) %>%
    rgeoda::local_moran(., cell_vals) %>%
    rgeoda::lisa_clusters(.) %>%
    data.frame(cellID = terra::cells(x), group = .)

  cells_to_fill <- terra::rowColFromCell(x, terra::cells(x))
  x_cluster <- terra::rast(nrow = terra::nrow(x), ncol = terra::ncol(x), crs = terra::crs(x), extent = terra::ext(x))
  index = 1
  for (i in sequence(nrow(lisa_clusters))) {
    x_cluster[cells_to_fill[i,1], cells_to_fill[i,2]] <- lisa_clusters[index, 2]
    index = index + 1
  }

  x_cluster <- terra::as.polygons(x_cluster, na.rm = TRUE) %>%
    terra::buffer(., 0.01) %>%
    terra::disagg(.) %>%
    sf::st_as_sf(.) %>%
    sf::st_area(.)

  cell_size <- terra::cellSize(x, transform = T) %>%
    terra::values(.) %>%
    mean(.)
  corr_range <- sqrt(as.numeric(median(x_cluster)))/sqrt(cell_size)

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

  terra::plot(x_sim)

  # if input values are all integers, make all simulated values integers
  if (all(na.omit(terra::values(x)) %% 1 == 0)) {terra::values(x_sim) <- round(terra::values(x_sim))}

  # output
  return(x_sim)

}
