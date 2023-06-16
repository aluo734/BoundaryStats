# null distribution for overlap statistics ----
#' @name overlap_null_distrib
#' @title Null distribution for boundary overlap statistics
#' @description
#' Given two RasterLayer objects with the same extent and resolution, this function will simulate new
#' RasterLayer objects from neutral models for a specified number of iterations. In each iteration,
#' boundary statistics will be computed for the simulated surface. The function outputs probability
#' functions for three boundary statistics (directly overlapping boundary elements, minimum distance
#' between boundary elements in x to y, and minimum distance between elements in x and y).
#'
#' @param x A RasterLayer object. If rand_both = F, only this raster will be modeled.
#' @param y A RasterLayer object. If rand_both = F, this raster does not change.
#' @param rand_both TRUE if distribution of traits in x and y should be modeled.
#' @param x_convert TRUE if x contains numeric trait data that needs to be converted to boundary intensities. default = FALSE.
#' @param y_convert TRUE if y contains numeric trait data that needs to be converted to boundary intensities. default = FALSE.
#' @param x_cat TRUE if x contains a categorical variable. default = FALSE.
#' @param y_cat TRUE if y contains a categorical variable. default = FALSE.
#' @param threshold A value between 0 and 1. The proportion of cells to keep as
#' boundary elements. default = 0.2.
#' @param n_iterations An integer indicating the number of iterations for the function. A value of 100 or 1000
#' is recommended to produce sufficient resolution for downstream statistical tests. default = 10.
#' @param x_model Neutral model to use. Options: 'random', 'gaussian_random' (Gaussian random field),
#' 'autoregressive' (simultaneous autoregressive model), and 'random_cluster' (modified random clusters method)
#' @param y_model Neutral model to use for y.
#' @param projection Numeric. EPSG code of input raster layer.
#' @param progress If progress = TRUE (default) a progress bar will be displayed.
#'
#' @return A list of probability distribution functions for boundary overlap statistics.
#' @examples
#' song_dialect_raster <- raster('song_dialect_boundaries.asc')
#' genetic_raster <- raster('genetic_group_assignments.asc')
#' null_overlap <- overlap_null_distrib(song_dialect_raster, genetic_raster, y_convert = T,
#'                                      threshold = 0.1, n_iterations = 1000,
#'                                      x_model = 'gaussian_random',
#'                                      y_model = 'gaussian_random',
#'                                      projection = 4326)
#'
#' @author Amy Luo
#' @references
#' Jacquez, G.M., Maruca, I.S. & Fortin, M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' James, P. M. A., Fleming, R.A., & Fortin, M.-J. (2010) Identifying significant scale-specific spatial boundaries using wavelets and null models: Spruce budworm defoliation in Ontario, Canada as a case study. Landscape Ecology, 6, 873-887.
#' @export
overlap_null_distrib <- function(x, y, rand_both, x_convert = F, y_convert = F, x_cat = F, y_cat = F, threshold = 0.2,
                                 n_iterations = 10, x_model = 'random', y_model = 'random', projection, progress = T) {
  # progress bar
  if (progress == T) {progress_bar = txtProgressBar(min = 0, max = n_iterations, initial = 0, char = '+', style = 3)}
  
  suppressWarnings({
    # make output vectors for repeat loop
    Odirect = c(); Ox = c(); Oxy = c()
    
    # if not randomizing x model, find boundaries in x
    if (rand_both == F) {
      if (y_cat == F) {
        y_boundary <- define_boundary(y, threshold, y_convert)
      } else {
        y_boundary <- categorical_boundary(y, projection)
      }
    }
    
    # n iterations of random boundaries + stats
    rep = 0
    repeat {
      
      # make random rasters with boundary elements (only for x if not randomizing y) ----
      repeat {
        x_sim <- simulate_rasters(x, x_model, projection)
        
        if (x_cat == F) {
          x_boundary <- define_boundary(x_sim, threshold, x_convert)
        } else {
          x_boundary <- categorical_boundary(x_sim, projection)
        }
        
        if (sum(x_boundary@data@values == 1, na.rm = T) >= 1) {break}
      }
      
      if (rand_both == T) {
        repeat {
          y_sim <- simulate_rasters(y, y_model, projection)
          
          if (y_cat == F) {
            y_boundary <- define_boundary(y_sim, threshold, y_convert)
          } else {
            y_boundary <- categorical_boundary(y_sim, projection)
          }
          
          if (sum(y_boundary@data@values == 1, na.rm = T) >= 1) {break}
        }
      }
      
      # make matrices to calculate stats ----
      x_mat <- raster::as.matrix(x_boundary)
      y_mat <- raster::as.matrix(y_boundary)
      
      # direct overlap statistic ----
      count = 0
      for (row in 1:nrow(x_mat)) {
        for (col in 1:ncol(x_mat)) {
          if (!is.na(x_mat[row, col]) & !is.na(y_mat[row, col])) {
            if (x_mat[row, col] == 1 & y_mat[row, col] == 1) {count = count + 1}
          }
        }
      }
      
      Odirect <- append(Odirect, count)
      
      # average minimum distance statistic (Ox; one boundary influences the other) ----
      min_distances <- c()
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
      
      Ox <- mean(min_distances) %>%
        append(Ox, .)
      
      # average minimum distance statistic (Oxy; both boundaries influence each other) ----
      min_distances <- c()
      x_present <- which(x_mat == 1, arr.ind = T)
      
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
      
      Oxy <- mean(min_distances) %>%
        append(Oxy, .)
      
      # loop count and break ----
      rep = rep + 1
      if (progress == T) {setTxtProgressBar(progress_bar, rep)}
      if (rep >= n_iterations) {
        if (progress == T) {close(progress_bar)}
        break
      }
    }
    
    # output
    Odirect <- pdqr::new_p(Odirect, type = 'continuous')
    Ox <- pdqr::new_p(Ox, type = 'continuous')
    Oxy <- pdqr::new_p(Oxy, type = 'continuous')
    distribs <- list(Odirect, Ox, Oxy)
    names(distribs) <- c('Odirect', 'Ox', 'Oxy')
    
    cat('DONE')
  })
  
  return(distribs)

}

# null distribution for boundary statistics ----
#' @name boundary_null_distrib
#' @title Null distribution for overlap statistics
#' @description
#' Given a RasterLayer object, this function will simulate a new RasterLayer from a neutral model for a
#' specified number of iterations. In each iteration, boundary statistics will be computed for the simulated
#' surface. The function outputs probability functions for two boundary statistics (number of subgraphs
#' and length of the longest subgraph).
#'
#' @param x A RasterLayer object.
#' @param convert TRUE if x contains numeric trait data that needs to be converted to boundary intensities. default = FALSE.
#' @param cat TRUE if the input RasterLayer contains a categorical variable. default = FALSE.
#' @param threshold A value between 0 and 1. The proportion of cells to keep as boundary elements. default = 0.2.
#' @param n_iterations An integer indicating the number of iterations for the function.
#' A value of 100 or 1000 is recommended to produce sufficient resolution for downstream
#' statistical tests. default = 10.
#' @param model Neutral model to use. Options: 'random', 'gaussian_random' (Gaussian random field),
#' 'autoregressive' (simultaneous autoregressive model), and 'random_cluster' (modified random clusters method)
#' @param projection Numeric. EPSG code of input raster layer.
#' @param progress If progress = TRUE (default) a progress bar will be displayed.
#'
#' @return A list of two probability distribution functions for boundary statistics.
#' @examples
#' song_dialect_raster <- raster('song_dialect_boundaries.asc')
#' null_dialect_boundary <- boundary_null_distrib(song_raster, threshold = 0.2,
#'                          n_iterations = 1000, model = 'gaussian_random',
#'                          projection = 4326)
#'
#' soil_type_raster <- raster('soil_types.asc')
#' boundary_null_distrib(soil_type_raster, cat = T, n_iterations = 100,
#'                       model = 'random_cluster', projection = 4210)
#'
#' @author Amy Luo
#' @references
#' Jacquez, G.M., Maruca, I.S. & Fortin, M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' James, P. M. A., Fleming, R.A., & Fortin, M.-J. (2010) Identifying significant scale-specific spatial boundaries using wavelets and null models: Spruce budworm defoliation in Ontario, Canada as a case study. Landscape Ecology, 6, 873-887.
#' @export
boundary_null_distrib <- function(x, convert = F, cat = F, threshold = 0.2, n_iterations = 10, model = 'random',
                                  projection, progress = T) {
  # progress bar
  if (progress == T) {progress_bar = txtProgressBar(min = 0, max = n_iterations, initial = 0, char = '+', style = 3)}

  suppressWarnings({
    # make output vectors for repeat loop
    n_subgraph = c(); longest_subgraph = c()
    
    # n iterations of random boundaries + stats
    rep = 0
    repeat {
      # simulate random raster with boundary elements ----
      repeat {
        x_sim <- simulate_rasters(x, model, projection)
        
        if (cat == F) {
          x_boundary <- define_boundary(x_sim, threshold, convert)
        } else {
          x_boundary <- categorical_boundary(x_sim, projection)
        }
        
        if (sum(x_boundary@data@values == 1, na.rm = T) >= 1) {break}
      }
      
      # number of subgraphs ----
      x_subgraphs <- raster::clump(x_boundary)
      n_subgraph = max(x_subgraphs@data@values, na.rm = T)
      
      # maximum length of a subgraph ----
      for (row in 1:nrow(x_boundary)) {
        for (col in 1:ncol(x_boundary)) {
          if(!is.na(x_boundary[row, col])) {
            if (x_boundary[row, col] == 0) {x_boundary[row, col] = NA}
          }
        }
      }
      
      xpolygon <- raster::rasterToPolygons(x_boundary, na.rm = T) %>%
        .[.@data$layer == 1,]
      raster::crs(xpolygon) <- paste0('+init=EPSG:', 4210)
      xpolygon <- terra::buffer(xpolygon, width = 0.001, dissolve = T) %>%
        sf::st_as_sf(.)
      
      lengths <- c()
      for (i in 1:nrow(xpolygon)) {
        lengths <- sf::st_geometry(xpolygon[i,]) %>%
          sf::st_cast(., 'POINT') %>%
          sf::st_distance(.) %>%
          max(.) %>%
          append(lengths, .) %>%
          as.numeric(.)
      }
      
      longest_subgraph <- max(lengths) %>%
        append(longest_subgraph, .)
      
      # loop count and break ----
      rep = rep + 1
      if (progress == T) {setTxtProgressBar(progress_bar, rep)}
      if (rep >= n_iterations) {
        if (progress == T) {close(progress_bar)}
        break
      }
    }
    
    # output
    n_subgraph <- pdqr::new_p(n_subgraph, type = 'continuous')
    longest_subgraph <- pdqr::new_p(longest_subgraph, type = 'continuous')
    distribs <- list(n_subgraph, longest_subgraph)
    names(distribs) <- c('n_subgraph', 'longest_subgraph')
    
    cat('DONE')
  })
  
  return(distribs)

}
