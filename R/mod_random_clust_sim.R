#' @name mod_random_clust_sim
#' @title Modified random cluster neutral landscape model
#' @description Simulates a neutral landscape of the same extent and resolution as the input raster, with the same
#' distribution of values.
#'
#' @param x A SpatRaster object.
#' @return A SpatRaster object with boundary elements.
#' 
#' @examples
#' data(grassland)
#' grassland <- terra::rast(grassland_matrix, crs = grassland_crs)
#' terra::ext(grassland) <- grassland_ext
#' 
#' simulation <- mod_random_clust_sim(grassland)
#' terra::plot(simulation)
#' 
#' @author Amy Luo
#' @references
#' Saura, S. & Martínez-Millán, J. (2000) Landscape patterns simulation with a modified random clusters method. Landscape Ecology, 15, 661 – 678.
#' @export
mod_random_clust_sim <- function (x) {
  # A: make percolated raster, where proportion of filled cells = p_cluster
  x_sim <- matrix(nrow = terra::nrow(x), ncol = terra::ncol(x))
  for (i in 1:length(x_sim)) {
    if (stats::runif(1) <= 0.5) {x_sim[i] = 1}
  }
  
  opp <- matrix(nrow = terra::nrow(x), ncol = terra::ncol(x)) # matrix with opposite cells marked for clumping
  opp[is.na(x_sim[])] <- 1
  
  # B: find clusters in raster
  x_sim <- terra::rast(x_sim, crs = terra::crs(x), ext = terra::ext(x)) %>%
    terra::patches(.)
  
  opp <- terra::rast(opp, crs = terra::crs(x), ext = terra::ext(x)) %>%
    terra::patches(.)
  terra::values(opp) <- terra::values(opp) + 0.5
  
  x_sim <- terra::merge(x_sim, opp)
  terra::values(x_sim)  <- terra::values(x) %>%
    na.omit(.) %>%
    max(.) + terra::values(x_sim)
  
  # C: assign clusters to categories
  prop <- terra::freq(x, digits = 2) %>% # proportions of cells in category
    as.data.frame(.) %>%
    tibble::add_column(., p = .$count/sum(.$count)) %>%
    .[order(-.$p),-1]
  
  clump_selection_order <- sample(terra::unique(x_sim))[,1]
  areas <- terra::cellSize(x_sim)
  total_area <- terra::expanse(x_sim)
  
  row = 1
  for (i in clump_selection_order) {
    assignment = prop[row, 1]; max_prop = prop[row, 3]
    terra::values(x_sim)[terra::values(x_sim) == i] <- assignment
    
    current_prop <- which(terra::values(x_sim) == assignment) %>%
      areas[.] %>%
      sum(.)/total_area %>%
      .[,2]
    if (current_prop >= max_prop) {row = row + 1}
    
    if (row > nrow(prop)) {break}
  }
  
  # crop extent of filled values to input data range
  x_sim <- terra::mask(x_sim, x)
  
  return(x_sim)
  
}