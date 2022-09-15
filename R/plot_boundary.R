#' @name plot_boundary
#' @title Map the boundary elements of two raster layers
#' @description
#' This is a wrapper function for ggplot2 that will produce a map of boundary
#' elements for two traits and show where boundary elements intersect.
#'
#' @param x A RasterLayer object with boundary elements.
#' @param y A RasterLayer object with boundary elements.
#' @param color Optional. A character vector of up to three colors (x boundary, y boundary, and overlapping elements).
#' @param trait_names Optional. A character vector with up to two elements (legend name for x and legend name for y).
#'
#' @return A ggplot object.
#' @examples
#' plot_boundary(soil_boundary, genetic_boundary, col = c('pink', 'orange', 'red'),
#'               trait_names = c('Soil Type', 'Genetic Group'))
#'
#' @author Amy Luo
#' @export
plot_boundary <- function(x, y, color = NA, trait_names = NA) {
  # prep boundary layers to plot
  x_layer <- raster::as.data.frame(x, xy = T, na.rm = T) %>%
    .[.[,3] != 0,]
  colnames(x_layer) <- c('lon', 'lat', 'values')
  y_layer <- raster::as.data.frame(y, xy = T, na.rm = T) %>%
    .[.[,3] != 0,]
  colnames(y_layer) <- c('lon', 'lat', 'values')

  # make overlap layer
  x_mat <- as.matrix(x)
  y_mat <- as.matrix(y)
  overlap <- matrix(NA, nrow = nrow(x_mat), ncol = ncol(x_mat))

  for (row in 1:nrow(x_mat)) {
    for (col in 1:ncol(x_mat)) {
      if (!is.na(x_mat[row, col]) && !is.na(y_mat[row, col])) {
        if(x_mat[row, col] & y_mat[row, col]) {overlap[row, col] = 1}
        }
    }
  }

  overlap <- raster::raster(overlap)
  raster::extent(overlap) <- raster::extent(x)
  overlap <- raster::as.data.frame(overlap, xy = T, na.rm = T)
  colnames(overlap) <- c('lon', 'lat', 'values')

# if there are inputs for colors and layer names, change the colors from default
  fill_col <- c('Trait 1' = '#6EC6CA', 'Trait 2' = '#CCABD8', 'Overlap' = '#055B5C')
  if (all(is.na(color))) {
    fill_col <- fill_col
  } else if (length(color) > 1) {
    if(!is.na(color[1])) {fill_col[1] = color[1]}
    if(!is.na(color[2])) {fill_col[2] = color[2]}
    if(!is.na(color[3])) {fill_col[3] = color[3]}
  }

  if(!is.na(trait_names[1])) {names(fill_col)[1] <- trait_names[1]}
  if(!is.na(trait_names[2])) {names(fill_col)[2] <- trait_names[2]}

  # make plot
  range <- raster::extent(x)
  lat_range <- c(range@ymin, range@ymax)
  lon_range <- c(range@xmin, range@xmax)

  p <- ggplot2::ggplot(mapping = ggplot2::aes(lon, lat)) +
    ggplot2::geom_sf(data = rnaturalearth::ne_countries(returnclass = 'sf'), ggplot2::aes(x = NULL, y = NULL),
            fill = '#f0f0f0', color = NA) +
    ggplot2::geom_raster(data = x_layer, ggplot2::aes(fill = names(fill_col)[1])) +   # first boundary layer
    ggplot2::geom_raster(data = y_layer, ggplot2::aes(fill = names(fill_col)[2])) +   # second boundary layer
    ggplot2::geom_raster(data = overlap, ggplot2::aes(fill = names(fill_col)[3])) +   # overlapping boundary elements
    ggplot2::scale_fill_manual(values = fill_col) +
    ggplot2::coord_sf(xlim = lon_range, ylim = lat_range) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = 'Longitude', y = 'Latitude', fill = 'Boundary Type')
  return(p)
}
