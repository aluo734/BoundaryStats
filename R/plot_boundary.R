#' @name plot_boundary
#' @title Map the boundary elements of two raster layers
#' @description
#' This is a wrapper function for ggplot2 that will produce a map of boundary
#' elements for two traits and show where boundary elements intersect.
#'
#' @param x A SpatRaster object with boundary elements.
#' @param y A SpatRaster object with boundary elements.
#' @param color Optional. A character vector of up to three colors (x boundary, y boundary, and overlapping elements).
#' @param trait_names Optional. A character vector with up to two elements (legend name for x and legend name for y).
#' @param output_raster Returns a SpatRaster object with the boundary elements of each trait and overlapping
#' boundary elements together in a single layer.
#'
#' @return A ggplot2 object.
#' 
#' @examples
#' data(T.cristatus)
#' T.cristatus <- terra::rast(T.cristatus_matrix, crs = T.cristatus_crs)
#' terra::ext(T.cristatus) <- T.cristatus_ext
#'
#' data(grassland)
#' grassland <- terra::rast(grassland_matrix, crs = grassland_crs)
#' terra::ext(grassland) <- grassland_ext
#'
#' Tcrist_boundaries <- define_boundary(T.cristatus, cat = TRUE)
#' grassland_boundaries <- define_boundary(grassland, threshold = 0.1)
#'
#' plot_boundary(Tcrist_boundaries, grassland_boundaries)
#'
#' @author Amy Luo
#' @export
plot_boundary <- function(x, y, color = NA, trait_names = NA, output_raster = FALSE) {
  # prep boundary layers to plot
  x_layer <- terra::as.data.frame(x, xy = TRUE, na.rm = TRUE) %>%
    .[.[,3] != 0,]
  colnames(x_layer) <- c('long', 'lat', 'values')
  y_layer <- terra::as.data.frame(y, xy = TRUE, na.rm = TRUE) %>%
    .[.[,3] != 0,]
  colnames(y_layer) <- c('long', 'lat', 'values')

  # make overlap layer
  y2 <- y*2
  overlap <- x+y2

  overlap_layer <- terra::as.data.frame(overlap, xy = TRUE, na.rm = TRUE)
  colnames(overlap_layer) <- c('long', 'lat', 'values')
  overlap_layer <- dplyr::mutate(overlap_layer, values = as.character(values),
                                 values = dplyr::case_when(values == '1' ~ 'Trait 1',
                                                           values == '2' ~ 'Trait 2',
                                                           values == '3' ~ 'Overlap')) %>%
    dplyr::filter(!is.na(values))

  # if there are inputs for colors and layer names, change the colors from default
  fill_col <- c('Trait 1' = '#6EC6CA', 'Trait 2' = '#CCABD8', 'Overlap' = '#055B5C')
  if (all(is.na(color))) {
    fill_col <- fill_col
  } else if (length(color) > 1) {
    if(!is.na(color[1])) {fill_col[1] = color[1]}
    if(!is.na(color[2])) {fill_col[2] = color[2]}
    if(!is.na(color[3])) {fill_col[3] = color[3]}
  }

  if(!is.na(trait_names[1])) {
    names(fill_col)[1] <- trait_names[1]
    overlap_layer <- dplyr::mutate(overlap_layer, values = ifelse(values == 'Trait 1', trait_names[1], values))
  }
  
  if(!is.na(trait_names[2])) {
    names(fill_col)[2] <- trait_names[2]
    overlap_layer <- dplyr::mutate(overlap_layer, values = ifelse(values == 'Trait 2', trait_names[2], values))
  }
  
  # make plot
  world_map <- ggplot2::map_data('world')
  p <- ggplot2::ggplot(mapping = ggplot2::aes(long, lat)) +
    ggplot2::geom_map(data = world_map, map = world_map, ggplot2::aes(map_id = region), fill = '#e8e8e8') +
    ggplot2::geom_tile(data = overlap_layer, map = ggplot2::aes(fill = values)) +
    ggplot2::scale_fill_manual(values = fill_col, breaks = names(fill_col)) +
    ggplot2::coord_sf(xlim = terra::ext(x)[1:2], ylim = terra::ext(x)[3:4]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::labs(x = 'Longitude', y = 'Latitude', fill = 'Boundary Type')

  print(p)
  
  # export raster
  if (output_raster == TRUE) {
    overlap <- terra::as.factor(overlap)
    levels(overlap) <- data.frame(value = c(0, 1, 2, 3), boundary = c('No Boundary', trait_names[1], trait_names[2], 'Overlap'))
    return(overlap)
  }
}
