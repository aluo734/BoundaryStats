#' @name define_boundary
#' @title Define the boundary elements of a SpatRaster with numeric data or boundary intensities
#' @description
#'
#' Defines boundary elements in a SpatRaster object.
#'
#' For categorical traits, boundary elements are defined as cells that neighbor a patch with a different value.
#'
#' If the trait is quantitative, function will keep a proportion of the cells with the highest boundary intensity values,
#' with a threshold chosen by the user. If the SpatRaster contains trait values, boundary intensity values can be
#' calculated (calculate_intensity = T) using a Sobel-Feldman operator.
#'
#' In some cases, there may be many redundant values for boundary intensity. If the threshold is set to cut off the values
#' at those points, the actual proportion of cells categorized as boundary elements would differ from the intended
#' threshold. The function reports the proportion of cells that are boundary elements, so that users can choose a
#' different value, if necessary.
#'
#' @param x A SpatRaster object.
#' @param cat TRUE if the input SpatRaster contains a categorical variable. default = FALSE.
#' @param threshold A value between 0 and 1. The proportion of cells to keep as boundary elements. default = 0.2.
#' @param calculate_intensity logical. If TRUE, calculate boundary intensity at each cell from trait data. default = FALSE.
#'
#' @return A SpatRaster object with cell values 1 for boundary elements and 0 for other cells
#' 
#' @examples \donttest{
#' data(grassland)
#' grassland <- terra::rast(grassland_matrix, crs = grassland_crs)
#' terra::ext(grassland) <- grassland_ext
#' 
#' grassland_boundaries <- define_boundary(grassland, threshold = 0.1)
#' }
#' 
#' @author Amy Luo
#' @references
#' Fortin, M.J. et al. (2000) Issues related to the detection of boundaries. Landscape Ecology, 15, 453-466.
#' Jacquez, G.M., Maruca,I S. & Fortin M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' @export
define_boundary <- function (x, cat = FALSE, threshold = 0.2, calculate_intensity = FALSE) {
  # if the use enters cat and calculate_intensity as TRUE, ignore calculate_intensity with a message
  if(cat == TRUE & calculate_intensity == TRUE) {message('Variable is categorical. Skipping boundary intensity calculation.')}
  
  # for quantitative traits
  if (cat == FALSE) {
    # if the raster contains trait values, estimate first partial derivatives of the cells in lon and lat directions
    if (calculate_intensity == TRUE) {x <- sobel_operator(x)}
    
    # sort the cell values from highest to lowest, then find the value above which only the threshold
    # proportion of cells would be kept
    threshold_value <- terra::values(x) %>%
      na.omit %>%
      sort(decreasing = TRUE) %>%
      head(., round(length(.) * threshold)) %>%
      min
    
    # some of the cell values are probably redundant, especially if they're integers, so ther exact proportion of
    # boundary elements is not the same as the given threshold. This message informs user of the actual proportion
    # of boundary elements, so they can change it if there are too many redundant values near that threshold
    proportion <- length(terra::cells(x)[na.omit(terra::values(x)) > threshold_value])/length(terra::cells(x))
    message(paste('Proportion of boundary elements:', round(proportion, 5)))

    # make raster with boundary elements (filter out cells with values below threshold)
    boundaries <- terra::rast(nrows = nrow(x), ncols = ncol(x), extent = terra::ext(x), crs = terra::crs(x))
    for (i in terra::cells(x)) {
      if (terra::values(x)[i] > threshold_value) {
        boundaries[i] = 1
      } else {
        boundaries[i] = 0
      }
    }
    
  # for categorical trait values
  } else if (cat == TRUE) {
    poly <- terra::as.polygons(x)
    mask <- terra::aggregate(poly) %>%
      terra::as.lines() %>%
      terra::rasterize(terra::rast(x))
    fill <- terra::aggregate(poly) %>%
      terra::rasterize(terra::rast(x), field = 0)
    
    boundaries <- terra::as.lines(poly) %>%
      terra::rasterize(terra::rast(x)) %>%
      terra::merge(fill) %>%
      terra::mask(mask, maskvalues = 1, updatevalue = 0)
    
    terra::crs(boundaries) <- terra::crs(x)
  } else {
    stop('cat should be TRUE or FALSE.')
  }
  
  # remove singleton boundary elements
  boundary_patches <- terra::patches(boundaries, directions = 8, zeroAsNA = TRUE) %>%
    terra::cells()
  neighbor_cells <- terra::adjacent(boundaries, directions = 'queen', cells = boundary_patches)
  neighbor_values <- matrix(nrow = nrow(neighbor_cells), ncol = ncol(neighbor_cells))
  rownames(neighbor_values) <- rownames(neighbor_cells)
  for (i in 1:length(neighbor_cells)) {neighbor_values[i] <- terra::values(boundaries)[neighbor_cells[i]]}
  neighbor_values[is.na(neighbor_values)] <- 0
  singletons <- rowSums(neighbor_values) %>%
    .[. == 0] %>%
    names %>%
    as.numeric
  if(length(intersect(boundary_patches, rownames(neighbor_cells))) < length(boundary_patches)) {singletons = singletons + 1}
  terra::values(boundaries)[singletons] <- 0
  
  return(boundaries)
}
