# define boundary with numeric data ----
#' @name define_boundary
#' @title Define the boundary elements of a RasterLayer with numeric data
#' @description
#'
#' Defines boundaries in raster objects by keeping a proportion of the cells with the highest
#' boundary intensity values. If cells have trait values, they can be converted to boundary
#' values based on approximations of the first partial derivatives of each cell, calculated
#' using a sobel operator.
#'
#' @importFrom raster values as.matrix extent crs
#'
#' @param x A RasterLayer object.
#' @param threshold A value between 0 and 1. The proportion of cells to keep as boundary elements. default = 0.2.
#' @param convert logical. If TRUE, will convert values of each cell from trait values to boundary intensities. default = FALSE.
#' value to the estimated rate of change at that cell using a sobel operator. default = FALSE.
#'
#' @return A RasterLayer object with cell values 1 for boundary element and NA for other cell
#' @examples
#' genetic_assignments <- raster('genetic_assignments.asc')
#' genetic_boundary <- define_boundary(genetic_assignments, threshold = 0.2, convert = T)
#'
#' genetic_boundary <- raster('genetic_boundary.asc')
#' genetic_boundary <- define_boundary(genetic_boundary, threshold = 0.1, convert = F)
#'
#' @author Amy Luo
#' @references
#' Fortin, M.J. et al. (2000) Issues related to the detection of boundaries. Landscape Ecology, 15, 453-466.
#' Jacquez, G.M., Maruca,I S. & Fortin M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' @export
define_boundary <- function (x, threshold = 0.2, convert = F) {
  # if the raster contains trait values, estimate first partial derivatives of the cells in lon and lat directions
  if (convert == T) {x <- spatialEco::sobal(x)}

  # sort the cell values from highest to lowest, then find the value above which only the threshold
  # proportion (default = 0.2) of cells would be kept
  threshold_value <- values(x) %>%
    na.omit(.) %>%
    sort(., decreasing = T) %>%
    utils::head(., round(length(.) * threshold)) %>%
    min(.)

  # Check proportion of cells above threshold value. Sometimes there are a lot of redundant cell values, so if there
  # are redundancies at the threshold value, then the proportion of cells kept will be higher than the threshold.
  prop <- values(x) %>%
    na.omit(.) %>%
    .[. >= threshold_value] %>%
    length(.)/length(na.omit(values(x)))

  # so we'll keep reducing the threshold value until the proportion of cells kept is below the threshold
  if (threshold_value %% 1 == 0) {                # sometimes the boundary values are integers
    while (prop > threshold) {
      threshold_value = threshold_value + 1       # so increase the threshold value by 1
      prop <- values(x) %>%
        na.omit(.) %>%
        .[. >= threshold_value] %>%
        length(.)/length(na.omit(values(x)))
      }
    } else {                                      # sometimes the boundary values are float
    rep = 0
    while (prop > threshold) {
      threshold_value = threshold_value * 1.1     # so increase the threshold value by 10%
      prop <- values(x) %>%                       # until the proportion of cells kept < threshold
        na.omit(.) %>%
        .[. >= threshold_value] %>%
        length(.)/length(na.omit(values(x)))
      rep = rep + 1
      if (rep >= 10) {break}                      # break after 10 reps if it gets that far
      }
    }

  # make raster with boundary elements (filter out cells with values below threshold)
  boundaries <- as.matrix(x)
  for (row in 1:nrow(boundaries)) {
    for (col in 1:ncol(boundaries))
      if(!is.na(boundaries[row, col])) {
        if (boundaries[row, col] < threshold_value) {boundaries[row, col] = 0} else {boundaries[row, col] = 1}
      }
  }

  # remove singleton boundary elements
  for (row in 1:nrow(boundaries)) {
    for (col in 1:ncol(boundaries)) {
      if (!is.na(boundaries[row, col])) {
        if(boundaries[row, col] == 1) {
          # name all neighboring cells
          if (row > 1) {up = boundaries[row - 1, col]} else {up = NA}
          if (row < nrow(boundaries)) {down = boundaries[row + 1, col]} else {down = NA}
          if (col > 1) {left = boundaries[row, col - 1]} else {left = NA}
          if (col < ncol(boundaries)) {right = boundaries[row, col + 1]} else {right = NA}

          if (row > 1 & col > 1) {upleft = boundaries[row - 1, col - 1]} else {upleft = NA}
          if (row > 1 & col < ncol(boundaries)) {upright = boundaries[row - 1, col + 1]} else {upright = NA}
          if (row < nrow(boundaries) & col > 1) {downleft = boundaries[row + 1, col - 1]} else {downleft = NA}
          if (row < nrow(boundaries) & col < ncol(boundaries)) {downright = boundaries[row + 1, col + 1]} else {downright = NA}

          # test whether any neighbors are boundary elements
          neighbors <- c(up, down, left, right, upleft, upright, downleft, downright)

          # if not, then set the focal cell to 0 (change focal cell to not boundary element)
          if(1 %in% neighbors == F) {boundaries[row, col] = 0}
        }
      }
    }
  }

  boundaries <- raster(boundaries)
  raster::extent(boundaries) <- raster::extent(x)
  return(boundaries)

}

# define boundary with categorical data ----
#' @name categorical_boundary
#' @title Define the boundary elements of a RasterLayer with categorical data
#' @description Creates boundary element cells where patches of two categories meet.
#'
#' @param x A Rasterlayer object.
#' @param projection The EPSG code of the input layer CRS. Numeric.
#'
#' @return A RasterLayer object with cell values 1 for boundary elements and 0 for other cells
#' @examples
#' soil_boundaries <- categorical_boundary(soil_type_raster, projection = 4326)
#'
#' @author Amy Luo
#' @export
categorical_boundary <- function(x, projection) {
  x_poly <- terra::rast(x) %>%
    terra::as.polygons(.)
  mask <- terra::aggregate(x_poly) %>%
    terra::as.lines(.) %>%
    terra::rasterize(., terra::rast(x))
  fill <- terra::aggregate(x_poly) %>%
    terra::rasterize(., terra::rast(x), field = 0)

  x_boundary <- terra::as.lines(x_poly) %>%
    terra::rasterize(., terra::rast(x)) %>%
    terra::merge(., fill) %>%
    terra::mask(., mask, maskvalues = 1, updatevalue = 0) %>%
    raster::raster(.)

  raster::crs(x_boundary) <- paste('+init=epsg:', projection, sep = '')
  return(x_boundary)
}
