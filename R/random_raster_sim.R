#' @name random_raster_sim
#' @title Simulates random raster with boundary elements.
#' @description Internal function that simulates a random raster of the same extent and resolution as
#' the input raster using a neutral model.
#'
#' @param x A SpatRaster object.
#' @return A SpatRaster object with boundary elements.
#'
#' @author Amy Luo
#' @references
#' James, P. M. A., Fleming, R.A., & Fortin, M.-J. (2010) Identifying significant scale-specific spatial boundaries using wavelets and null models: Spruce budworm defoliation in Ontario, Canada as a case study. Landscape Ecology, 6, 873-887.
#' @export
random_raster_sim <- function (x) {
  values <- terra::values(x) %>%
    na.omit(.) %>%
    sample(.)
  cells_to_fill <- terra::cells(x) %>%
    terra::rowColFromCell(x, .) %>%
    t(.) %>%
    as.data.frame(.)
  
  x_sim <- terra::rast(nrow = terra::nrow(x), ncol = terra::ncol(x), crs = terra::crs(x), extent = terra::ext(x))
  index = 1
  for (i in cells_to_fill) {
    x_sim[i[1], i[2]] <- values[index]
    index = index + 1
  }

  return(x_sim)
  
}
