#' @name sobel_operator
#' @title Sobel-Feldman operator for edge detection
#' @description Uses a Sobel-Feldman operator (3x3 kernel) to detect internal edges in a SpatRaster object.
#' 
#' @param x A SpatRaster object.
#' @return A SpatRaster object with boundary values.
#' 
#' @examples
#' data(T.cristatus)
#' T.cristatus <- terra::rast(T.cristatus_matrix, crs = T.cristatus_crs)
#' terra::ext(T.cristatus) <- T.cristatus_ext
#' 
#' edges <- sobel_operator(T.cristatus)
#' terra::plot(edges)
#' 
#' @author Amy Luo
#' @export
sobel_operator <- function(x) {
  # filters
  gx <- matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3)
  gy <- matrix(c(1, 0, -1, 2, 0, -2, 1, 0, -1), nrow = 3)
  
  # empty raster with same resolution, crs, and extent as input raster
  out <- terra::rast(nrows = nrow(x), ncols = ncol(x), extent = terra::ext(x), crs = terra::crs(x))
  
  # transform values from input raster to boundary intensities
  for (i in terra::cells(x)) {
    neighbors <- terra::adjacent(out, i, 'queen', include = TRUE) %>%
      terra::values(x)[.]
    
    if (length(na.omit(neighbors)) == 9) {
      neighbors <- matrix(neighbors, 3, 3)
      gxa <- sum(neighbors*gx, na.rm = T)
      gya <- sum(neighbors*gy, na.rm = T)
      out[i] <- sqrt(gxa ^ 2 + gya ^ 2)
    } else {
      out[i] = NA
    }
  }

  return(out)
}
