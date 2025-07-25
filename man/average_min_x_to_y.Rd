% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlap_statistics.R
\name{average_min_x_to_y}
\alias{average_min_x_to_y}
\title{Average minimum distance from x boundary elements to nearest y boundary element.}
\usage{
average_min_x_to_y(x, y, null_distrib)
}
\arguments{
\item{x}{A SpatRaster object with boundary elements.}

\item{y}{A SpatRaster object with boundary elements.}

\item{null_distrib}{A list of probability functions output from overlap_null_distrib().}
}
\value{
The average minimum distance and a p-value.
}
\description{
Statistical test for the average minimum distance between each boundary element in raster x
and the nearest boundary element in raster y. Uses Euclidean distance. The boundaries of
trait x depend on the boundaries of trait y.
}
\examples{
\donttest{
data(T.cristatus)
T.cristatus <- terra::rast(T.cristatus_matrix, crs = T.cristatus_crs)
terra::ext(T.cristatus) <- T.cristatus_ext

data(grassland)
grassland <- terra::rast(grassland_matrix, crs = grassland_crs)
terra::ext(grassland) <- grassland_ext

Tcrist_ovlp_null <- overlap_null_distrib(T.cristatus, grassland, rand_both = FALSE,
  x_cat = TRUE, n_iterations = 100, x_model = 'random_cluster')
Tcrist_boundaries <- define_boundary(T.cristatus, cat = TRUE)
grassland_boundaries <- define_boundary(grassland, cat = FALSE, threshold = 0.1)

average_min_x_to_y(Tcrist_boundaries, grassland_boundaries, Tcrist_ovlp_null)
}

}
\references{
Jacquez, G.M., Maruca,I S. & Fortin,M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
Fortin, M.-J., Drapeau, P. & Jacquez, G.M. (1996) Quantification of the Spatial Co-Occurrences of Ecological Boundaries. Oikos, 77, 51-60.
}
\author{
Amy Luo
}
